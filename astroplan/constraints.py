# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Specify and constraints to determine which targets are observable for
an observer.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from abc import ABCMeta, abstractmethod
import numpy as np
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import get_sun

DEFAULT_TIME_RESOLUTION = 0.5*u.hour

__all__ = ["AltitudeConstraint", "AirmassConstraint", "AtNight",
           "is_observable", "is_always_observable", "time_grid_from_range"]

@u.quantity_input(time_resolution=u.hour)
def time_grid_from_range(time_range, time_resolution=DEFAULT_TIME_RESOLUTION):
    """
    Get linearly-spaced sequence of times.

    Parameters
    ----------
    time_range : `~astropy.time.Time` of length=2
        Lower and upper bounds on time sequence

    time_resolution : `~astropy.units.quantity` (optional)
        Time-grid spacing

    Returns
    -------
    times : `~astropy.time.Time`
        Linearly-spaced sequence of times
    """
    return Time(np.arange(time_range[0].jd, time_range[1].jd,
                          time_resolution.to(u.day).value), format='jd')

class Constraint(object):
    """
    Abstract class for objects defining observational constraints.
    """
    __metaclass__ = ABCMeta

    def __call__(self, time_range, observer, targets):
        cons = self._compute_constraint(time_range, observer, targets)
        return cons

    def _compute_constraint(self, time_range, observer, targets):
        # Should be implemented on each subclass of Constraint
        raise NotImplementedError

    def _get_altaz(self, time_range, observer, targets,
                   time_resolution=DEFAULT_TIME_RESOLUTION,
                   force_zero_pressure=False):
        """
        Calculate alt/az for ``target`` at times linearly spaced between
        the two times in ``time_range`` with grid spacing ``time_resolution``
        for ``observer``.

        Cache the result on the ``observer`` object.

        Parameters
        ----------
        time_range : `~astropy.time.Time` with length=2
            Time range on which to compute these constraints

        targets : {list, `~astropy.coordinates.SkyCoord`, `~astroplan.FixedTarget`}
            Target or list of targets

        observer : `~astroplan.Observer`
            The observer who has constraints ``constraints``

        time_resolution : `~astropy.units.Quantity` (optional)
            Set the time resolution in calculations of the altitude/azimuth

        Returns
        -------
        altaz_dict : dict
            Dictionary containing two key-value pairs. (1) 'times' contains the
            times for the alt/az computations, (2) 'altaz' contains the
            corresponding alt/az coordinates at those times.
        """
        if not hasattr(observer, '_altaz_cache'):
            observer._altaz_cache = {}

        times = time_grid_from_range(time_range)

        # convert times, targets to tuple for hashing
        aakey = (tuple(times.jd), tuple(targets))

        if aakey not in observer._altaz_cache:
            if force_zero_pressure:
                observer_old_pressure = observer.pressure
                observer.pressure = 0

            altaz = observer.altaz(times, targets)
            observer._altaz_cache[aakey] = dict(times=times,
                                                altaz=altaz)
            if force_zero_pressure:
                observer.pressure = observer_old_pressure

        return observer._altaz_cache[aakey]

class AltitudeConstraint(Constraint):
    """
    Constrain the altitude of the target.

    Note: this will misbehave if you try to constrain negative altitudes, as
    the `~astropy.coordinates.AltAz` frame tends to mishandle negative
    altitudes.
    """
    def __init__(self, min=None, max=None):
        """
        Parameters
        ----------
        min : `~astropy.units.Quantity` or `None`
            Minimum altitude of the target. `None` indicates no limit.

        max : `~astropy.units.Quantity` or `None`
            Maximum altitude of the target. `None` indicates no limit.
        """
        if min is None:
            self.min = 0*u.deg
        else:
            self.min = min
        if max is None:
            self.max = 90*u.deg
        else:
            self.max = max

    def _compute_constraint(self, time_range, observer, targets,
                            time_resolution=DEFAULT_TIME_RESOLUTION):

        cached_altaz = self._get_altaz(time_range, observer, targets)
        altaz = cached_altaz['altaz']
        lowermask = self.min < altaz.alt
        uppermask = altaz.alt < self.max
        return lowermask & uppermask

class AirmassConstraint(AltitudeConstraint):
    """
    Constrain the airmass of a target.

    The airmass is approximated by the secant of the zenith angle.
    """
    def __init__(self, max=None, min=None):
        """
        Note: the ``max`` and ``min`` arguments appear in the order (max, min)
        in this initializer to support the common case for users who care about
        the upper limit on the airmass (``max``) and not the lower limit.

        Parameters
        ----------
        max : float or `None`
            Maximum airmass of the target. `None` indicates no limit.

        min : float or `None`
            Minimum airmass of the target. `None` indicates no limit.
        """
        self.min = min
        self.max = max

    def _compute_constraint(self, time_range, observer, targets):
        cached_altaz = self._get_altaz(time_range, observer, targets)
        altaz = cached_altaz['altaz']
        if self.min is None:
            mask = altaz.secz < self.max
        elif self.max is None:
            mask = self.min < altaz.secz
        elif self.min is not None and self.max is not None:
            mask = (self.min < altaz.secz) & (altaz.secz < self.max)
        return mask

class AtNight(Constraint):
    """
    Constrain the Sun to be below ``horizon``.
    """
    @u.quantity_input(horizon=u.deg)
    def __init__(self, max_solar_altitude=0*u.deg, force_pressure_zero=True):
        """
        Parameters
        ----------
        max_solar_altitude : `~astropy.units.Quantity`
            Define "night" as when the sun is below ``max_solar_altitude``.
            Default is zero degrees altitude.
        """
        self.max_solar_altitude = max_solar_altitude
        self.force_pressure_zero = force_pressure_zero

    @classmethod
    def twilight_civil(cls, **kwargs):
        return cls(max_solar_altitude=-6*u.deg, **kwargs)

    @classmethod
    def twilight_nautical(cls, **kwargs):
        return cls(max_solar_altitude=-12*u.deg, **kwargs)

    @classmethod
    def twilight_astronomical(cls, **kwargs):
        return cls(max_solar_altitude=-18*u.deg, **kwargs)

    def _get_solar_altitudes(self, time_range, observer, targets,
                             time_resolution=DEFAULT_TIME_RESOLUTION):
        if not hasattr(observer, '_altaz_cache'):
            observer._altaz_cache = {}

        times = Time(np.arange(time_range[0].jd, time_range[1].jd,
                               time_resolution.to(u.day).value), format='jd')

        aakey = (tuple(times.jd), 'sun')

        if aakey not in observer._altaz_cache:
            if self.force_pressure_zero:
                observer_old_pressure = observer.pressure
                observer.pressure = 0

            # Broadcast the solar altitudes for the number of targets:
            altaz = observer.altaz(times, get_sun(times))
            altitude = altaz.alt
            altitude.resize(1, len(altitude))
            altitude = altitude + np.zeros((len(targets), 1))

            observer._altaz_cache[aakey] = dict(times=times,
                                                altitude=altitude)
            if self.force_pressure_zero:
                observer.pressure = observer_old_pressure

        return observer._altaz_cache[aakey]

    def _compute_constraint(self, time_range, observer, targets):
        sun_altaz = self._get_solar_altitudes(time_range, observer, targets)
        solar_altitude = sun_altaz['altitude']
        mask = solar_altitude < self.max_solar_altitude
        return mask

def is_always_observable(constraints, time_range, targets, observer):
    """
    Are the ``targets`` always observable throughout ``time_range`` given
    constraints in ``constraints_list`` for ``observer``?

    Parameters
    ----------
    constraints : list or `~astroplan.constraints.Constraint`
        Observational constraint(s)

    time_range : `~astropy.time.Time` with length=2
        Time range on which to compute these constraints

    targets : {list, `~astropy.coordinates.SkyCoord`, `~astroplan.FixedTarget`}
        Target or list of targets

    observer : `~astroplan.Observer`
        The observer who has constraints ``constraints``

    Returns
    -------
    ever_observable : list
        List of booleans of same length as ``targets`` for whether or not each
        target is observable in the time range given the constraints.
    """
    if not hasattr(constraints, '__len__'):
        constraints = [constraints]

    contraint_arr = np.logical_and.reduce([constraint(time_range, observer,
                                                      targets)
                                          for constraint in constraints])
    return np.all(contraint_arr, axis=1)

def is_observable(constraints, time_range, targets, observer):
    """
    Are the ``targets`` observable during ``time_range`` given constraints in
    ``constraints_list`` for ``observer``?

    Parameters
    ----------
    constraints : list or `~astroplan.constraints.Constraint`
        Observational constraint(s)

    time_range : `~astropy.time.Time` with length=2
        Time range on which to compute these constraints

    targets : {list, `~astropy.coordinates.SkyCoord`, `~astroplan.FixedTarget`}
        Target or list of targets

    observer : `~astroplan.Observer`
        The observer who has constraints ``constraints``

    Returns
    -------
    ever_observable : list
        List of booleans of same length as ``targets`` for whether or not each
        target is ever observable in the time range given the constraints.
    """
    if not hasattr(constraints, '__len__'):
        constraints = [constraints]

    contraint_arr = np.logical_and.reduce([constraint(time_range, observer,
                                                      targets)
                                          for constraint in constraints])
    return np.any(contraint_arr, axis=1)
