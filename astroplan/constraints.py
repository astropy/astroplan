# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Specify and constraints to determine which targets are observable for
an observer.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# Standard library
from abc import ABCMeta, abstractmethod
import datetime

# Third-party
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import get_sun, Angle, SkyCoord
from astropy.coordinates.angle_utilities import angular_separation
from astropy import table
import numpy as np

# Package
from .moon import get_moon, moon_illumination
from .utils import time_grid_from_range
from .target import FixedTarget

__all__ = ["AltitudeConstraint", "AirmassConstraint", "AtNightConstraint",
           "is_observable", "is_always_observable", "time_grid_from_range",
           "SunSeparationConstraint", "MoonSeparationConstraint",
           "MoonIlluminationConstraint", "LocalTimeConstraint", "Constraint",
           "observability_table"]


def _get_altaz(times, observer, targets,
               force_zero_pressure=False):
    """
    Calculate alt/az for ``target`` at times linearly spaced between
    the two times in ``time_range`` with grid spacing ``time_resolution``
    for ``observer``.

    Cache the result on the ``observer`` object.

    Parameters
    ----------
    times : `~astropy.time.Time`
        Array of times on which to test the constraint

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

    # convert times, targets to tuple for hashing
    aakey = (tuple(times.jd), tuple(targets))

    if aakey not in observer._altaz_cache:
        try:
            if force_zero_pressure:
                observer_old_pressure = observer.pressure
                observer.pressure = 0

            altaz = observer.altaz(times, targets)
            observer._altaz_cache[aakey] = dict(times=times,
                                                altaz=altaz)
        finally:
            if force_zero_pressure:
                observer.pressure = observer_old_pressure

    return observer._altaz_cache[aakey]


@abstractmethod
class Constraint(object):
    """
    Abstract class for objects defining observational constraints.
    """
    __metaclass__ = ABCMeta

    def __call__(self, observer, targets, times=None,
                 time_range=None, time_grid_resolution=0.5*u.hour):

        if times is None and time_range is not None:
            times = time_grid_from_range(time_range,
                                         time_resolution=time_grid_resolution)
        elif not isinstance(times, Time):
            times = Time(times)

        if times.isscalar:
            times = Time([times])

        if hasattr(targets, '__len__'):
            targets = [FixedTarget(coord=target) if isinstance(target, SkyCoord)
                       else target for target in targets]
        else:
            if isinstance(targets, SkyCoord):
                targets = FixedTarget(coord=targets)

        cons = self.compute_constraint(times, observer, targets)
        return cons

    def compute_constraint(self, times, observer, targets):
        # Should be implemented on each subclass of Constraint
        raise NotImplementedError


class AltitudeConstraint(Constraint):
    """
    Constrain the altitude of the target.

    .. note::
        This will misbehave if you try to constrain negative altitudes, as
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

    def compute_constraint(self, times, observer, targets):

        cached_altaz = _get_altaz(times, observer, targets)
        altaz = cached_altaz['altaz']
        lowermask = self.min < altaz.alt
        uppermask = altaz.alt < self.max
        return lowermask & uppermask


class AirmassConstraint(AltitudeConstraint):
    """
    Constrain the airmass of a target.

    In the current implementation the airmass is approximated by the secant of
    the zenith angle.
    """
    def __init__(self, max=None, min=None):
        """
        .. note::
            The ``max`` and ``min`` arguments appear in the order (max, min)
            in this initializer to support the common case for users who care
            about the upper limit on the airmass (``max``) and not the lower
            limit.

        Parameters
        ----------
        max : float or `None`
            Maximum airmass of the target. `None` indicates no limit.

        min : float or `None`
            Minimum airmass of the target. `None` indicates no limit.

        Examples
        --------
        To create a constraint that requires the airmass be "better than 2",
        i.e. at a higher altitude than airmass=2::

            AirmassConstraint(2)
        """
        self.min = min
        self.max = max

    def compute_constraint(self, times, observer, targets):
        cached_altaz = _get_altaz(times, observer, targets)
        altaz = cached_altaz['altaz']
        if self.min is None and self.max is not None:
            mask = altaz.secz < self.max
        elif self.max is None and self.min is not None:
            mask = self.min < altaz.secz
        elif self.min is not None and self.max is not None:
            mask = (self.min < altaz.secz) & (altaz.secz < self.max)
        else:
            raise ValueError("No max and/or min specified in "
                             "AirmassConstraint.")
        return mask


class AtNightConstraint(Constraint):
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

        force_pressure_zero : bool (optional)
            Force the pressure to zero for solar altitude calculations. This
            avoids errors in the altitude of the Sun that can occur when the
            Sun is below the horizon and the corrections for atmospheric
            refraction return nonsense values. Default is `True`.
        """
        self.max_solar_altitude = max_solar_altitude
        self.force_pressure_zero = force_pressure_zero

    @classmethod
    def twilight_civil(cls, **kwargs):
        """
        Consider nighttime as time between civil twilights (-6 degrees).
        """
        return cls(max_solar_altitude=-6*u.deg, **kwargs)

    @classmethod
    def twilight_nautical(cls, **kwargs):
        """
        Consider nighttime as time between nautical twilights (-12 degrees).
        """
        return cls(max_solar_altitude=-12*u.deg, **kwargs)

    @classmethod
    def twilight_astronomical(cls, **kwargs):
        """
        Consider nighttime as time between astronomical twilights (-18 degrees).
        """
        return cls(max_solar_altitude=-18*u.deg, **kwargs)

    def _get_solar_altitudes(self, times, observer, targets):
        if not hasattr(observer, '_altaz_cache'):
            observer._altaz_cache = {}

        aakey = (tuple(times.jd), 'sun')

        if aakey not in observer._altaz_cache:
            try:
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
            finally:
                if self.force_pressure_zero:
                    observer.pressure = observer_old_pressure

        return observer._altaz_cache[aakey]

    def compute_constraint(self, times, observer, targets):
        sun_altaz = self._get_solar_altitudes(times, observer, targets)
        solar_altitude = sun_altaz['altitude']
        mask = solar_altitude < self.max_solar_altitude
        return mask


class SunSeparationConstraint(Constraint):
    """
    Constrain the distance between the Sun and some targets.
    """
    def __init__(self, min=None, max=None):
        """
        Parameters
        ----------
        min : `~astropy.units.Quantity` or `None` (optional)
            Minimum acceptable separation between Sun and target. `None`
            indicates no limit.
        max : `~astropy.units.Quantity` or `None` (optional)
            Minimum acceptable separation between Sun and target. `None`
            indicates no limit.
        """
        self.min = min
        self.max = max

    def compute_constraint(self, times, observer, targets):
        sun = get_sun(times)
        targets = [target.coord if hasattr(target, 'coord') else target
                   for target in targets]
        solar_separation = Angle([sun.separation(target) for target in targets])
        if self.min is None and self.max is not None:
            mask = self.max > solar_separation
        elif self.max is None and self.min is not None:
            mask = self.min < solar_separation
        elif self.min is not None and self.max is not None:
            mask = ((self.min < solar_separation) &
                    (solar_separation < self.max))
        else:
            raise ValueError("No max and/or min specified in "
                             "SunSeparationConstraint.")
        return mask


class MoonSeparationConstraint(Constraint):
    """
    Constrain the distance between the Earth's moon and some targets.
    """
    def __init__(self, min=None, max=None):
        """
        Parameters
        ----------
        min : `~astropy.units.Quantity` or `None` (optional)
            Minimum acceptable separation between moon and target. `None`
            indicates no limit.
        max : `~astropy.units.Quantity` or `None` (optional)
            Minimum acceptable separation between moon and target. `None`
            indicates no limit.
        """
        self.min = min
        self.max = max

    def compute_constraint(self, times, observer, targets):
        moon = get_moon(times, observer.location, observer.pressure)
        targets = [target.coord if hasattr(target, 'coord') else target
                   for target in targets]

        moon_separation = Angle([angular_separation(moon.spherical.lon,
                                                    moon.spherical.lat,
                                                    target.spherical.lon,
                                                    target.spherical.lat)
                                 for target in targets])
        # The line below should have worked, but needs a workaround.
        # TODO: once bug has been fixed, replace workaround with simpler version.
        # Relevant PR: https://github.com/astropy/astropy/issues/4033
#        moon_separation = Angle([moon.separation(target) for target in targets])
        if self.min is None and self.max is not None:
            mask = self.max > moon_separation
        elif self.max is None and self.min is not None:
            mask = self.min < moon_separation
        elif self.min is not None and self.max is not None:
            mask = ((self.min < moon_separation) &
                    (moon_separation < self.max))
        else:
            raise ValueError("No max and/or min specified in "
                             "MoonSeparationConstraint.")
        return mask


class MoonIlluminationConstraint(Constraint):
    """
    Constrain the fractional illumination of the Earth's moon.
    """
    def __init__(self, min=None, max=None):
        """
        Parameters
        ----------
        min : float or `None` (optional)
            Minimum acceptable fractional illumination. `None` indicates no
            limit.
        max : float or `None` (optional)
            Maximum acceptable fractional illumination. `None` indicates no
            limit.
        """
        self.min = min
        self.max = max

    def compute_constraint(self, times, observer, targets):
        illumination = np.array(moon_illumination(times,
                                                  observer.location))
        if self.min is None and self.max is not None:
            mask = self.max > illumination
        elif self.max is None and self.min is not None:
            mask = self.min < illumination
        elif self.min is not None and self.max is not None:
            mask = ((self.min < illumination) &
                    (illumination < self.max))
        else:
            raise ValueError("No max and/or min specified in "
                             "MoonSeparationConstraint.")
        return mask


class LocalTimeConstraint(Constraint):
    """
    Constrain the observable hours.
    """
    def __init__(self, min=None, max=None):
        """
        Parameters
        ----------
        min : `~datetime.time`
            Earliest local time. `None` indicates no limit.

        max : `~datetime.time`
            Latest local time. `None` indicates no limit.

        Examples
        --------
        Constrain the observations to targets that are observable between
        23:50 and 04:08 local time:

        >>> from astroplan import Observer
        >>> from astroplan.constraints import LocalTimeConstraint
        >>> import datetime as dt
        >>> subaru = Observer.at_site("Subaru", timezone="US/Hawaii")
        >>> constraint = LocalTimeConstraint(min=dt.time(23,50), max=dt.time(4,8)) # bound times between 23:50 and 04:08 local Hawaiian time
        """

        self.min = min
        self.max = max

        if self.min is None and self.max is None:
            raise ValueError("You must at least supply either a minimum or a maximum time.")

        if self.min is not None:
            if not isinstance(self.min, datetime.time):
                raise TypeError("Time limits must be specified as datetime.time objects.")

        if self.max is not None:
            if not isinstance(self.max, datetime.time):
                raise TypeError("Time limits must be specified as datetime.time objects.")

    def compute_constraint(self, times, observer, targets):

        timezone = None

        # get timezone from time objects, or from observer
        if self.min is not None:
            timezone = self.min.tzinfo

        elif self.max is not None:
            timezone = self.max.tzinfo

        if timezone is None:
            timezone = observer.timezone

        if self.min is not None:
            min_time = self.min
        else:
            min_time = self.min = datetime.time(0, 0, 0)

        if self.max is not None:
            max_time = self.max
        else:
            max_time = datetime.time(23, 59, 59)

        # If time limits occur on same day:
        if self.min < self.max:
            mask = [min_time < t.datetime.time() < max_time for t in times]

        # If time boundaries straddle midnight:
        else:
            mask = [(t.datetime.time() > min_time) or
                    (t.datetime.time() < max_time) for t in times]

        return mask


def is_always_observable(constraints, observer, targets, times=None,
                         time_range=None, time_grid_resolution=0.5*u.hour):
    """
    Are the ``targets`` always observable throughout ``time_range`` given
    constraints in ``constraints_list`` for ``observer``?

    Parameters
    ----------
    constraints : list or `~astroplan.constraints.Constraint`
        Observational constraint(s)

    observer : `~astroplan.Observer`
        The observer who has constraints ``constraints``

    targets : {list, `~astropy.coordinates.SkyCoord`, `~astroplan.FixedTarget`}
        Target or list of targets

    times : `~astropy.time.Time` (optional)
        Array of times on which to test the constraint

    time_range : `~astropy.time.Time` (optional)
        Lower and upper bounds on time sequence, with spacing
        ``time_resolution``. This will be passed as the first argument into
        `~astroplan.time_grid_from_range`.

    time_resolution : `~astropy.units.Quantity` (optional)
        If ``time_range`` is specified, determine whether constraints are met
        between test times in ``time_range`` by checking constraint at
        linearly-spaced times separated by ``time_resolution``. Default is 0.5
        hours.

    Returns
    -------
    ever_observable : list
        List of booleans of same length as ``targets`` for whether or not each
        target is observable in the time range given the constraints.
    """
    if not hasattr(constraints, '__len__'):
        constraints = [constraints]

    applied_constraints = [constraint(observer, targets, times=times,
                                      time_range=time_range,
                                      time_grid_resolution=time_grid_resolution)
                           for constraint in constraints]
    contraint_arr = np.logical_and.reduce(applied_constraints)
    return np.all(contraint_arr, axis=1)


def is_observable(constraints, observer, targets, times=None,
                  time_range=None, time_grid_resolution=0.5*u.hour):
    """
    Are the ``targets`` observable during ``time_range`` given constraints in
    ``constraints_list`` for ``observer``?

    Parameters
    ----------
    constraints : list or `~astroplan.constraints.Constraint`
        Observational constraint(s)

    observer : `~astroplan.Observer`
        The observer who has constraints ``constraints``

    targets : {list, `~astropy.coordinates.SkyCoord`, `~astroplan.FixedTarget`}
        Target or list of targets

    times : `~astropy.time.Time` (optional)
        Array of times on which to test the constraint

    time_range : `~astropy.time.Time` (optional)
        Lower and upper bounds on time sequence, with spacing
        ``time_resolution``. This will be passed as the first argument into
        `~astroplan.time_grid_from_range`.

    time_resolution : `~astropy.units.Quantity` (optional)
        If ``time_range`` is specified, determine whether constraints are met
        between test times in ``time_range`` by checking constraint at
        linearly-spaced times separated by ``time_resolution``. Default is 0.5
        hours.

    Returns
    -------
    ever_observable : list
        List of booleans of same length as ``targets`` for whether or not each
        target is ever observable in the time range given the constraints.
    """
    if not hasattr(constraints, '__len__'):
        constraints = [constraints]

    applied_constraints = [constraint(observer, targets, times=times,
                                      time_range=time_range,
                                      time_grid_resolution=time_grid_resolution)
                           for constraint in constraints]
    contraint_arr = np.logical_and.reduce(applied_constraints)
    return np.any(contraint_arr, axis=1)

def observability_table(constraints, observer, targets, times=None,
                        time_range=None, time_grid_resolution=0.5*u.hour):
    """
    Creates a table with information about observablity for all  the ``targets``
    over the requeisted ``time_range``, given the constraints in
    ``constraints_list`` for ``observer``.

    Parameters
    ----------
    constraints : list or `~astroplan.constraints.Constraint`
        Observational constraint(s)

    observer : `~astroplan.Observer`
        The observer who has constraints ``constraints``

    targets : {list, `~astropy.coordinates.SkyCoord`, `~astroplan.FixedTarget`}
        Target or list of targets

    times : `~astropy.time.Time` (optional)
        Array of times on which to test the constraint

    time_range : `~astropy.time.Time` (optional)
        Lower and upper bounds on time sequence, with spacing
        ``time_resolution``. This will be passed as the first argument into
        `~astroplan.time_grid_from_range`.

    time_resolution : `~astropy.units.Quantity` (optional)
        If ``time_range`` is specified, determine whether constraints are met
        between test times in ``time_range`` by checking constraint at
        linearly-spaced times separated by ``time_resolution``. Default is 0.5
        hours.

    Returns
    -------
    observability_table : `~astropy.table.Table`
        A Table containing the observability information for each of the
        ``targets``. The table contains four columns with information about the
        target and it's observability: ``'target name'``, ``'ever observable'``,
        ``'always observable'``, and ``'fraction of time observable'``.  It also
        contains metadata entries ``'times'`` (with an array of all the times),
        ``'observer'`` (the `~astroplan.Observer` object), and ``'constraints'``
        (containing the supplied ``constraints``).
    """
    if not hasattr(constraints, '__len__'):
        constraints = [constraints]

    applied_constraints = [constraint(observer, targets, times=times,
                                      time_range=time_range,
                                      time_grid_resolution=time_grid_resolution)
                           for constraint in constraints]
    contraint_arr = np.logical_and.reduce(applied_constraints)

    colnames = ['target name', 'ever observable', 'always observable',
                'fraction of time observable']

    target_names = [target.name for target in targets]
    ever_obs = np.any(contraint_arr, axis=1)
    always_obs = np.all(contraint_arr, axis=1)
    frac_obs = np.sum(contraint_arr, axis=1) / contraint_arr.shape[1]

    tab = table.Table(names=colnames, data=[target_names, ever_obs, always_obs, frac_obs])

    if times is None and time_range is not None:
        times = time_grid_from_range(time_range,
                                     time_resolution=time_grid_resolution)

    tab.meta['times'] = times.datetime
    tab.meta['observer'] = observer
    tab.meta['constraints'] = constraints

    return tab
