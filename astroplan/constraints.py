from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from abc import ABCMeta, abstractmethod
import numpy as np
from astropy.time import Time
import astropy.units as u

DEFAULT_TIME_RESOLUTION = 0.5*u.hour

class Constraint(object):
    """
    Abstract class for objects defining observational constraints.
    """
    __metaclass__ = ABCMeta

    def __call__(self, time_range, observer, targets):
        cons = self._compute_constraint(time_range, observer, targets)
        return cons

    def _compute_constraint(self, time_range, observer, targets):
        raise NotImplementedError

    def _get_altaz(self, time_range, observer, targets,
                   time_resolution=DEFAULT_TIME_RESOLUTION):
        """
        Calculate alt/az for ``target`` at times linearly spaced between
        the two times in ``time_range`` with grid spacing ``time_resolution``
        for ``observer``.

        Cache the result on the ``observer`` object.
        """
        if not hasattr(observer, '_altaz_cache'):
            observer._altaz_cache = {}

        times = Time(np.arange(time_range[0].jd, time_range[1].jd,
                               time_resolution.to(u.day).value), format='jd')

        # convert targets to tuple for hashing
        aakey = (times, tuple(targets))

        if aakey not in observer._altaz_cache:
            observer._altaz_cache[aakey] = dict(times=times,
                                                altaz=observer.altaz(times,
                                                                     targets))
        return observer._altaz_cache[aakey]

class AltitudeConstraint(Constraint):
    """
    Constrain the altitude of the target.

    Note: this will misbehave if you try to constrain negative altitudes, as
    the `~astropy.coordinates.AltAz` frame tends to.
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
                            time_resolution=0.5*u.hour):

        cached_altaz = self._get_altaz(time_range, observer, targets)
        altaz = cached_altaz['altaz']
        lowermsk = self.min < altaz.alt
        uppermsk = altaz.alt < self.max
        return lowermsk & uppermsk

class AirmassConstraint(AltitudeConstraint):
    """
    Constrain the airmass of a target.

    The airmass is approximated by the secant of the zenith angle.
    """
    def __init__(self, max=None, min=None):
        """
        Note: the ``max`` and ``min`` arguments appear in the order (max, min)
        in this initializer to support the common case for users who care about
        the upper limit on the airmass (max) and not the lower limit.

        Parameters
        ----------
        max : float or `None`
            Maximum airmass of the target. `None` indicates no limit.

        min : float or `None`
            Minimum airmass of the target. `None` indicates no limit.
        """
        if min is None:
            self.min = 0
        else:
            self.min = min
        if max is None:
            self.max = 100
        else:
            self.max = max

    def _compute_constraint(self, time_range, observer, targets):
        cached_altaz = self._get_altaz(time_range, observer, targets)
        altaz = cached_altaz['altaz']
        lowermsk = self.min < altaz.secz
        uppermsk = altaz.secz < self.max
        return lowermsk & uppermsk

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
