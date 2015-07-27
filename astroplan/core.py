# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from astropy.coordinates import (EarthLocation, SkyCoord, AltAz, get_sun,
                                 Angle, Latitude, Longitude, UnitSphericalRepresentation)

import astropy.units as u
import datetime
from astropy.time import Time
import pytz
import numpy as np
################################################################################
# TODO: Temporary solution to IERS tables problems
from astropy.utils.data import download_file
from astropy.utils import iers

iers.IERS.iers_table = iers.IERS_A.open(download_file(iers.IERS_A_URL,
                                                      cache=True))
################################################################################


from astropy.extern.six import string_types
from .exceptions import TargetNeverUpWarning, TargetAlwaysUpWarning
from .sites import get_site
from .moon import get_moon, moon_illumination, moon_phase_angle
import warnings

from abc import ABCMeta, abstractmethod

import numpy as np

__all__ = ["Observer", "Target", "FixedTarget", "NonFixedTarget",
           "Constraint", "TimeWindow", "AltitudeRange",
           "AboveAirmass", "Observation"]

#__doctest_requires__ = {'*': ['scipy.integrate']}

def _generate_24hr_grid(t0, start, end, N, for_deriv=False):
    """
    Generate a nearly linearly spaced grid of time durations.

    The midpoints of these grid points will span times from ``t0``+``start``
    to ``t0``+``end``, including the end points, which is useful when taking
    numerical derivatives.

    Parameters
    ----------
    t0 : `~astropy.time.Time`
        Time queried for, grid will be built from or up to this time.

    start : float
        Number of days before/after ``t0`` to start the grid.

    end : float
        Number of days before/after ``t0`` to end the grid.

    N : int
        Number of grid points to generate

    for_deriv : bool
        Generate time series for taking numerical derivative (modify
        bounds)?

    Returns
    -------
    `~astropy.time.Time`
    """

    if for_deriv:
        time_grid = np.concatenate([[start - 1/(N-1)],
                                    np.linspace(start, end, N)[1:-1],
                                    [end + 1/(N-1)]])*u.day
    else:
        time_grid = np.linspace(start, end, N)*u.day

    return t0 + time_grid

def transform_target_list_to_altaz(times, targets, location):
    """
    Workaround for transforming a list of coordinates ``targets`` to
    altitudes and azimuths.

    Parameters
    ----------
    times : `~astropy.time.Time` or list of `~astropy.time.Time` objects
        Time of observation

    targets : `~astropy.coordinates.SkyCoord` or list of `~astropy.coordinates.SkyCoord` objects
        List of target coordinates

    location : `~astropy.coordinates.EarthLocation`
        Location of observer

    Returns
    -------
    altitudes : list
        List of altitudes for each target, at each time
    """
    if times.isscalar:
        times = Time([times])

    if not isinstance(targets, list) and targets.isscalar:
        targets = [targets]

    repeated_times = np.tile(times, len(targets))
    repeated_targets = np.repeat(targets, len(times))
    target_SkyCoord = SkyCoord(SkyCoord(repeated_targets).data.represent_as(
                               UnitSphericalRepresentation),
                               representation=UnitSphericalRepresentation)

    transformed_coord = target_SkyCoord.transform_to(AltAz(location=location,
                                                           obstime=repeated_times))
    return transformed_coord

def _target_is_vector(target):
    if hasattr(target, '__iter__'):
        return True
    else:
        return False

def list_FixedTarget_to_SkyCoord(list_of_FixedTargets):
    """
    Convert a list of `~astroplan.core.FixedTarget` objects to a vector
    `~astropy.coordinates.SkyCoord` object.

    Parameters
    ----------
    list_of_FixedTargets : list
        `~astroplan.core.FixedTarget` objects

    Returns
    -------
    sc : `~astropy.coordinates.SkyCoord`
    """
    coord_list = [target.coord for target in list_of_FixedTargets]
    sc = SkyCoord(SkyCoord(coord_list).data.represent_as(
                  UnitSphericalRepresentation),
                  representation=UnitSphericalRepresentation)
    return sc

class Observer(object):
    """
    A container class for information about an observer's location and
    environment.

    TODO: write this docstring
    """
    @u.quantity_input(elevation=u.m)
    def __init__(self, location=None, timezone='UTC', name=None, latitude=None,
                 longitude=None, elevation=0*u.m, pressure=None,
                 relative_humidity=None, temperature=None, description=None):
        """
        Parameters
        ----------
        name : str
            A short name for the telescope, observatory or location.

        location : `~astropy.coordinates.EarthLocation`
            The location (latitude, longitude, elevation) of the observatory.

        longitude : float, str, `~astropy.units.Quantity` (optional)
            The longitude of the observing location. Should be valid input for
            initializing a `~astropy.coordinates.Longitude` object.

        latitude : float, str, `~astropy.units.Quantity` (optional)
            The latitude of the observing location. Should be valid input for
            initializing a `~astropy.coordinates.Latitude` object.

        elevation : `~astropy.units.Quantity` (optional), default = 0 meters
            The elevation of the observing location, with respect to sea
            level. Defaults to sea level.

        pressure : `~astropy.units.Quantity` (optional)
            The ambient pressure. Defaults to zero (i.e. no atmosphere).

        relative_humidity : float (optional)
            The ambient relative humidity.

        temperature : `~astropy.units.Quantity` (optional)
            The ambient temperature.

        timezone : str or `datetime.tzinfo` (optional)
            The local timezone to assume. If a string, it will be passed through
            `pytz.timezone()` to produce the timezone object.

        description : str (optional)
            A short description of the telescope, observatory or observing
            location.
        """

        self.name = name
        self.pressure = pressure
        self.temperature = temperature
        self.relative_humidity = relative_humidity

        # If lat/long given instead of EarthLocation, convert them
        # to EarthLocation
        if location is None and (latitude is not None and
                                 longitude is not None):
            self.location = EarthLocation.from_geodetic(longitude, latitude,
                                                        elevation)

        elif isinstance(location, EarthLocation):
            self.location = location

        else:
            raise TypeError('Observatory location must be specified with '
                            'either (1) an instance of '
                            'astropy.coordinates.EarthLocation or (2) '
                            'latitude and longitude in degrees as '
                            'accepted by astropy.coordinates.Latitude and '
                            'astropy.coordinates.Latitude.')

        # Accept various timezone inputs, default to UTC
        if isinstance(timezone, datetime.tzinfo):
            self.timezone = timezone
        elif isinstance(timezone, string_types):
            self.timezone = pytz.timezone(timezone)
        else:
            raise TypeError('timezone keyword should be a string, or an '
                            'instance of datetime.tzinfo')

    @classmethod
    def at_site(cls, site_name, **kwargs):
        """
        Initialize an `~astroplan.core.Observer` object with a site name.

        Extra keyword arguments are passed to the `~astroplan.core.Observer`
        constructor (see `~astroplan.core.Observer` for available keyword
        arguments).

        Parameters
        ----------
        site_name : str
            Observatory name, must be resolvable with
            `~astroplan.sites.get_site`.

        Returns
        -------
        `~astroplan.core.Observer`
            Observer object.

        Examples
        --------
        Initialize an observer at Kitt Peak National Observatory:

        >>> from astroplan import Observer
        >>> import astropy.units as u
        >>> kpno_generic = Observer.at_site('kpno')
        >>> kpno_today = Observer.at_site('kpno', pressure=1*u.bar, temperature=0*u.deg_C)
        """
        name = kwargs.pop('name', site_name)
        if 'location' in kwargs:
            raise ValueError("Location kwarg should not be used if "
                             "initializing an Observer with Observer.at_site()")
        return cls(location=get_site(site_name), name=name, **kwargs)

    def astropy_time_to_datetime(self, astropy_time):
        """
        Convert the `~astropy.time.Time` object ``astropy_time`` to a
        localized `~datetime.datetime` object.

        Timezones localized with `~pytz`.

        Parameters
        ----------
        astropy_time : `~astropy.time.Time`
            Scalar or list-like.

        Returns
        -------
        `~datetime.datetime`
            Localized datetime, where the timezone of the datetime is
            set by the ``timezone`` keyword argument of the
            `~astroplan.Observer` constructor.
        """

        if not astropy_time.isscalar:
            return [self.astropy_time_to_datetime(t) for t in astropy_time]

        # Convert astropy.time.Time to a UTC localized datetime (aware)
        utc_datetime = pytz.utc.localize(astropy_time.utc.datetime)

        # Convert UTC to local timezone
        return self.timezone.normalize(utc_datetime)

    def datetime_to_astropy_time(self, date_time):
        """
        Convert the `~datetime.datetime` object ``date_time`` to a
        `~astropy.time.Time` object.

        Timezones localized with `~pytz`. If the ``date_time`` is naive, the
        implied timezone is the ``timezone`` structure of ``self``.

        Parameters
        ----------
        date_time : `~datetime.datetime` or list-like

        Returns
        -------
        `~astropy.time.Time`
            Astropy time object (no timezone information preserved).
        """

        if hasattr(date_time, '__iter__'):
            return Time([self.datetime_to_astropy_time(t) for t in date_time])

        # For timezone-naive datetimes, assign local timezone
        if date_time.tzinfo is None:
            date_time = self.timezone.localize(date_time)

        return Time(date_time, location=self.location)

    def altaz(self, time, target=None, obswl=None):
        """
        Get an `~astropy.coordinates.AltAz` frame or coordinate.

        If ``target`` is None, generates an altitude/azimuth frame. Otherwise,
        calculates the transformation to that frame for the requested ``target``.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            The time at which the observation is taking place. Will be used as
            the ``obstime`` attribute in the resulting frame or coordinate. This
            will be passed in as the first argument to the `~astropy.time.Time`
            initializer, so it can be anything that `~astropy.time.Time` will
            accept (including a `~astropy.time.Time` object)

        target : `~astroplan.FixedTarget`, `~astropy.coordinates.SkyCoord`, or list; defaults to `None` (optional)
            Celestial object(s) of interest. If ``target`` is `None`, return the
            `~astropy.coordinates.AltAz` frame without coordinates.

        obswl : `~astropy.units.Quantity` (optional)
            Wavelength of the observation used in the calculation.

        Returns
        -------
        `~astropy.coordinates.AltAz`
            If ``target`` is `None`, returns `~astropy.coordinates.AltAz` frame.
            If ``target`` is not `None`, returns the ``target`` transformed to
            the `~astropy.coordinates.AltAz` frame.
        """
        if not isinstance(time, Time):
            time = Time(time)

        altaz_frame = AltAz(location=self.location, obstime=time,
                            pressure=self.pressure, obswl=obswl,
                            temperature=self.temperature,
                            relative_humidity=self.relative_humidity)
        if target is None:
            # Return just the frame
            return altaz_frame
        else:
            # If target is a list of targets:
            if _target_is_vector(target):
                get_coord = lambda x: x.coord if hasattr(x, 'coord') else x
                transformed_coords = transform_target_list_to_altaz(time,
                                                                    list(map(get_coord,
                                                                        target)),
                                                                    self.location)
                n_targets = len(target)
                new_shape = (n_targets, int(len(transformed_coords)/n_targets))

                for comp in transformed_coords.data.components:
                    getattr(transformed_coords.data, comp).resize(new_shape)
                return transformed_coords

            # If single target is a FixedTarget or a SkyCoord:
            if hasattr(target, 'coord'):
                coordinate = target.coord
            else:
                coordinate = target
            return coordinate.transform_to(altaz_frame)

    def parallactic_angle(self, time, target):
        '''
        Calculate the parallactic angle.

        Parameters
        ----------
        time : `~astropy.time.Time`
            Observation time.

        target : `~astroplan.FixedTarget` or `~astropy.coordinates.SkyCoord`
            Celestial object of interest.

        Returns
        -------
        `~astropy.coordinates.Angle`
            Parallactic angle.

        Notes
        -----
        The parallactic angle is the angle between the great circle that
        intersects a celestial object and the zenith, and the object's hour
        circle [1]_.

        .. [1] https://en.wikipedia.org/wiki/Parallactic_angle

        '''
        if not isinstance(time, Time):
            time = Time(time)

        if _target_is_vector(target):
            get_coord = lambda x: x.coord if hasattr(x, 'coord') else x
            coordinate = SkyCoord(list(map(get_coord, target)))
        else:
            if hasattr(target, 'coord'):
                coordinate = target.coord
            else:
                coordinate = target

        # Eqn (14.1) of Meeus' Astronomical Algorithms
        LST = time.sidereal_time('mean', longitude=self.location.longitude)
        H = (LST - coordinate.ra).radian
        q = np.arctan(np.sin(H) /
                      (np.tan(self.location.latitude.radian)*
                       np.cos(coordinate.dec.radian) -
                       np.sin(coordinate.dec.radian)*np.cos(H)))*u.rad

        return Angle(q)

    # Sun-related methods.
    @u.quantity_input(horizon=u.deg)
    def _horiz_cross(self, t, alt, rise_set, horizon=0*u.degree):
        """
        Find time ``t`` when values in array ``a`` go from
        negative to positive or positive to negative (exclude endpoints)

        ``return_limits`` will return nearest times to zero-crossing.

        Parameters
        ----------
        t : `~astropy.time.Time`
            Grid of times
        alt : `~astropy.units.Quantity`
            Grid of altitudes
        rise_set : {"rising",  "setting"}
            Calculate either rising or setting across the horizon
        horizon : float
            Number of degrees above/below actual horizon to use
            for calculating rise/set times (i.e.,
            -6 deg horizon = civil twilight, etc.)
        Returns
        -------
        Returns the lower and upper limits on the time and altitudes
        of the horizon crossing.
        """
        alt = Latitude(alt)

        if len(np.shape(alt)) == 1:
            alt = alt[np.newaxis, :]
        n_targets = alt.shape[0]

        if rise_set == 'rising':
            # Find index where altitude goes from below to above horizon
            condition = (alt[:, :-1] < horizon) * (alt[:, 1:] > horizon)
        elif rise_set == 'setting':
            # Find index where altitude goes from above to below horizon
            condition = (alt[:, :-1] > horizon) * (alt[:, 1:] < horizon)

        target_inds, time_inds = np.nonzero(condition)

        if np.count_nonzero(condition) < n_targets:
            target_inds, _ = np.nonzero(condition)
            noncrossing_target_ind = np.setdiff1d(np.arange(n_targets),
                                                  target_inds,
                                                  assume_unique=True)#[0]

            warnmsg = ('Target(s) index {} does not cross horizon={} within '
                       '24 hours'.format(noncrossing_target_ind, horizon))

            if (alt[noncrossing_target_ind, :] > horizon).all():
                warnings.warn(warnmsg, TargetAlwaysUpWarning)
            else:
                warnings.warn(warnmsg, TargetNeverUpWarning)

            # Fill in missing time with np.nan
            target_inds = np.insert(target_inds, noncrossing_target_ind,
                                    noncrossing_target_ind)
            time_inds = np.insert(time_inds.astype(float),
                                  noncrossing_target_ind,
                                  np.nan)
        elif np.count_nonzero(condition) > n_targets:
            dup_target_inds = list(set([target_ind for target_ind in target_inds
                                        if np.sum(target_inds == target_ind) > 1]))
            old_target_inds = np.copy(target_inds)
            old_time_inds = np.copy(time_inds)

            time_inds = []
            target_inds = []
            for tgt, tm in zip(old_target_inds, old_time_inds):
                if tgt not in target_inds:
                    time_inds.append(tm)
                    target_inds.append(tgt)
            target_inds = np.array(target_inds)
            time_inds = np.array(time_inds)

        times = [t[i:i+2] if not np.isnan(i) else i for i in time_inds]
        altitudes = [alt[i, j:j+2] if not np.isnan(j) else j
                     for i, j in zip(target_inds, time_inds)]

        return times, altitudes

    @u.quantity_input(horizon=u.deg)
    def _two_point_interp(self, times, altitudes, horizon=0*u.deg):
        """
        Do linear interpolation between two ``altitudes`` at
        two ``times`` to determine the time where the altitude
        goes through zero.

        Parameters
        ----------
        times : `~astropy.time.Time`
            Two times for linear interpolation between

        altitudes : array of `~astropy.units.Quantity`
            Two altitudes for linear interpolation between

        horizon : `~astropy.units.Quantity`
            Solve for the time when the altitude is equal to
            reference_alt.

        Returns
        -------
        t : `~astropy.time.Time`
            Time when target crosses the horizon

        """
        if not isinstance(times, Time) and np.any(np.isnan(times)):
            return np.nan
        else:
            slope = (altitudes[1] - altitudes[0])/(times[1].jd - times[0].jd)
            return Time(times[1].jd - ((altitudes[1] - horizon)/slope).value,
                        format='jd')

    def _altitude_trig(self, LST, target):
        """
        Calculate the altitude of ``target`` at local sidereal times ``LST``.

        This method provides a factor of ~3 speed up over calling `altaz`, and
        inherently does *not* take the atmosphere into account.

        Parameters
        ----------
        LST : `~astropy.time.Time`
            Local sidereal times (array)

        target : {`~astropy.coordinates.SkyCoord`, `FixedTarget`} or similar
            Target celestial object's coordinates.

        Returns
        -------
        alt : `~astropy.unit.Quantity`
            Array of altitudes
        """
        alt = np.arcsin(np.sin(self.location.latitude.radian) *
                        np.sin(target.dec) +
                        np.cos(self.location.latitude.radian) *
                        np.cos(target.dec) *
                        np.cos(LST.radian - target.ra.radian))
        return alt

    def _calc_riseset(self, time, target, prev_next, rise_set, horizon, N=150):
        """
        Time at next rise/set of ``target``.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observation. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object)

        target : `~astropy.coordinates.SkyCoord`
            Position of target or multiple positions of that target
            at multiple times (if target moves, like the Sun)

        prev_next : str - either 'previous' or 'next'
            Test next rise/set or previous rise/set

        rise_set : str - either 'rising' or 'setting'
            Compute prev/next rise or prev/next set

        location : `~astropy.coordinates.EarthLocation`
            Location of observer

        horizon : `~astropy.units.Quantity`
            Degrees above/below actual horizon to use
            for calculating rise/set times (i.e.,
            -6 deg horizon = civil twilight, etc.)

        N : int
            Number of altitudes to compute when searching for
            rise or set.

        Returns
        -------
        ret1 : `~astropy.time.Time`
            Time of rise/set
        """

        if not isinstance(time, Time):
            time = Time(time)

        target_is_vector = _target_is_vector(target)

        if prev_next == 'next':
            times = _generate_24hr_grid(time, 0, 1, N)
        else:
            times = _generate_24hr_grid(time, -1, 0, N)

        altaz = self.altaz(times, target)
        if target_is_vector:
            altitudes = [aa.alt for aa in altaz]
        else:
            altitudes = altaz.alt

        time_limits, altitude_limits = self._horiz_cross(times, altitudes, rise_set,
                                                    horizon)
        if not target_is_vector:
            return self._two_point_interp(time_limits[0], altitude_limits[0],
                                          horizon=horizon)
        else:
            return Time([self._two_point_interp(time_limit, altitude_limit,
                                                horizon=horizon)
                         for time_limit, altitude_limit in
                         zip(time_limits, altitude_limits)])

    def _calc_transit(self, time, target, prev_next, antitransit=False, N=150):
        """
        Time at next transit of the meridian of `target`.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observation. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object)

        target : `~astropy.coordinates.SkyCoord`
            Position of target or multiple positions of that target
            at multiple times (if target moves, like the Sun)

        prev_next : str - either 'previous' or 'next'
            Test next rise/set or previous rise/set

        antitransit : bool
            Toggle compute antitransit (below horizon, equivalent to midnight
            for the Sun)

        location : `~astropy.coordinates.EarthLocation`
            Location of observer

        N : int
            Number of altitudes to compute when searching for
            rise or set.

        Returns
        -------
        ret1 : `~astropy.time.Time`
            Time of transit/antitransit
        """
        if not isinstance(time, Time):
            time = Time(time)

        target_is_vector = _target_is_vector(target)

        if prev_next == 'next':
            times = _generate_24hr_grid(time, 0, 1, N, for_deriv=True)
        else:
            times = _generate_24hr_grid(time, -1, 0, N, for_deriv=True)

        # The derivative of the altitude with respect to time is increasing
        # from negative to positive values at the anti-transit of the meridian
        if antitransit:
            rise_set = 'rising'
        else:
            rise_set = 'setting'

        altaz = self.altaz(times, target)
        if target_is_vector:
            #d_altitudes = [each_altaz.alt.diff() for each_altaz in altaz]
            d_altitudes = [each_alt.diff() for each_alt in altaz.alt]
        else:
            altitudes = altaz.alt
            d_altitudes = altitudes.diff()

        print('alts:', [d_altitudes])
        dt = Time((times.jd[1:] + times.jd[:-1])/2, format='jd')

        horizon = 0*u.degree # Find when derivative passes through zero
        time_limits, altitude_limits = self._horiz_cross(dt, d_altitudes,
                                                         rise_set, horizon)
        if not target_is_vector:
            return self._two_point_interp(time_limits[0], altitude_limits[0],
                                          horizon=horizon)
        else:
            return Time([self._two_point_interp(time_limit, altitude_limit,
                                                horizon=horizon)
                         for time_limit, altitude_limit in
                         zip(time_limits, altitude_limits)])

    @u.quantity_input(horizon=u.deg)
    def target_rise_time(self, time, target, which='nearest', horizon=0*u.degree):
        """
        Calculate rise time.

        Compute time of the next/previous/nearest rise of the ``target``
        object, where "rise" is defined as the time when the ``target``
        transitions from altitudes below the ``horizon`` to above the
        ``horizon``.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observation. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object)

        target : coordinate object (i.e. `~astropy.coordinates.SkyCoord`, `~astroplan.FixedTarget`)
            Target celestial object

        which : {'next', 'previous', 'nearest'}
            Choose which sunrise relative to the present ``time`` would you
            like to calculate

        horizon : `~astropy.units.Quantity` (optional), default = zero degrees
            Degrees above/below actual horizon to use
            for calculating rise/set times (i.e.,
            -6 deg horizon = civil twilight, etc.)

        Returns
        -------
        `~astropy.time.Time`
            Rise time of target
        """
        if not isinstance(time, Time):
            time = Time(time)

        if which == 'next' or which == 'nearest':
            next_rise = self._calc_riseset(time, target, 'next', 'rising',
                                           horizon)
            if which == 'next':
                return next_rise

        if which == 'previous' or which == 'nearest':
            previous_rise = self._calc_riseset(time, target, 'previous',
                                               'rising', horizon)
            if which == 'previous':
                return previous_rise

        if which == 'nearest':
            if _target_is_vector(target):
                return_times = []
                for next_r, prev_r in zip(next_rise, previous_rise):
                    if abs(time - prev_r) < abs(time - next_r):
                        return_times.append(prev_r)
                    else:
                        return_times.append(next_r)
                return Time(return_times)
            else:
                if abs(time - previous_rise) < abs(time - next_rise):
                    return previous_rise
                else:
                    return next_rise

        raise ValueError('"which" kwarg must be "next", "previous" or '
                         '"nearest".')

    @u.quantity_input(horizon=u.deg)
    def target_set_time(self, time, target, which='nearest', horizon=0*u.degree):
        """
        Calculate set time.

        Compute time of the next/previous/nearest set of ``target``, where
        "set" is defined as when the ``target`` transitions from altitudes
        above ``horizon`` to below ``horizon``.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observation. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object)

        target : coordinate object (i.e. `~astropy.coordinates.SkyCoord`, `~astroplan.FixedTarget`)
            Target celestial object

        which : {'next', 'previous', 'nearest'}
            Choose which sunset relative to the present ``time`` would you
            like to calculate

        horizon : `~astropy.units.Quantity` (optional), default = zero degrees
            Degrees above/below actual horizon to use
            for calculating rise/set times (i.e.,
            -6 deg horizon = civil twilight, etc.)

        Returns
        -------
        `~astropy.time.Time`
            Set time of target.
        """
        if not isinstance(time, Time):
            time = Time(time)

        if which == 'next' or which == 'nearest':
            next_set = self._calc_riseset(time, target, 'next', 'setting',
                                          horizon)
            if which == 'next':
                return next_set

        if which == 'previous' or which == 'nearest':
            previous_set = self._calc_riseset(time, target, 'previous',
                                              'setting', horizon)
            if which == 'previous':
                return previous_set

        if which == 'nearest':
            if _target_is_vector(target):
                return_times = []
                for next_s, prev_s in zip(next_set, previous_set):
                    if abs(time - prev_s) < abs(time - next_s):
                        return_times.append(prev_s)
                    else:
                        return_times.append(next_s)
                return Time(return_times)
            else:
                if abs(time - previous_set) < abs(time - next_set):
                    return previous_set
                else:
                    return next_set

        raise ValueError('"which" kwarg must be "next", "previous" or '
                         '"nearest".')

    def target_meridian_transit_time(self, time, target, which='nearest'):
        """
        Calculate time at the transit of the meridian.

        Compute time of the next/previous/nearest transit of the ``target``
        object.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observation. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object)

        target : coordinate object (i.e. `~astropy.coordinates.SkyCoord`, `~astroplan.FixedTarget`)
            Target celestial object

        which : {'next', 'previous', 'nearest'}
            Choose which sunrise relative to the present ``time`` would you
            like to calculate

        Returns
        -------
        `~astropy.time.Time`
            Transit time of target
        """
        if not isinstance(time, Time):
            time = Time(time)

        if which == 'next' or which == 'nearest':
            next_transit = self._calc_transit(time, target, 'next')
            if which == 'next':
                return next_transit

        if which == 'previous' or which == 'nearest':
            previous_transit = self._calc_transit(time, target, 'previous')
            if which == 'previous':
                return previous_transit

        if which == 'nearest':
            if _target_is_vector(target):
                return_times = []
                for next_t, prev_t in zip(next_transit, previous_transit):
                    print(prev_t, next_t,
                          abs(time - prev_t) < abs(time - next_t),
                          abs(time - prev_t), abs(time - next_t))
                    if abs(time - prev_t) < abs(time - next_t):
                        return_times.append(prev_t)
                    else:
                        return_times.append(next_t)
                return Time(return_times)
            else:
                print(previous_transit, next_transit,
                      abs(time - previous_transit) < abs(time - next_transit),
                      abs(time - previous_transit), abs(time - next_transit))
                if abs(time - previous_transit) < abs(time - next_transit):
                    return previous_transit
                else:
                    return next_transit
        raise ValueError('"which" kwarg must be "next", "previous" or '
                         '"nearest".')

    def target_meridian_antitransit_time(self, time, target, which='nearest'):
        """
        Calculate time at the antitransit of the meridian.

        Compute time of the next/previous/nearest antitransit of the ``target``
        object.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observation. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        target : coordinate object (i.e. `~astropy.coordinates.SkyCoord`, `~astroplan.FixedTarget`)
            Target celestial object

        which : {'next', 'previous', 'nearest'}
            Choose which sunrise relative to the present ``time`` would you
            like to calculate

        Returns
        -------
        `~astropy.time.Time`
            Antitransit time of target
        """
        if not isinstance(time, Time):
            time = Time(time)

        if which == 'next' or which == 'nearest':
            next_antitransit = self._calc_transit(time, target, 'next',
                                                  antitransit=True)
            if which == 'next':
                return next_antitransit

        if which == 'previous' or which == 'nearest':
            previous_antitransit = self._calc_transit(time, target, 'previous',
                                                      antitransit=True)
            if which == 'previous':
                return previous_antitransit

        if which == 'nearest':
            if _target_is_vector(target):
                return_times = []
                for next_at, prev_at in zip(next_antitransit,
                                          previous_antitransit):
                    if abs(time - prev_at) < abs(time - next_at):
                        return_times.append(prev_at)
                    else:
                        return_times.append(next_at)
                return Time(return_times)
            else:
                if (abs(time - previous_antitransit) <
                        abs(time - next_antitransit)):
                    return previous_antitransit
                else:
                    return next_antitransit

        raise ValueError('"which" kwarg must be "next", "previous" or '
                         '"nearest".')

    @u.quantity_input(horizon=u.deg)
    def sun_rise_time(self, time, which='nearest', horizon=0*u.degree):
        """
        Time of sunrise.

        Compute time of the next/previous/nearest sunrise, where
        sunrise is defined as when the Sun transitions from altitudes
        below ``horizon`` to above ``horizon``.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observation. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        which : {'next', 'previous', 'nearest'}
            Choose which sunrise relative to the present ``time`` would you
            like to calculate.

        horizon : `~astropy.units.Quantity` (optional), default = zero degrees
            Degrees above/below actual horizon to use
            for calculating rise/set times (i.e.,
            -6 deg horizon = civil twilight, etc.)

        Returns
        -------
        `~astropy.time.Time`
            Time of sunrise
        """
        return self.target_rise_time(time, get_sun(time), which, horizon)

    @u.quantity_input(horizon=u.deg)
    def sun_set_time(self, time, which='nearest', horizon=0*u.degree):
        """
        Time of sunset.

        Compute time of the next/previous/nearest sunset, where
        sunset is defined as when the Sun transitions from altitudes
        below ``horizon`` to above ``horizon``.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observation. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        which : {'next', 'previous', 'nearest'}
            Choose which sunset relative to the present ``time`` would you
            like to calculate

        horizon : `~astropy.units.Quantity` (optional), default = zero degrees
            Degrees above/below actual horizon to use
            for calculating rise/set times (i.e.,
            -6 deg horizon = civil twilight, etc.)

        Returns
        -------
        `~astropy.time.Time`
            Time of sunset
        """
        return self.target_set_time(time, get_sun(time), which, horizon)

    def noon(self, time, which='nearest'):
        """
        Time at solar noon.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observation. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        which : {'next', 'previous', 'nearest'}
            Choose which noon relative to the present ``time`` would you
            like to calculate

        Returns
        -------
        `~astropy.time.Time`
            Time at solar noon
        """
        return self.target_meridian_transit_time(time, get_sun(time), which)

    def midnight(self, time, which='nearest'):
        """
        Time at solar midnight.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observation. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        which : {'next', 'previous', 'nearest'}
            Choose which noon relative to the present ``time`` would you
            like to calculate

        Returns
        -------
        `~astropy.time.Time`
            Time at solar midnight
        """
        return self.target_meridian_antitransit_time(time, get_sun(time), which)

    # Twilight convenience functions

    def evening_astronomical(self, time, which='nearest'):
        """
        Time at evening astronomical (-18 degree) twilight.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observations. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        which : {'next', 'previous', 'nearest'}
            Choose which twilight relative to the present ``time`` would you
            like to calculate. Default is nearest.

        Returns
        -------
        `~astropy.time.Time`
            Time of twilight
        """
        return self.sun_set_time(time, which, horizon=-18*u.degree)

    def evening_nautical(self, time, which='nearest'):

        """
        Time at evening nautical (-12 degree) twilight.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observations. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        which : {'next', 'previous', 'nearest'}
            Choose which twilight relative to the present ``time`` would you
            like to calculate. Default is nearest.

        Returns
        -------
        `~astropy.time.Time`
            Time of twilight
        """
        return self.sun_set_time(time, which, horizon=-12*u.degree)

    def evening_civil(self, time, which='nearest'):
        """
        Time at evening civil (-6 degree) twilight.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observations. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        which : {'next', 'previous', 'nearest'}
            Choose which twilight relative to the present ``time`` would you
            like to calculate. Default is nearest.

        Returns
        -------
        `~astropy.time.Time`
            Time of twilight
        """
        return self.sun_set_time(time, which, horizon=-6*u.degree)

    def morning_astronomical(self, time, which='nearest'):
        """
        Time at morning astronomical (-18 degree) twilight.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observations. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        which : {'next', 'previous', 'nearest'}
            Choose which twilight relative to the present ``time`` would you
            like to calculate

        Returns
        -------
        `~astropy.time.Time`
            Time of twilight
        """
        return self.sun_rise_time(time, which, horizon=-18*u.degree)

    def morning_nautical(self, time, which='nearest'):
        """
        Time at morning nautical (-12 degree) twilight.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observations. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        which : {'next', 'previous', 'nearest'}
            Choose which twilight relative to the present ``time`` would you
            like to calculate. Default is nearest.

        Returns
        -------
        `~astropy.time.Time`
            Time of twilight
        """
        return self.sun_rise_time(time, which, horizon=-12*u.degree)

    def morning_civil(self, time, which='nearest'):
        """
        Time at morning civil (-6 degree) twilight.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observations. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        which : {'next', 'previous', 'nearest'}
            Choose which twilight relative to the present ``time`` would you
            like to calculate. Default is nearest.

        Returns
        -------
        `~astropy.time.Time`
            Time of sunset
        """
        return self.sun_rise_time(time, which, horizon=-6*u.degree)

    # Moon-related methods.

    def moon_rise_time(self, time, **kwargs):
        """
        Returns the local moonrise time.

        The default moonrise returned is the next one to occur.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)

        Keywords: str, optional
            previous
            next
        """
        raise NotImplementedError()

    def moon_set_time(self, time, **kwargs):
        """
        Returns the local moonset time.

        The default moonset returned is the next one to occur.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        Keywords: str, optional
            previous
            next
        """
        raise NotImplementedError()

    def moon_illumination(self, time):
        """
        Calculate the illuminated fraction of the moon

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        moon : `~astropy.coordinates.SkyCoord` or `None` (default)
            Position of the moon at time ``time``. If `None`, will calculate
            the position of the moon with `~astroplan.moon.get_moon`.

        sun : `~astropy.coordinates.SkyCoord` or `None` (default)
            Position of the sun at time ``time``. If `None`, will calculate
            the position of the Sun with `~astropy.coordinates.get_sun`.

        Returns
        -------
        float
            Fraction of lunar surface illuminated
        """
        if not isinstance(time, Time):
            time = Time(time)

        return moon_illumination(time, self.location)

    def moon_phase(self, time=None, moon=None, sun=None):
        """
        Calculate lunar orbital phase.

        For example, phase=0 is "new", phase=1 is "full".

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        moon : `~astropy.coordinates.SkyCoord` or `None` (default)
            Position of the moon at time ``time``. If `None`, will calculate
            the position of the moon with `~astroplan.moon.get_moon`.

        sun : `~astropy.coordinates.SkyCoord` or `None` (default)
            Position of the sun at time ``time``. If `None`, will calculate
            the position of the Sun with `~astropy.coordinates.get_sun`.
        """
        if time is not None and not isinstance(time, Time):
            time = Time(time)

        return moon_phase_angle(time, self.location)

    def moon_altaz(self, time):
        """
        Returns the position of the moon in alt/az.

        TODO: Currently `moon_altaz` uses PyEphem to calculate the position
        of the moon.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        Returns
        -------
        altaz : `~astropy.coordinates.SkyCoord`
            Position of the moon transformed to altitude and azimuth
        """
        if not isinstance(time, Time):
            time = Time(time)

        try:
            import ephem
        except ImportError:
            raise ImportError("The moon_altaz function currently requires "
                              "PyEphem to compute the position of the moon.")

        moon = ephem.Moon()
        obs = ephem.Observer()
        obs.lat = self.location.latitude.to(u.degree).to_string(sep=':')
        obs.lon = self.location.longitude.to(u.degree).to_string(sep=':')
        obs.elevation = self.location.height.to(u.m).value
        if self.pressure is not None:
            obs.pressure = self.pressure.to(u.bar).value*1000.0

        if time.isscalar:
            obs.date = time.datetime
            moon.compute(obs)
            moon_alt = float(moon.alt)
            moon_az = float(moon.az)
        else:
            moon_alt = []
            moon_az = []
            for t in time:
                obs.date = t.datetime
                moon.compute(obs)
                moon_alt.append(float(moon.alt))
                moon_az.append(float(moon.az))
        return SkyCoord(alt=moon_alt*u.rad, az=moon_az*u.rad,
                        frame=self.altaz(time))

    @u.quantity_input(horizon=u.deg)
    def target_is_up(self, time, target, horizon=0*u.degree, return_altaz=False):
        """
        Is ``target`` above ``horizon`` at this ``time``?

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observation. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object)

        target : coordinate object (i.e. `~astropy.coordinates.SkyCoord`, `~astroplan.FixedTarget`)
            Target celestial object

        horizon : `~astropy.units.Quantity` (optional), default = zero degrees
            Degrees above/below actual horizon to use
            for calculating rise/set times (i.e.,
            -6 deg horizon = civil twilight, etc.)

        return_altaz : bool (optional)
            Also return the '~astropy.coordinates.AltAz' coordinate.
        """
        if not isinstance(time, Time):
            time = Time(time)

        altaz = self.altaz(time, target)
        if _target_is_vector(target):
            observable = [alt > horizon for alt in altaz.alt]
        else:
            altitudes = altaz.alt
            observable = altitudes > horizon

        if not return_altaz:
            return observable
        else:
            return observable, altaz

    @u.quantity_input(horizon=u.deg)
    def is_night(self, time, horizon=0*u.deg, obswl=None):
        """
        Is the Sun below ``horizon`` at ``time``?

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observation. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object)

        horizon : `~astropy.units.Quantity` (optional), default = zero degrees
            Degrees above/below actual horizon to use
            for calculating day/night (i.e.,
            -6 deg horizon = civil twilight, etc.)

        obswl : `~astropy.units.Quantity` (optional)
            Wavelength of the observation used in the calculation

        Returns
        -------
        sun_below_horizon : bool
            `True` if sun is below ``horizon`` at ``time``, else `False`.
        """
        if not isinstance(time, Time):
            time = Time(time)

        solar_altitude = self.altaz(time, target=get_sun(time), obswl=obswl).alt
        return solar_altitude < horizon

class Target(object):
    """
    This is an abstract base class -- you can't instantiate
    examples of this class, but must work with one of its
    subclasses such as ``FixedTarget`` or ``NonFixedTarget``.

    Would need to import six, abc to make this a metaclass?
    """
    __metaclass__ = ABCMeta

    def __init__(self, name=None, ra=None, dec=None, marker=None):
        """
        Defines a single observation target.

        Parameters
        ----------
        name : str, optional

        ra : WHAT TYPE IS ra ?

        dec : WHAT TYPE IS dec ?

        marker : str, optional
            User-defined markers to differentiate between different types
            of targets (e.g., guides, high-priority, etc.).
        """
        raise NotImplementedError()

    @property
    def ra(self):
        """
        Right ascension.
        """
        if isinstance(self, FixedTarget):
            return self.coord.ra
        raise NotImplementedError()

    @property
    def dec(self):
        """
        Declination.
        """
        if isinstance(self, FixedTarget):
            return self.coord.dec
        raise NotImplementedError()


class FixedTarget(Target):
    """
    An object that is "fixed" with respect to the celestial sphere.
    """
    def __init__(self, coord, name=None, **kwargs):
        """
        TODO: Docstring.
        """
        if not (hasattr(coord, 'transform_to') and
                hasattr(coord, 'represent_as')):
            raise TypeError('`coord` must be a coordinate object.')

        self.name = name
        self.coord = coord

    @classmethod
    def from_name(cls, query_name, name=None, **kwargs):
        """
        Initialize a `FixedTarget` by querying for a name, using the machinery
        in `~astropy.coordinates.SkyCoord.from_name`.
        """
        # Allow manual override for name keyword so that the target name can
        # be different from the query name, otherwise assume name=queryname.
        if name is None:
            name = query_name
        return cls(SkyCoord.from_name(query_name), name=name, **kwargs)

class NonFixedTarget(Target):
    """
    Placeholder for future function.
    """


class Constraint(object):
    """
    An object containing observational constraints.

    A Constraints object is used in conjunction with a Target
    and an Observer object (via the apply_constraints method) to find out
    if a particular target is visible to the observer.
    """
    __metaclass__ = ABCMeta

    @abstractmethod
    def apply_constraints(self, target, observer, constraint_list):
        """
        Returns information on a target's visibility.

        Finds out if a Target is observable by an Observer given a list
        of Constraint objects.  The list must contain at least one
        Constraint object.

        Parameters
        ----------
        target : WHAT TYPE IS Target OBJECT ?

        constraint_list : WHAT TYPE IS constraint_list ? `numpy.array` ??
        """
        raise NotImplementedError


class TimeWindow(Constraint):
    """
    An object containing start and end times for an observation.
    """

    def __init__(self, start, end):
        """
        Initializes a TimeWindow object.

        Parameters
        ----------
        start : STRING OR astropy.time OBJECT ?

        end : STRING OR astropy.time OBJECT ?
        """
        raise NotImplementedError


class AltitudeRange(Constraint):
    """
    An object containing upper and lower altitude limits.
    """

    def __init__(self, low, high):
        """
        Initializes an AltitudeRange object.

        Parameters
        ----------
        low : `~astropy.units.Quantity`

        high : `~astropy.units.Quantity`
        """
        raise NotImplementedError


class AboveAirmass(Constraint):
    """
    An object containing an airmass lower limit.
    """

    def __init__(self, low):
        """
        Initializes an AboveAirmass object.

        Parameters
        ----------
        low : float
        """
        raise NotImplementedError


class Observation(object):
    """
    Comments.
    """

    def __init__(self, target, time):
        """
        Initializes an Observation object.

        Parameters
        ----------
        target : WHAT TYPE IS Target OBJECT ?

        date : WHAT TYPE IS date OBJECT ?
        """
        raise NotImplementedError()

    # Observability properties.

    @property
    def alt(self):
        """
        Altitude at time of observation.
        """
        raise NotImplementedError()

    @property
    def az(self):
        """
        Azimuth at time of observation.
        """
        raise NotImplementedError()

    @property
    def airmass(self):
        """
        Airmass.
        """
        raise NotImplementedError()

    @property
    def pang(self):
        """
        Parallactic angle.
        """
        raise NotImplementedError()

    @property
    def ha(self):
        """ 
        Hour angle.
        """
        raise NotImplementedError()

    @property
    def moon_sep(self):
        """
        Separation between moon and object at time of observation.
        """
        raise NotImplementedError()

    # Time properties.

    @property
    def ut(self):
        """
        Time of observation in UST.
        """
        raise NotImplementedError()

    @property
    def lt(self):
        """
        Time of observation in local time.
        """
        raise NotImplementedError()

    @property
    def gmst(self):
        """
        Time of observation in GMST.
        """
        raise NotImplementedError()

    @property
    def lmst(self):
        """
        Time of observation in local mean sidereal time.
        """
        raise NotImplementedError()
