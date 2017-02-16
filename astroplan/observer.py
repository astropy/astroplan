# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# Standard library
import datetime
import warnings

# Third-party
from astropy.coordinates import (EarthLocation, SkyCoord, AltAz, get_sun,
                                 get_moon, Angle, Latitude, Longitude,
                                 UnitSphericalRepresentation)
from astropy.extern.six import string_types
import astropy.units as u
from astropy.time import Time
from astropy.utils import isiterable
import numpy as np
import pytz

# Package
from .exceptions import TargetNeverUpWarning, TargetAlwaysUpWarning
from .moon import moon_illumination, moon_phase_angle

__all__ = ["Observer", "MAGIC_TIME"]

MAGIC_TIME = Time(-999, format='jd')


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


class Observer(object):

    """
    A container class for information about an observer's location and
    environment.

    Examples
    --------
    We can create an observer at Subaru Observatory in Hawaii two ways. First,
    locations for some observatories are stored in astroplan, and these can be
    accessed by name, like so:

    >>> from astroplan import Observer
    >>> subaru = Observer.at_site("Subaru", timezone="US/Hawaii")

    To find out which observatories can be accessed by name, check out
    `~astropy.coordinates.EarthLocation.get_site_names`.

    Next, you can initialize an observer by specifying the location with
    `~astropy.coordinates.EarthLocation`:

    >>> from astropy.coordinates import EarthLocation
    >>> import astropy.units as u
    >>> location = EarthLocation.from_geodetic(-155.4761*u.deg, 19.825*u.deg,
    ...                                        4139*u.m)
    >>> subaru = Observer(location=location, name="Subaru", timezone="US/Hawaii")

    You can also create an observer without an
    `~astropy.coordinates.EarthLocation`:

    >>> from astroplan import Observer
    >>> import astropy.units as u
    >>> subaru = Observer(longitude=-155.4761*u.deg, latitude=19.825*u.deg,
    ...                   elevation=0*u.m, name="Subaru", timezone="US/Hawaii")

    """
    @u.quantity_input(elevation=u.m)
    def __init__(self, location=None, timezone='UTC', name=None, latitude=None,
                 longitude=None, elevation=0*u.m, pressure=None,
                 relative_humidity=None, temperature=None, description=None):
        """
        Parameters
        ----------
        location : `~astropy.coordinates.EarthLocation`
            The location (latitude, longitude, elevation) of the observatory.

        timezone : str or `datetime.tzinfo` (optional)
            The local timezone to assume. If a string, it will be passed
            through ``pytz.timezone()`` to produce the timezone object.

        name : str
            A short name for the telescope, observatory or location.

        latitude : float, str, `~astropy.units.Quantity` (optional)
            The latitude of the observing location. Should be valid input for
            initializing a `~astropy.coordinates.Latitude` object.

        longitude : float, str, `~astropy.units.Quantity` (optional)
            The longitude of the observing location. Should be valid input for
            initializing a `~astropy.coordinates.Longitude` object.

        elevation : `~astropy.units.Quantity` (optional), default = 0 meters
            The elevation of the observing location, with respect to sea
            level. Defaults to sea level.

        pressure : `~astropy.units.Quantity` (optional)
            The ambient pressure. Defaults to zero (i.e. no atmosphere).

        relative_humidity : float (optional)
            The ambient relative humidity.

        temperature : `~astropy.units.Quantity` (optional)
            The ambient temperature.

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

    def __repr__(self):
        """
        String representation of the `~astroplan.Observer` object.

        Examples
        --------

        >>> from astroplan import Observer
        >>> keck = Observer.at_site("Keck", timezone="US/Hawaii")
        >>> print(keck)                                    # doctest: +FLOAT_CMP
        <Observer: name='Keck',
            location (lon, lat, el)=(-155.478333333 deg, 19.8283333333 deg, 4160.0 m),
            timezone=<DstTzInfo 'US/Hawaii' LMT-1 day, 13:29:00 STD>>
        """
        class_name = self.__class__.__name__
        attr_names = ['name', 'location', 'timezone', 'pressure', 'temperature',
                      'relative_humidity']
        attr_values = [getattr(self, attr) for attr in attr_names]
        attributes_strings = []
        for name, value in zip(attr_names, attr_values):
            if value is not None:
                # Format location for easy readability
                if name == 'location':
                    formatted_loc = ["{} {}".format(i.value, i.unit)
                                     for i in value.to_geodetic()]
                    attributes_strings.append(
                        "{} (lon, lat, el)=({})".format(name,
                                                        ", ".join(formatted_loc)))
                else:
                    if name != 'name':
                        value = repr(value)
                    else:
                        value = "'{}'".format(value)
                    attributes_strings.append("{}={}".format(name, value))
        return "<{}: {}>".format(class_name, ",\n    ".join(attributes_strings))

    @classmethod
    def at_site(cls, site_name, **kwargs):
        """
        Initialize an `~astroplan.observer.Observer` object with a site name.

        Extra keyword arguments are passed to the `~astroplan.Observer`
        constructor (see `~astroplan.Observer` for available keyword
        arguments).

        Parameters
        ----------
        site_name : str
            Observatory name, must be resolvable with
            `~astropy.coordinates.EarthLocation.get_site_names`.

        Returns
        -------
        `~astroplan.observer.Observer`
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
        return cls(location=EarthLocation.of_site(site_name), name=name, **kwargs)

    def astropy_time_to_datetime(self, astropy_time):
        """
        Convert the `~astropy.time.Time` object ``astropy_time`` to a
        localized `~datetime.datetime` object.

        Timezones localized with `pytz`_.

        .. _pytz: https://pypi.python.org/pypi/pytz/

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

        Examples
        --------
        Convert an astropy time to a localized `~datetime.datetime`:

        >>> from astroplan import Observer
        >>> from astropy.time import Time
        >>> subaru = Observer.at_site("Subaru", timezone="US/Hawaii")
        >>> astropy_time = Time('1999-12-31 06:00:00')
        >>> print(subaru.astropy_time_to_datetime(astropy_time))
        1999-12-30 20:00:00-10:00
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

        Timezones localized with `pytz`_. If the ``date_time`` is naive, the
        implied timezone is the ``timezone`` structure of ``self``.

        Parameters
        ----------
        date_time : `~datetime.datetime` or list-like

        Returns
        -------
        `~astropy.time.Time`
            Astropy time object (no timezone information preserved).

        Examples
        --------
        Convert a localized `~datetime.datetime` to a `~astropy.time.Time`
        object. Non-localized datetimes are assumed to be UTC.
        <Time object: scale='utc' format='datetime' value=1999-12-31 06:00:00>

        >>> from astroplan import Observer
        >>> import datetime
        >>> import pytz
        >>> subaru = Observer.at_site("Subaru", timezone="US/Hawaii")
        >>> hi_date_time = datetime.datetime(2005, 6, 21, 20, 0, 0, 0)
        >>> subaru.datetime_to_astropy_time(hi_date_time)
        <Time object: scale='utc' format='datetime' value=2005-06-22 06:00:00>
        >>> utc_date_time = datetime.datetime(2005, 6, 22, 6, 0, 0, 0,
        ...                                   tzinfo=pytz.timezone("UTC"))
        >>> subaru.datetime_to_astropy_time(utc_date_time)
        <Time object: scale='utc' format='datetime' value=2005-06-22 06:00:00>
        """

        if hasattr(date_time, '__iter__'):
            return Time([self.datetime_to_astropy_time(t) for t in date_time])

        # For timezone-naive datetimes, assign local timezone
        if date_time.tzinfo is None:
            date_time = self.timezone.localize(date_time)

        return Time(date_time, location=self.location)

    def _transform_target_list_to_altaz(self, times, targets):
        """
        Workaround for transforming a list of coordinates ``targets`` to
        altitudes and azimuths.

        Parameters
        ----------
        times : `~astropy.time.Time` or list of `~astropy.time.Time` objects
            Time of observation

        targets : `~astropy.coordinates.SkyCoord` or list of `~astropy.coordinates.SkyCoord` objects
            List of target coordinates

        Returns
        -------
        altitudes : list
            List of altitudes for each target, at each time
        """
        if times.isscalar:
            times = Time([times])

        if not isinstance(targets, list) and targets.isscalar:
            targets = [targets]

        targets_is_unitsphericalrep = [x.data.__class__ is
                                       UnitSphericalRepresentation for x in targets]
        if all(targets_is_unitsphericalrep) or not any(targets_is_unitsphericalrep):
            repeated_times = np.tile(times, len(targets))
            ra_list = Longitude([x.icrs.ra for x in targets])
            dec_list = Latitude([x.icrs.dec for x in targets])
            repeated_ra = np.repeat(ra_list, len(times))
            repeated_dec = np.repeat(dec_list, len(times))
            inner_sc = SkyCoord(ra=repeated_ra, dec=repeated_dec)
            target_SkyCoord = SkyCoord(inner_sc.data.represent_as(UnitSphericalRepresentation),
                                       representation=UnitSphericalRepresentation)
            transformed_coord = target_SkyCoord.transform_to(AltAz(location=self.location,
                                                                   obstime=repeated_times))
        else:
            # TODO: This is super slow.
            repeated_times = np.tile(times, len(targets))
            repeated_targets = np.repeat(targets, len(times))
            target_SkyCoord = SkyCoord(SkyCoord(repeated_targets).data.represent_as(
                                       UnitSphericalRepresentation),
                                       representation=UnitSphericalRepresentation)

            transformed_coord = target_SkyCoord.transform_to(AltAz(location=self.location,
                                                                   obstime=repeated_times))
        return transformed_coord

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

        Examples
        --------
        Create an instance of the `~astropy.coordinates.AltAz` frame for an
        observer at Apache Point Observatory at a particular time:

        >>> from astroplan import Observer
        >>> from astropy.time import Time
        >>> from astropy.coordinates import SkyCoord
        >>> apo = Observer.at_site("APO")
        >>> time = Time('2001-02-03 04:05:06')
        >>> target = SkyCoord(0*u.deg, 0*u.deg)
        >>> altaz_frame = apo.altaz(time)

        Now transform the target's coordinates to the alt/az frame:

        >>> target_altaz = target.transform_to(altaz_frame) # doctest: +SKIP

        Alternatively, construct an alt/az frame and transform the target to
        that frame all in one step:

        >>> target_altaz = apo.altaz(time, target) # doctest: +SKIP
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
            if isiterable(target) and not isinstance(target, SkyCoord):
                get_coord = lambda x: x.coord if hasattr(x, 'coord') else x
                transformed_coords = self._transform_target_list_to_altaz(time,
                                                                          list(map(get_coord, target)))
                n_targets = len(target)
                new_shape = (n_targets, int(len(transformed_coords)/n_targets))

                for comp in transformed_coords.data.components:
                    getattr(transformed_coords.data, comp).resize(new_shape)
                return transformed_coords

            # If target is a FixedTarget or a SkyCoord:
            if hasattr(target, 'coord'):
                coordinate = target.coord
            else:
                coordinate = target
            return coordinate.transform_to(altaz_frame)

    def parallactic_angle(self, time, target):
        """
        Calculate the parallactic angle.

        Parameters
        ----------
        time : `~astropy.time.Time`
            Observation time.

        target : `~astroplan.FixedTarget` or `~astropy.coordinates.SkyCoord` or list
            Target celestial object(s).

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

        """
        if not isinstance(time, Time):
            time = Time(time)

        if isiterable(target):
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
                      (np.tan(self.location.latitude.radian) *
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
        alt = np.atleast_2d(Latitude(alt))
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
                                                  assume_unique=True)  # [0]

            warnmsg = ('Target(s) index {} does not cross horizon={} within '
                       '24 hours'.format(noncrossing_target_ind, horizon))

            if (alt[noncrossing_target_ind, :] > horizon).all():
                warnings.warn(warnmsg, TargetAlwaysUpWarning)
            else:
                warnings.warn(warnmsg, TargetNeverUpWarning)

            # Fill in missing time with MAGIC_TIME
            target_inds = np.insert(target_inds, noncrossing_target_ind,
                                    noncrossing_target_ind)
            time_inds = np.insert(time_inds.astype(float),
                                  noncrossing_target_ind,
                                  np.nan)
        elif np.count_nonzero(condition) > n_targets:
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

        times = [t[int(i):int(i+2)] if not np.isnan(i) else np.nan for i in time_inds]
        altitudes = [alt[int(i), int(j):int(j+2)] if not np.isnan(j) else np.nan
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
        if not isinstance(times, Time):
            return MAGIC_TIME
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

        target_is_vector = isiterable(target)

        if prev_next == 'next':
            times = _generate_24hr_grid(time, 0, 1, N)
        else:
            times = _generate_24hr_grid(time, -1, 0, N)

        altaz = self.altaz(times, target)
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

        target_is_vector = isiterable(target)

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
            d_altitudes = [each_alt.diff() for each_alt in altaz.alt]
        else:
            altitudes = altaz.alt
            d_altitudes = altitudes.diff()

        dt = Time((times.jd[1:] + times.jd[:-1])/2, format='jd')

        horizon = 0*u.degree  # Find when derivative passes through zero
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

    def _determine_which_event(self, function, args_dict):
        """
        Run through the next/previous/nearest permutations of the solutions
        to `function(time, ...)`, and return the previous/next/nearest one
        specified by the args stored in args_dict.
        """
        time = args_dict.pop('time', None)
        target = args_dict.pop('target', None)
        which = args_dict.pop('which', None)
        horizon = args_dict.pop('horizon', None)
        rise_set = args_dict.pop('rise_set', None)
        antitransit = args_dict.pop('antitransit', None)

        # Assemble arguments for function, depending on the function.
        if function == self._calc_riseset:
            args = lambda w: (time, target, w, rise_set, horizon)
        elif function == self._calc_transit:
            args = lambda w: (time, target, w, antitransit)
        else:
            raise ValueError('Function {} not supported in '
                             '_determine_which_event.'.format(function))

        if not isinstance(time, Time):
            time = Time(time)

        if which == 'next' or which == 'nearest':
            next_event = function(*args('next'))
            if which == 'next':
                return next_event

        if which == 'previous' or which == 'nearest':
            previous_event = function(*args('previous'))
            if which == 'previous':
                return previous_event

        if which == 'nearest':
            if isiterable(target):
                return_times = []
                for next_e, prev_e in zip(next_event, previous_event):
                    if abs(time - prev_e) < abs(time - next_e):
                        return_times.append(prev_e)
                    else:
                        return_times.append(next_e)
                return Time(return_times)
            else:
                if abs(time - previous_event) < abs(time - next_event):
                    return previous_event
                else:
                    return next_event

        raise ValueError('"which" kwarg must be "next", "previous" or '
                         '"nearest".')

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

        target : coordinate object (i.e. `~astropy.coordinates.SkyCoord`, `~astroplan.FixedTarget`) or list
            Target celestial object(s)

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

        Examples
        --------
        Calculate the rise time of Rigel at Keck Observatory:

        >>> from astroplan import Observer, FixedTarget
        >>> from astropy.time import Time
        >>> time = Time("2001-02-03 04:05:06")
        >>> target = FixedTarget.from_name("Rigel")
        >>> keck = Observer.at_site("Keck")
        >>> rigel_rise_time = keck.target_rise_time(time, target, which="next")
        >>> print("ISO: {0.iso}, JD: {0.jd}".format(rigel_rise_time)) # doctest: +FLOAT_CMP
        ISO: 2001-02-04 00:51:23.330, JD: 2451944.53569
        """
        return self._determine_which_event(self._calc_riseset,
                                           dict(time=time, target=target,
                                                which=which, rise_set='rising',
                                                horizon=horizon))

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

        target : coordinate object (i.e. `~astropy.coordinates.SkyCoord`, `~astroplan.FixedTarget`) or list
            Target celestial object(s)

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

        Examples
        --------
        Calculate the set time of Rigel at Keck Observatory:

        >>> from astroplan import Observer, FixedTarget
        >>> from astropy.time import Time
        >>> time = Time("2001-02-03 04:05:06")
        >>> target = FixedTarget.from_name("Rigel")
        >>> keck = Observer.at_site("Keck")
        >>> rigel_set_time = keck.target_set_time(time, target, which="next")
        >>> print("ISO: {0.iso}, JD: {0.jd}".format(rigel_set_time)) # doctest: +FLOAT_CMP
        ISO: 2001-02-03 12:29:34.768, JD: 2451944.02054
        """
        return self._determine_which_event(self._calc_riseset,
                                           dict(time=time, target=target,
                                                which=which, rise_set='setting',
                                                horizon=horizon))

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

        target : coordinate object (i.e. `~astropy.coordinates.SkyCoord`, `~astroplan.FixedTarget`) or list
            Target celestial object(s)

        which : {'next', 'previous', 'nearest'}
            Choose which sunrise relative to the present ``time`` would you
            like to calculate

        Returns
        -------
        `~astropy.time.Time`
            Transit time of target

        Examples
        --------
        Calculate the meridian transit time of Rigel at Keck Observatory:

        >>> from astroplan import Observer, FixedTarget
        >>> from astropy.time import Time
        >>> time = Time("2001-02-03 04:05:06")
        >>> target = FixedTarget.from_name("Rigel")
        >>> keck = Observer.at_site("Keck")
        >>> rigel_transit_time = keck.target_meridian_transit_time(time, target,
        ...                                                        which="next")
        >>> print("ISO: {0.iso}, JD: {0.jd}".format(rigel_transit_time)) # doctest: +FLOAT_CMP
        ISO: 2001-02-03 06:42:26.863, JD: 2451943.77948
        """
        return self._determine_which_event(self._calc_transit,
                                           dict(time=time, target=target,
                                                which=which,
                                                rise_set='setting'))

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

        target : coordinate object (i.e. `~astropy.coordinates.SkyCoord`, `~astroplan.FixedTarget`) or list
            Target celestial object(s)

        which : {'next', 'previous', 'nearest'}
            Choose which sunrise relative to the present ``time`` would you
            like to calculate

        Returns
        -------
        `~astropy.time.Time`
            Antitransit time of target

        Examples
        --------
        Calculate the meridian anti-transit time of Rigel at Keck Observatory:

        >>> from astroplan import Observer, FixedTarget
        >>> from astropy.time import Time
        >>> time = Time("2001-02-03 04:05:06")
        >>> target = FixedTarget.from_name("Rigel")
        >>> keck = Observer.at_site("Keck")
        >>> rigel_antitransit_time = keck.target_meridian_antitransit_time(time, target,
        ...                                                                which="next")
        >>> print("ISO: {0.iso}, JD: {0.jd}".format(rigel_antitransit_time)) # doctest: +FLOAT_CMP
        ISO: 2001-02-03 18:40:29.761, JD: 2451944.27812

        """
        return self._determine_which_event(self._calc_transit,
                                           dict(time=time, target=target,
                                                which=which, antitransit=True,
                                                rise_set='setting'))

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

        Examples
        --------
        Calculate the time of the previous sunrise at Apache Point Observatory:

        >>> from astroplan import Observer
        >>> from astropy.time import Time
        >>> apo = Observer.at_site("APO")
        >>> time = Time('2001-02-03 04:05:06')
        >>> sun_rise = apo.sun_rise_time(time, which="previous")
        >>> print("ISO: {0.iso}, JD: {0.jd}".format(sun_rise)) # doctest: +SKIP
        ISO: 2001-02-02 14:02:50.554, JD: 2451943.08531
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

        Examples
        --------
        Calculate the time of the next sunset at Apache Point Observatory:

        >>> from astroplan import Observer
        >>> from astropy.time import Time
        >>> apo = Observer.at_site("APO")
        >>> time = Time('2001-02-03 04:05:06')
        >>> sun_set = apo.sun_set_time(time, which="next")
        >>> print("ISO: {0.iso}, JD: {0.jd}".format(sun_set)) # doctest: +SKIP
        ISO: 2001-02-04 00:35:42.102, JD: 2451944.52479
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

    def twilight_evening_astronomical(self, time, which='nearest'):
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

    def twilight_evening_nautical(self, time, which='nearest'):
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

    def twilight_evening_civil(self, time, which='nearest'):
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

    def twilight_morning_astronomical(self, time, which='nearest'):
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

    def twilight_morning_nautical(self, time, which='nearest'):
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

    def twilight_morning_civil(self, time, which='nearest'):
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

        Keywords : str, optional
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

        Keywords : str, optional
            previous
            next
        """
        raise NotImplementedError()

    def moon_illumination(self, time):
        """
        Calculate the illuminated fraction of the moon.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        Returns
        -------
        float
            Fraction of lunar surface illuminated

        Examples
        --------
        How much of the lunar surface is illuminated at 2015-08-29 18:35 UTC,
        which we happen to know is the time of a full moon?

        >>> from astroplan import Observer
        >>> from astropy.time import Time
        >>> apo = Observer.at_site("APO")
        >>> time = Time("2015-08-29 18:35")
        >>> apo.moon_illumination(time) # doctest: +SKIP
        array([ 0.99972487])
        """
        if not isinstance(time, Time):
            time = Time(time)

        return moon_illumination(time)

    def moon_phase(self, time=None):
        """
        Calculate lunar orbital phase.

        For example, phase=2*pi is "new", phase=0 is "full".

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        Returns
        -------
        moon_phase_angle : float
            Orbital phase angle of the moon where 2*pi corresponds to new moon,
            zero corresponds to full moon.

        Examples
        --------
        Calculate the phase of the moon at 2015-08-29 18:35 UTC. Near zero
        radians corresponds to a nearly full moon.

        >>> from astroplan import Observer
        >>> from astropy.time import Time
        >>> apo = Observer.at_site('APO')
        >>> time = Time('2015-08-29 18:35')
        >>> apo.moon_phase(time) # doctest: +SKIP
        <Quantity [ 0.03317537] rad>
        """
        if time is not None and not isinstance(time, Time):
            time = Time(time)

        return moon_phase_angle(time)

    def moon_altaz(self, time, ephemeris=None):
        """
        Returns the position of the moon in alt/az.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        ephemeris : str, optional
            Ephemeris to use.  If not given, use the one set with
            ``astropy.coordinates.solar_system_ephemeris.set`` (which is
            set to 'builtin' by default).


        Returns
        -------
        altaz : `~astropy.coordinates.SkyCoord`
            Position of the moon transformed to altitude and azimuth

        Examples
        --------
        Calculate the altitude and azimuth of the moon at Apache Point
        Observatory:

        >>> from astroplan import Observer
        >>> from astropy.time import Time
        >>> apo = Observer.at_site("APO")
        >>> time = Time("2015-08-29 18:35")
        >>> altaz_moon = apo.moon_altaz(time) # doctest: +SKIP
        >>> print("alt: {0.alt}, az: {0.az}".format(altaz_moon)) # doctest: +SKIP
        alt: -63.72706397691421 deg, az: 345.3640380598265 deg
        """
        if not isinstance(time, Time):
            time = Time(time)

        # TODO: when astropy/astropy#5069 is resolved, replace this workaround which
        # handles scalar and non-scalar time inputs differently

        if time.isscalar:
            altaz_frame = AltAz(location=self.location, obstime=time)
            sun = get_sun(time).transform_to(altaz_frame)

            moon = get_moon(time, location=self.location, ephemeris=ephemeris).transform_to(altaz_frame)
            return moon

        else:
            moon_coords = []
            for t in time:
                altaz_frame = AltAz(location=self.location, obstime=t)
                moon_coord = get_moon(t, location=self.location, ephemeris=ephemeris).transform_to(altaz_frame)
                moon_coords.append(moon_coord)
            obstime = [coord.obstime for coord in moon_coords]
            alts = u.Quantity([coord.alt for coord in moon_coords])
            azs = u.Quantity([coord.az for coord in moon_coords])
            dists = u.Quantity([coord.distance for coord in moon_coords])
            return SkyCoord(AltAz(azs, alts, dists, obstime=obstime, location=self.location))

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

        target : coordinate object (i.e. `~astropy.coordinates.SkyCoord`, `~astroplan.FixedTarget`) or list
            Target celestial object(s)

        horizon : `~astropy.units.Quantity` (optional), default = zero degrees
            Degrees above/below actual horizon to use
            for calculating rise/set times (i.e.,
            -6 deg horizon = civil twilight, etc.)

        return_altaz : bool (optional)
            Also return the '~astropy.coordinates.AltAz' coordinate.

        Returns
        -------
        observable : boolean
            True if ``target`` is above ``horizon`` at ``time``, else False.

        Examples
        --------
        Are Aldebaran and Vega above the horizon at Apache Point Observatory
        at 2015-08-29 18:35 UTC?

        >>> from astroplan import Observer, FixedTarget
        >>> from astropy.time import Time
        >>> apo = Observer.at_site("APO")
        >>> time = Time("2015-08-29 18:35")
        >>> aldebaran = FixedTarget.from_name("Aldebaran")
        >>> vega = FixedTarget.from_name("Vega")
        >>> apo.target_is_up(time, aldebaran)
        True
        >>> apo.target_is_up(time, [aldebaran, vega])
        [True, False]
        """
        if not isinstance(time, Time):
            time = Time(time)

        altaz = self.altaz(time, target)
        if isiterable(target):
            observable = [bool(alt > horizon) for alt in altaz.alt]
        else:
            observable = bool(altaz.alt > horizon)

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

        Examples
        --------
        Is it "nighttime" (i.e. is the Sun below ``horizon``) at Apache Point
        Observatory at 2015-08-29 18:35 UTC?

        >>> from astroplan import Observer
        >>> from astropy.time import Time
        >>> apo = Observer.at_site("APO")
        >>> time = Time("2015-08-29 18:35")
        >>> apo.is_night(time) # doctest: +SKIP
        False
        """
        if not isinstance(time, Time):
            time = Time(time)

        solar_altitude = self.altaz(time, target=get_sun(time), obswl=obswl).alt
        return bool(solar_altitude < horizon)

    def local_sidereal_time(self, time, kind='apparent', model=None):
        """
        Convert ``time`` to local sidereal time for observer.

        This is a thin wrapper around the `~astropy.time.Time.sidereal_time`
        method.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observation. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object)

        kind : {'mean', 'apparent'} (optional)
            Passed to the ``kind`` argument of
            `~astropy.time.Time.sidereal_time`

        model : str or `None`; optional
            The precession/nutation model to assume - see
            `~astropy.time.Time.sidereal_time` for more details.

        Returns
        -------
        `~astropy.coordinates.Longitude`
            Local sidereal time.
        """
        if not isinstance(time, Time):
            time = Time(time)

        return time.sidereal_time(kind, longitude=self.location.longitude,
                                  model=model)

    def target_hour_angle(self, time, target):
        """
        Calculate the local hour angle of ``target`` at ``time``.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observation. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object)

        target : coordinate object (i.e. `~astropy.coordinates.SkyCoord`, `~astroplan.FixedTarget`) or list
            Target celestial object(s)

        Returns
        -------
        hour_angle : `~astropy.coordinates.Angle`
            The hour angle(s) of the target(s) at ``time``
        """
        if not isinstance(time, Time):
            time = Time(time)

        if isiterable(target):
            coords = [t.coord if hasattr(t, 'coord') else t
                      for t in target]

            hour_angle = Longitude([self.local_sidereal_time(time) - coord.ra
                                    for coord in coords])

        else:
            coord = target.coord if hasattr(target, 'coord') else target
            hour_angle = Longitude(self.local_sidereal_time(time) - coord.ra)

        return hour_angle

    @u.quantity_input(horizon=u.degree)
    def tonight(self, time=None, horizon=0 * u.degree, obswl=None):
        """
        Return a time range corresponding to the nearest night

        This will return a range of `~astropy.time.Time` corresponding to the
        beginning and ending of the night. If in the middle of a given night,
        return times from `~astropy.time.Time.now` until the nearest
        `~astroplan.Observer.sun_rise_time`

        Parameters
        ----------
        time : `~astropy.time.Time` (optional), default = `~astropy.time.Time.now`
            The start time for tonight, which is allowed to be arbitrary. See description
            above for behavior
        horizon : `~astropy.units.Quantity` (optional), default = zero degrees
            Degrees above/below actual horizon to use for calculating rise/set times
            (e.g., -6 deg horizon = civil twilight, etc.)
        obswl : `~astropy.units.Quantity` (optional)
            Wavelength of the observation used in the calculation

        Returns
        -------
        times : `~astropy.time.Time`
            A tuple of times corresponding to the start and end of current night
        """
        current_time = Time.now() if time is None else time
        if self.is_night(current_time, horizon=horizon, obswl=obswl):
            start_time = current_time
        else:
            start_time = self.sun_set_time(current_time, which='next', horizon=horizon)

        end_time = self.sun_rise_time(current_time, which='next', horizon=horizon)

        return start_time, end_time
