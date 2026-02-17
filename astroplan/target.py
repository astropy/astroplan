# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Standard library
from abc import ABCMeta
import warnings

# Third-party
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, ICRS, UnitSphericalRepresentation, AltAz
try:
    from sgp4.io import twoline2rv
    from sgp4.earth_gravity import wgs84 as sgp4_wgs84
    from skyfield.api import load, wgs84
    from skyfield.sgp4lib import EarthSatellite
    skyfield_available = True
except ImportError:
    skyfield_available = False

# Package
from .exceptions import InvalidTLEDataWarning


__all__ = ["Target", "FixedTarget", "TLETarget", "NonFixedTarget"]

# Docstring code examples include printed SkyCoords, but the format changed
# in astropy 1.3. Thus the doctest needs astropy >=1.3 and this is the
# easiest way to make it work.

__doctest_requires__ = {'FixedTarget.*': ['astropy.modeling.Hermite1D']}


class Target:
    """
    Abstract base class for target objects.

    This is an abstract base class -- you can't instantiate
    examples of this class, but must work with one of its
    subclasses such as `~astroplan.target.FixedTarget` or
    `~astroplan.target.NonFixedTarget`.
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
    Coordinates and metadata for an object that is "fixed" with respect to the
    celestial sphere.

    Examples
    --------
    Create a `~astroplan.FixedTarget` object for Sirius:

    >>> from astroplan import FixedTarget
    >>> from astropy.coordinates import SkyCoord
    >>> import astropy.units as u
    >>> sirius_coord = SkyCoord(ra=101.28715533*u.deg, dec=16.71611586*u.deg)
    >>> sirius = FixedTarget(coord=sirius_coord, name="Sirius")

    Create an equivalent `~astroplan.FixedTarget` object for Sirius by querying
    for the coordinates of Sirius by name:

    >>> from astroplan import FixedTarget
    >>> sirius = FixedTarget.from_name("Sirius")  # doctest: +REMOTE_DATA
    """

    def __init__(self, coord, name=None, **kwargs):
        """
        Parameters
        ----------
        coord : `~astropy.coordinates.SkyCoord`
            Coordinate of the target

        name : str (optional)
            Name of the target, used for plotting and representing the target
            as a string
        """
        if not (hasattr(coord, 'transform_to') and
                hasattr(coord, 'represent_as')):
            raise TypeError('`coord` must be a coordinate object.')

        self.name = name
        self.coord = coord

    @classmethod
    def from_name(cls, query_name, name=None, **kwargs):
        """
        Initialize a `FixedTarget` by querying for a name from the CDS name
        resolver, using the machinery in
        `~astropy.coordinates.SkyCoord.from_name`.

        This

        Parameters
        ----------
        query_name : str
            Name of the target used to query for coordinates.

        name : string or `None`
            Name of the target to use within astroplan. If `None`, query_name
            is used as ``name``.

        Examples
        --------
        >>> from astroplan import FixedTarget
        >>> sirius = FixedTarget.from_name("Sirius")  # doctest: +REMOTE_DATA
        >>> sirius.coord                              # doctest: +FLOAT_CMP +REMOTE_DATA
        <SkyCoord (ICRS): (ra, dec) in deg
            ( 101.28715533, -16.71611586)>
        """
        # Allow manual override for name keyword so that the target name can
        # be different from the query name, otherwise assume name=queryname.
        if name is None:
            name = query_name
        return cls(SkyCoord.from_name(query_name), name=name, **kwargs)

    def __repr__(self):
        """
        String representation of `~astroplan.FixedTarget`.

        Examples
        --------
        Show string representation of a `~astroplan.FixedTarget` for Vega:

        >>> from astroplan import FixedTarget
        >>> from astropy.coordinates import SkyCoord
        >>> vega_coord = SkyCoord(ra='279.23473479d', dec='38.78368896d')
        >>> vega = FixedTarget(coord=vega_coord, name="Vega")
        >>> print(vega)                             # doctest: +FLOAT_CMP
        <FixedTarget "Vega" at SkyCoord (ICRS): (ra, dec) in deg ( 279.23473479, 38.78368894)>
        """
        class_name = self.__class__.__name__
        fmt_coord = repr(self.coord).replace('\n   ', '')[1:-1]
        return '<{} "{}" at {}>'.format(class_name, self.name, fmt_coord)


class NonFixedTarget(Target):
    """
    Placeholder for future function.
    """


class TLETarget(Target):
    """
    A target defined by TLE (Two-Line Element set) for satellites.
    """
    def __init__(self, line1, line2, name=None, observer=None, skip_tle_check=False):
        """
        Parameters
        ----------
        line1 : str
            The first line of the TLE set

        line2 : str
            The second line of the TLE set

        name : str (optional)
            Name of the target, used for plotting and representing the target
            as a string

        observer : `~astroplan.Observer` (optional)
            The location of observer.
            If `None`, the observer is assumed to be at sea level at the equator.

        skip_tle_check : bool (optional)
            Whether to skip TLE validation
        """
        if not skyfield_available:
            raise ImportError("Please install the skyfield package to use the TLETarget class.")

        if not skip_tle_check:
            twoline2rv(line1, line2, sgp4_wgs84)    # Raises ValueError if TLE is invalid

        self.name = name
        self.satellite = EarthSatellite(line1, line2, name, load.timescale())

        if observer is None:
            # Prevent circular import and usually not used
            from .observer import Observer
            self.observer = Observer(latitude=0*u.deg, longitude=0*u.deg, elevation=0*u.m)
        else:
            self.observer = observer

        longitude, latitude, height = self.observer.location.to_geodetic()
        self.geographic_position = wgs84.latlon(latitude.to(u.deg).value,
                                                longitude.to(u.deg).value,
                                                height.to(u.m).value)

    @classmethod
    def from_string(cls, tle_string, name=None, *args, **kwargs):
        """
        Creates a TLETarget instance from a complete TLE string.

        Parameters
        ----------
        tle_string : str
            String to be parsed, expected to contain 2 or 3 newline-separated lines.

        name : str (optional)
            Name of the target. If not provided and the tle_string contains 3 lines,
            the first line will be used as the name.

        args, kwargs : tuple, dict (optional)
            Additional arguments and keyword arguments to be passed to the TLETarget constructor.
        """
        lines = tle_string.strip().splitlines()

        if len(lines) not in (2, 3):
            raise ValueError(f"Expected TLE string to contain 2 or 3 lines, got {len(lines)}")

        if len(lines) == 3:
            line1, line2, name = lines[1], lines[2], name or lines[0]
        else:   # len(lines) == 2
            line1, line2 = lines
        return cls(line1, line2, name, *args, **kwargs)

    def _compute_topocentric(self, times=None):
        """
        Compute the topocentric coordinates (relative to observer) at a particular time.

        Parameters
        ----------
        times : `~astropy.time.Time` (optional)
            The time(s) to use in the calculation.

        Returns
        -------
        topocentric : `skyfield.positionlib.ICRF`
            The topocentric object representing the relative coordinates of the target.
        """
        if times is None:
            times = Time.now()
        ts = load.timescale()
        t = ts.from_astropy(times)

        topocentric = (self.satellite - self.geographic_position).at(t)

        # Check for invalid TLE data. A non-None usually message means the computation went beyond
        # the physically sensible point. Details:
        # https://rhodesmill.org/skyfield/earth-satellites.html#detecting-propagation-errors
        message = topocentric.message
        if (
            (message is not None and not isinstance(message, list)) or
            (isinstance(message, list) and not all(x is None for x in message))
        ):
            warnings.warn(f"Invalid TLE Data: {message}", InvalidTLEDataWarning)

        return topocentric

    def get_skycoord(self, times=None):
        """
        Get the coordinates of the target at a particular time.

        Parameters
        ----------
        times : `~astropy.time.Time` (optional)
            The time(s) to use in the calculation.

        Returns
        -------
        coord : `~astropy.coordinates.SkyCoord`
            A single SkyCoord object, which may be non-scalar, representing the target's
            RA/Dec coordinates at the specified time(s). Might return np.nan and output a
            warning for times where the elements stop making physical sense.
        """
        topocentric = self._compute_topocentric(times)
        ra, dec, distance = topocentric.radec()
        # No distance, in SkyCoord, distance is from frame origin, but here, it's from observer.
        return SkyCoord(
            ra.hours*u.hourangle,
            dec.degrees*u.deg,
            obstime=times,
            frame='icrs',
            location=self.observer.location,
        )

    def altaz(self, times=None):
        """
        Get the altitude and azimuth of the target at a particular time.

        Parameters
        ----------
        times : `~astropy.time.Time` (optional)
            The time(s) to use in the calculation.

        Returns
        -------
        altaz_coord : `~astropy.coordinates.SkyCoord`
            SkyCoord object representing the target's altitude and azimuth at the specified time(s)
        """
        topocentric = self._compute_topocentric(times)

        temperature_C = None
        pressure_mbar = None
        if self.observer.temperature is not None:
            temperature_C = self.observer.temperature.to_value(u.deg_C)
        if self.observer.pressure is not None:
            pressure_mbar = self.observer.pressure.to_value(u.mbar)

        alt, az, distance = topocentric.altaz(
            temperature_C=temperature_C,
            pressure_mbar=pressure_mbar,
        )

        # 'relative_humidity' and 'obswl' were not used in coordinate calculation
        altaz_frame = AltAz(
            location=self.observer.location,
            obstime=times,
            pressure=self.observer.pressure,
            temperature=self.observer.temperature,
        )

        return SkyCoord(alt=alt.degrees*u.deg, az=az.degrees*u.deg, frame=altaz_frame)

    def __repr__(self):
        return f'<{self.__class__.__name__} "{self.name}">'

    def __str__(self):
        return self.name


def get_skycoord(targets, times=None):
    """
    Return an `~astropy.coordinates.SkyCoord` object.

    When performing calculations it is usually most efficient to have
    a single `~astropy.coordinates.SkyCoord` object, rather than a
    list of `Target` or `~astropy.coordinates.SkyCoord` objects.

    This is a convenience routine to do that, and it also supports targets
    that require evaluation at specific times (e.g., AltAz-defined targets).

    Parameters
    ----------
    targets : list, `~astropy.coordinates.SkyCoord`, `~astroplan.Target`
        Either a single target or a list of targets.

    times : `~astropy.time.Time` or time-like (optional)
        Times at which to evaluate time-dependent targets. Required if any
        target in ``targets`` needs evaluation at a time.

    Returns
    -------
    coord : `~astropy.coordinates.SkyCoord`
        A single SkyCoord object, which may be non-scalar. If ``times``
        is provided and any target is time-dependent, coordinates are broadcast
        or evaluated across time along subsequent axes.
    """
    if times is not None and not isinstance(times, Time):
        times = Time(times)

    def _is_time_dependent(obj):
        return callable(getattr(obj, "get_skycoord", None)) and not hasattr(obj, "coord")

    def _as_coord(obj):
        if hasattr(obj, "coord"):
            return obj.coord
        if callable(getattr(obj, "get_skycoord", None)):
            return obj.get_skycoord(times)
        return obj

    # Ignore non-scalar SkyCoords targets here
    # e.g. from get_body/get_sun, because they represent a single target
    is_multiple_targets = isinstance(targets, (list, tuple))
    if not is_multiple_targets:
        return _as_coord(targets)

    coords = [_as_coord(t) for t in targets]

    # If any target is time dependent, broadcast fixed coords to match times.shape
    time_dependent = (times is not None) and any(_is_time_dependent(t) for t in targets)
    times_shape = times.shape if time_dependent else None

    def _broadcast_quantity(q, shape):
        """Broadcast quantity to target shape if needed."""
        if shape is None or q.shape == shape:
            return q
        return u.Quantity(np.broadcast_to(q.to_value(q.unit), shape), q.unit)

    # Are all SkyCoord's in equivalent frames? If not, convert to ICRS
    convert_to_icrs = not all(
        [coord.frame.is_equivalent_frame(coords[0].frame) for coord in coords[1:]]
    )

    # we also need to be careful about handling mixtures of
    # UnitSphericalRepresentations and others
    targets_is_unitsphericalrep = [x.data.__class__ is
                                   UnitSphericalRepresentation for x in coords]

    longitudes = []
    latitudes = []
    distances = []
    get_distances = not all(targets_is_unitsphericalrep)
    if convert_to_icrs:
        # mixture of frames
        for coordinate in coords:
            icrs_coordinate = coordinate.icrs
            lon = icrs_coordinate.ra
            lat = icrs_coordinate.dec
            if times_shape is not None:
                lon = _broadcast_quantity(lon, times_shape)
                lat = _broadcast_quantity(lat, times_shape)
            longitudes.append(lon)
            latitudes.append(lat)
            if get_distances:
                dist = icrs_coordinate.distance
                if times_shape is not None:
                    dist = _broadcast_quantity(dist, times_shape)
                distances.append(dist)
        frame = ICRS()
    else:
        # all the same frame, get the longitude and latitude names
        try:
            # from astropy v2.0, keys are classes
            lon_name, lat_name = [
                mapping.framename for mapping in
                coords[0].frame_specific_representation_info[UnitSphericalRepresentation]]
        except BaseException:            # whereas prior to that they were strings.
            lon_name, lat_name = [mapping.framename for mapping in
                                  coords[0].frame_specific_representation_info['spherical']]

        frame = coords[0].frame
        for coordinate in coords:
            lon = getattr(coordinate, lon_name)
            lat = getattr(coordinate, lat_name)
            if times_shape is not None:
                lon = _broadcast_quantity(lon, times_shape)
                lat = _broadcast_quantity(lat, times_shape)
            longitudes.append(lon)
            latitudes.append(lat)
            if get_distances:
                dist = coordinate.distance
                if times_shape is not None:
                    dist = _broadcast_quantity(dist, times_shape)
                distances.append(dist)

    # Convert all longitude/latitude quantities to a common unit
    # and plain ndarrays before stacking (robust across units/Quantity subclasses).
    lon_unit = longitudes[0].unit
    lat_unit = latitudes[0].unit
    lon_vals = np.stack([lon.to_value(lon_unit) for lon in longitudes], axis=0)
    lat_vals = np.stack([lat.to_value(lat_unit) for lat in latitudes], axis=0)
    lon_q = u.Quantity(lon_vals, unit=lon_unit)
    lat_q = u.Quantity(lat_vals, unit=lat_unit)

    # now let's deal with the fact that we may have a mixture of coords with distances and
    # coords with UnitSphericalRepresentations
    if all(targets_is_unitsphericalrep):
        return SkyCoord(lon_q, lat_q, frame=frame)

    if not any(targets_is_unitsphericalrep):
        dist_unit = distances[0].unit
        dist_vals = np.stack([d.to_value(dist_unit) for d in distances], axis=0)
        dist_q = u.Quantity(dist_vals, unit=dist_unit)
        return SkyCoord(lon_q, lat_q, dist_q, frame=frame)

    # Mixture of coords with distances and without.
    # Assign large distances to UnitSphericalRepresentation objects.
    filled_distances = []
    for dist, is_unitspherical in zip(distances, targets_is_unitsphericalrep):
        if is_unitspherical:
            fill_vals = np.broadcast_to(100.0, dist.shape if dist.shape else ())
            filled_distances.append(u.Quantity(fill_vals, u.kpc))
        else:
            filled_distances.append(dist)

    dist_unit = filled_distances[0].unit
    dist_vals = np.stack([d.to_value(dist_unit) for d in filled_distances], axis=0)
    dist_q = u.Quantity(dist_vals, unit=dist_unit)
    return SkyCoord(lon_q, lat_q, dist_q, frame=frame)


class SpecialObjectFlag:
    """
    Flag this object as a special non-fixed target, which has a ``get_*`` method
    within astropy (like the Sun or Moon)
    """
    pass


class SunFlag(SpecialObjectFlag):
    """
    Flag for a computation with the Sun
    """
    approx_sidereal_drift = 5 * u.min


class MoonFlag(SpecialObjectFlag):
    """
    Flag for a computation with the Moon
    """
    approx_sidereal_drift = 60 * u.min
