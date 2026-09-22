# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Standard library
from abc import ABCMeta

# Third-party
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import (
    SkyCoord,
    ICRS,
    UnitSphericalRepresentation,
    AltAz,
    EarthLocation,
)

__all__ = ["Target", "FixedTarget", "AltAzTarget", "NonFixedTarget"]


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
    def is_time_dependent(self):
        """
        Whether this target requires evaluation at a specific time.

        Returns
        -------
        is_time_dependent : bool
            `True` for targets whose coordinates depend on ``obstime``,
            otherwise `False`.
        """
        return False

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

    @classmethod
    def _from_name_mock(cls, query_name, name=None):
        """
        Mock method to replace `FixedTarget.from_name` in tests without
        internet connection.
        """
        # The lowercase method will be run on names, so enter keys in lowercase:
        stars = {
            "rigel": {"ra": 78.63446707*u.deg, "dec": -8.20163837*u.deg},
            "sirius": {"ra": 101.28715533*u.deg, "dec": -16.71611586*u.deg},
            "vega": {"ra": 279.23473479*u.deg, "dec": 38.78368896*u.deg},
            "aldebaran": {"ra": 68.98016279*u.deg, "dec": 16.50930235*u.deg},
            "polaris": {"ra": 37.95456067*u.deg, "dec": 89.26410897*u.deg},
            "deneb": {"ra": 310.35797975*u.deg, "dec": 45.28033881*u.deg},
            "m13": {"ra": 250.423475*u.deg, "dec": 36.4613194*u.deg},
            "altair": {"ra": 297.6958273*u.deg, "dec": 8.8683212*u.deg},
            "hd 209458": {"ra": 330.79*u.deg, "dec": 18.88*u.deg}
        }

        if query_name.lower() in stars:
            return cls(coord=SkyCoord(**stars[query_name.lower()]),
                       name=query_name)
        else:
            raise ValueError("Target named {} not in mocked FixedTarget "
                             "method".format(query_name))


class AltAzTarget(Target):
    """
    Coordinates and metadata for a target defined in horizontal coordinates
    (altitude/azimuth) at a fixed observing location.

    Unlike `~astroplan.FixedTarget`, an `~astroplan.AltAzTarget` is *time-dependent*:
    the stored AltAz direction is evaluated at requested time(s) by transforming
    to ICRS, yielding an ICRS coordinate that varies with ``obstime``.

    This class is useful for targets defined by a local pointing direction (e.g.,
    “look at az=120°, alt=30° from this observatory”), rather than a fixed celestial
    coordinate.

    Notes
    -----
    The stored direction can be interpreted as geometric (vacuum) or apparent
    (refracted) depending on the atmospheric parameters provided. A pressure of
    ``0 hPa`` disables refraction.

    Downstream computations that transform the evaluated coordinate back to AltAz
    use the atmospheric parameters on the `~astroplan.Observer`. To preserve the
    exact apparent direction implied by this target's atmospheric parameters, use
    matching atmospheric parameters on the `~astroplan.Observer`.

    Examples
    --------
    Define a fixed horizontal direction at a given observatory:

    >>> import astropy.units as u
    >>> from astropy.coordinates import EarthLocation
    >>> from astroplan import AltAzTarget
    >>> location = EarthLocation.from_geodetic(-155.4761*u.deg, 19.825*u.deg,
    ...                                        4139*u.m)
    >>> t = AltAzTarget(alt=30*u.deg, az=120*u.deg, location=location, name="Pointing")
    """

    @property
    def is_time_dependent(self):
        """
        Whether this target requires evaluation at a specific time.
        """
        return True

    @u.quantity_input(alt=u.deg, az=u.deg)
    def __init__(
        self,
        alt,
        az,
        location,
        name=None,
        pressure=0 * u.hPa,
        temperature=None,
        relative_humidity=None,
        obswl=None,
        marker=None,
        **kwargs,
    ):
        """
        Parameters
        ----------
        alt : `~astropy.units.Quantity`
            Altitude angle. Must have angular units (e.g., ``u.deg``).

        az : `~astropy.units.Quantity`
            Azimuth angle. Must have angular units (e.g., ``u.deg``). By convention,
            azimuth is measured East of North.

        location : `~astropy.coordinates.EarthLocation`
            The observing location to which these AltAz coordinates apply.

        name : str, optional
            Name of the target, used for plotting and representing the target
            as a string.

        pressure : `~astropy.units.Quantity`, optional
            Atmospheric pressure used to interpret the stored AltAz direction.
            If set to ``0 hPa`` (default), the direction is treated as vacuum
            (geometric). If non-zero, the direction is treated as apparent
            (refracted) under the supplied atmospheric conditions.

        temperature : `~astropy.units.Quantity`, optional
            Ambient temperature for the refraction model (used when ``pressure``
            is non-zero). Default is ``0 deg_C``.

        relative_humidity : float, optional
            Relative humidity for the refraction model (used when ``pressure`` is
            non-zero). Must be in the interval [0, 1]. Default is 0.

        obswl : `~astropy.units.Quantity`, optional
            Observation wavelength for the refraction model (used when ``pressure``
            is non-zero). Default is ``1 micron``.

        marker : str, optional
            User-defined marker to differentiate between different types of targets
            (e.g., guides, high-priority, etc.).
        """
        if not isinstance(location, EarthLocation):
            raise TypeError("`location` must be an `astropy.coordinates.EarthLocation`.")

        self.name = name
        self.marker = marker

        self.alt = u.Quantity(alt).to(u.deg)
        self.az = u.Quantity(az).to(u.deg)
        self.location = location

        # Store atmosphere parameters for interpreting the stored AltAz direction.
        self.pressure = pressure
        self.temperature = temperature
        self.relative_humidity = relative_humidity
        self.obswl = obswl

    @classmethod
    def from_observer(cls, *, alt, az, observer, obswl=None, **kwargs):
        """
        Initialize an `~astroplan.AltAzTarget` from an `~astroplan.Observer`.

        This is a convenience constructor that uses the observer's location and
        atmospheric parameters to interpret the supplied AltAz direction.

        Parameters
        ----------
        alt : `~astropy.units.Quantity`
            Altitude angle.

        az : `~astropy.units.Quantity`
            Azimuth angle. By convention, azimuth is measured East of North.

        observer : `~astroplan.Observer`
            Observer that provides the location (and atmospheric parameters if
            present).

        obswl : `~astropy.units.Quantity`, optional
            Observation wavelength for the refraction model (used when ``pressure``
            is non-zero). Default is ``1 micron``.

        **kwargs
            Additional keywords passed to `~astroplan.AltAzTarget` (e.g., ``name``,
            ``marker``).

        Returns
        -------
        target : `~astroplan.AltAzTarget`
            The constructed target.
        """
        return cls(
            alt=alt, az=az, location=observer.location,
            pressure=observer.pressure, temperature=observer.temperature,
            relative_humidity=observer.relative_humidity, obswl=obswl,
            **kwargs
        )

    def __repr__(self):
        class_name = self.__class__.__name__
        alt = self.alt.to(u.deg).value
        az = self.az.to(u.deg).value
        return '<{} "{}" at (alt, az)=({:.6f} deg, {:.6f} deg)>'.format(
            class_name, self.name, alt, az
        )

    def get_skycoord(self, times, observer=None):
        """
        Evaluate this target to an ICRS `~astropy.coordinates.SkyCoord` at ``times``.

        Parameters
        ----------
        times : `~astropy.time.Time` or time-like
            Times at which to evaluate the target.

        observer : `~astroplan.Observer`, optional
            Observer supplied through the common `get_skycoord()` interface.
            This target's AltAz direction is defined by its stored location,
            so this argument is not used.

        Returns
        -------
        coord : `~astropy.coordinates.SkyCoord`
            ICRS coordinate evaluated at ``times`` (time-dependent).
        """
        if times is None:
            raise ValueError("`times` is required to evaluate an AltAzTarget.")
        if not isinstance(times, Time):
            times = Time(times)

        # Construct the AltAz frame used to interpret the stored alt/az direction.
        altaz_frame = AltAz(
            location=self.location,
            obstime=times,
            pressure=self.pressure,
            temperature=self.temperature,
            relative_humidity=self.relative_humidity,
            obswl=self.obswl,
        )

        # The stored alt/az are treated as scalar directions and broadcast to `times`.
        alt = u.Quantity(np.broadcast_to(self.alt.to_value(u.deg), times.shape), u.deg)
        az = u.Quantity(np.broadcast_to(self.az.to_value(u.deg), times.shape), u.deg)

        return SkyCoord(az=az, alt=alt, frame=altaz_frame).icrs


class NonFixedTarget(Target):
    """
    Placeholder for future function.
    """


def get_skycoord(targets, times=None, observer=None):
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

    observer : `~astroplan.Observer`, optional
        Observer to use when evaluating targets whose apparent sky position
        depends on the observing location. Ignored for coordinates and targets
        that do not require observer context.

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
        return isinstance(obj, Target) and obj.is_time_dependent

    def _as_coord(obj):
        if hasattr(obj, "coord"):
            return obj.coord
        if callable(getattr(obj, "get_skycoord", None)):
            return obj.get_skycoord(times, observer=observer)
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
