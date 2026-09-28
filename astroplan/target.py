# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Standard library
from abc import ABCMeta
from collections.abc import Mapping
from datetime import datetime, UTC
import csv
import io
import json
import re
from typing import ClassVar
import warnings

# Third-party
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import (
    CartesianDifferential,
    CartesianRepresentation,
    ICRS,
    ITRS,
    SkyCoord,
    TEME,
    UnitSphericalRepresentation,
)

# Package
from .exceptions import SatellitePropagationWarning


__all__ = ["Target", "FixedTarget", "SGP4SatelliteTarget", "NonFixedTarget"]

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


class NonFixedTarget(Target):
    """
    Placeholder for future function.
    """


class SGP4SatelliteTarget(Target):
    """
    A satellite target propagated from TLE or OMM elements using SGP4.

    Notes
    -----
    GP element sets are generally useful only near their epoch. Scheduling
    fast-moving satellites requires a time resolution appropriate to the
    orbit; for a LEO target ``time_resolution=10*u.s`` can be appropriate.
    The default arguments for rise, set, and meridian-transit searches,
    such as `n_grid_points` in `target_rise_time`, might not be precise enough
    for finding LEO passes.
    """

    _OMM_FORMATS: ClassVar[set[str]] = {"kvn", "xml", "json", "csv"}
    _OMM_CORE_FIELDS: ClassVar[set[str]] = {
        "EPOCH",
        "MEAN_MOTION",
        "ECCENTRICITY",
        "INCLINATION",
        "RA_OF_ASC_NODE",
        "ARG_OF_PERICENTER",
        "MEAN_ANOMALY",
        "BSTAR",
    }
    _OMM_METADATA: ClassVar[dict[str, str]] = {
        "CENTER_NAME": "EARTH",
        "REF_FRAME": "TEME",
        "TIME_SYSTEM": "UTC",
        "MEAN_ELEMENT_THEORY": "SGP4",
    }
    _OMM_DEFAULTS: ClassVar[dict[str, str]] = {
        "CLASSIFICATION_TYPE": "U",
        "OBJECT_ID": "",
        "EPHEMERIS_TYPE": "0",
        "ELEMENT_SET_NO": "0",
        "REV_AT_EPOCH": "0",
        "MEAN_MOTION_DOT": "0.0",
        "MEAN_MOTION_DDOT": "0.0",
        "NORAD_CAT_ID": "0",
    }
    _OMM_UNITS: ClassVar[dict[str, str]] = {
        "MEAN_MOTION": "rev/day",
        "INCLINATION": "deg",
        "RA_OF_ASC_NODE": "deg",
        "ARG_OF_PERICENTER": "deg",
        "MEAN_ANOMALY": "deg",
        "BSTAR": "1/ER",
        "MEAN_MOTION_DOT": "rev/day**2",
        "MEAN_MOTION_DDOT": "rev/day**3",
    }
    _OMM_USED_FIELDS: ClassVar[set[str]] = _OMM_CORE_FIELDS | set(_OMM_DEFAULTS)

    @property
    def is_time_dependent(self):
        return True

    def __init__(
        self,
        *,
        tle=None,
        omm=None,
        omm_format=None,
        name=None,
        gravity_model="wgs72",
        validate=True,
    ):
        """
        Parameters
        ----------
        tle : str or sequence of str (optional)
            Two-line TLE string, three-line TLE string including a name, or a
            two-element sequence containing the TLE element lines. For sequence
            input, the target name can be provided separately with name.
        omm : mapping, str, bytes, or file-like (optional)
            One OMM record. KVN, XML, JSON, and CSV text are supported.
        omm_format : {None, 'kvn', 'xml', 'json', 'csv'} (optional)
            OMM format for text, bytes, or file-like input. If `None`, the code tries
            to infer the format from the input. Ignored for mapping input.
        name : str (optional)
            Name of the target. This overrides a TLE line-zero name or OMM
            ``OBJECT_NAME``.
        gravity_model : {'wgs72', 'wgs72old', 'wgs84'} (optional)
            SGP4 gravity model. The default is WGS72, which normally gives the
            best agreement with standard GP element products.
        validate : bool (optional)
            Whether to strictly validate TLE fixed-width formatting.
        """
        if (tle is None) == (omm is None):
            raise ValueError("Provide exactly one of `tle` or `omm`.")
        try:
            import sgp4  # noqa: F401
        except ModuleNotFoundError as exc:
            raise ImportError(
                "Install astroplan with the `satellite` extra to use SGP4SatelliteTarget."
            ) from exc

        gravity_constant, validation_gravity_constant = self._gravity_constants(
            gravity_model
        )
        element_name = None
        if tle is not None:
            line1, line2, element_name = self._normalize_tle(tle)
            self._satrec = self._build_satrec_from_tle(
                line1, line2, gravity_constant, validation_gravity_constant, validate
            )
            self._source_format = "tle"
        else:
            fields, resolved_omm_format = self._parse_omm_record(omm, omm_format)
            fields = self._normalize_omm_fields(fields, resolved_omm_format)
            self._satrec = self._build_satrec_from_omm(fields, gravity_constant)
            self._source_format = "omm"
            element_name = fields.get("OBJECT_NAME")

        if name is not None:
            self.name = str(name)
        elif element_name:
            self.name = str(element_name)
        elif self.catalog_number:
            self.name = f"Satellite {self.catalog_number}"
        else:
            self.name = "Unnamed satellite"

    @staticmethod
    def _gravity_constants(gravity_model):
        """Return SGP4 gravity constants for propagation and TLE validation."""
        from sgp4.api import WGS72, WGS72OLD, WGS84
        from sgp4.earth_gravity import wgs72, wgs72old, wgs84

        models = {
            "wgs72": (WGS72, wgs72),
            "wgs72old": (WGS72OLD, wgs72old),
            "wgs84": (WGS84, wgs84),
        }
        try:
            return models[gravity_model.lower()]
        except (AttributeError, KeyError):
            raise ValueError(
                "`gravity_model` must be 'wgs72', 'wgs72old', or 'wgs84'."
            ) from None

    @staticmethod
    def _normalize_tle(tle):
        """Normalize TLE input to two element lines and an optional name"""
        name = None
        if isinstance(tle, str):
            lines = [line for line in tle.splitlines() if line.strip()]
            if len(lines) not in (2, 3):
                raise ValueError(
                    f"Expected TLE string to contain 2 or 3 lines, got {len(lines)}."
                )
            if len(lines) == 3:
                name, lines = lines[0].strip(), lines[1:]
        else:
            try:
                lines = list(tle)
            except TypeError:
                raise TypeError(
                    "`tle` must be a two- or three-line string or a two-element sequence."
                ) from None
            if len(lines) != 2 or not all(isinstance(line, str) for line in lines):
                raise ValueError("TLE sequences must contain exactly two strings.")
        line1, line2 = lines
        if not line1.startswith("1 ") or not line2.startswith("2 "):
            raise ValueError("TLE element lines must begin with '1 ' and '2 '.")
        return line1, line2, name

    @staticmethod
    def _is_kvn_keyword(keyword):
        """Return whether a string is a valid OMM KVN keyword."""
        return re.fullmatch(r"[A-Z][A-Z0-9_]*", keyword) is not None

    @classmethod
    def _detect_omm_format(cls, text):
        """Detect the OMM serialization format from its text."""
        # Remove UTF-8 Unicode BOM (byte-order mark)
        text = text.removeprefix("\ufeff").lstrip()
        if text.startswith("<"):
            return "xml"
        if text.startswith(("{", "[")):
            return "json"

        first_line = next((line for line in text.splitlines() if line.strip()), "")
        key, separator, _ = first_line.partition("=")
        if separator and cls._is_kvn_keyword(key.strip()):
            return "kvn"

        try:
            header = {field.strip().upper() for field in next(csv.reader([first_line]))}
        except (csv.Error, StopIteration):
            header = set()
        if cls._OMM_CORE_FIELDS <= header:
            return "csv"
        raise ValueError(
            "Could not determine OMM format; specify `omm_format` explicitly."
        )

    @classmethod
    def _parse_omm_record(cls, value, omm_format=None):
        """Parse one OMM record and return its fields and resolved format."""
        if omm_format is not None:
            try:
                omm_format = omm_format.lower()
            except AttributeError:
                raise ValueError(
                    "`omm_format` must be 'kvn', 'xml', 'json', 'csv', or None."
                ) from None
            if omm_format not in cls._OMM_FORMATS:
                raise ValueError(
                    "`omm_format` must be 'kvn', 'xml', 'json', 'csv', or None."
                )
        if isinstance(value, Mapping):
            return dict(value), "mapping"
        read = getattr(value, "read", None)
        if callable(read):
            value = read()
        if isinstance(value, bytes):
            value = value.decode("utf-8")
        if not isinstance(value, str):
            raise TypeError(
                "`omm` must be a mapping, text, UTF-8 bytes, or a file-like object."
            )
        # Remove UTF-8 Unicode BOM (byte-order mark) and leading whitespace
        text = value.removeprefix("\ufeff").lstrip()
        omm_format = omm_format or cls._detect_omm_format(text)
        if omm_format == "json":
            records = json.loads(text)
            if isinstance(records, Mapping):
                return dict(records), "json"
            if not isinstance(records, list):
                raise ValueError(
                    "Expected a JSON object or a one-element array containing an object."
                )
            if len(records) != 1:
                raise ValueError(
                    f"Expected exactly one OMM record; found {len(records)}."
                )
            if not isinstance(records[0], Mapping):
                raise ValueError(
                    "Expected a JSON object or a one-element array containing an object."
                )
            return dict(records[0]), "json"
        if omm_format in ("csv", "xml"):
            from sgp4 import omm as sgp4_omm

            parser = sgp4_omm.parse_csv if omm_format == "csv" else sgp4_omm.parse_xml
            records = list(parser(io.StringIO(text)))
            if len(records) != 1:
                raise ValueError(
                    f"Expected exactly one OMM record; found {len(records)}."
                )
            return records[0], omm_format
        return cls._parse_omm_kvn(text), "kvn"

    @classmethod
    def _parse_omm_kvn(cls, text):
        """Parse one OMM record serialized in CCSDS KVN format."""
        fields = {}
        for line in text.splitlines():
            line = line.strip()
            if not line or line == "COMMENT" or line.startswith("COMMENT "):
                continue
            if "=" not in line:
                raise ValueError(f"Invalid OMM KVN line: {line!r}.")
            key, value = (part.strip() for part in line.split("=", 1))
            if not cls._is_kvn_keyword(key):
                raise ValueError(f"Invalid OMM KVN keyword: {key!r}.")
            if key in fields:
                raise ValueError(f"Duplicate OMM KVN keyword: {key!r}.")

            # Split an optional trailing `[unit]` suffix from the field value
            match = re.fullmatch(r"(.*?)\s*\[([^]]+)\]\s*", value)
            if match and key in cls._OMM_UNITS:
                value, unit = match.group(1).strip(), match.group(2).strip()
                if unit != cls._OMM_UNITS[key]:
                    raise ValueError(
                        f"Unsupported unit {unit!r} for OMM field `{key}`; expected "
                        f"{cls._OMM_UNITS[key]!r}."
                    )
            elif match and key in cls._OMM_USED_FIELDS:
                raise ValueError(f"Units are not supported for OMM field `{key}`.")
            fields[key] = value
        return fields

    @staticmethod
    def _normalize_epoch(epoch):
        """Normalize an OMM epoch to the format expected by python-sgp4."""
        value = str(epoch).strip().removesuffix("Z")
        for format_ in ("%Y-%m-%dT%H:%M:%S.%f", "%Y-%m-%dT%H:%M:%S"):
            try:
                parsed = datetime.strptime(value, format_).replace(tzinfo=UTC)
            except ValueError:
                continue
            return parsed.strftime("%Y-%m-%dT%H:%M:%S.%f")
        raise ValueError(
            "Unsupported OMM `EPOCH`; expected YYYY-MM-DDTHH:MM:SS[.ffffff]."
        )

    @classmethod
    def _normalize_omm_fields(cls, fields, omm_format):
        """Normalize OMM fields to the format expected by python-sgp4."""
        normalized = {}
        for key, value in fields.items():
            key = str(key).strip().upper()
            if key in normalized:
                raise ValueError(
                    f"Duplicate OMM field after key normalization: {key!r}."
                )
            normalized[key] = value
        fields = normalized

        # CelesTrak XML may use null for an unavailable international designator
        if fields.get("OBJECT_ID") is None:
            fields["OBJECT_ID"] = ""

        if omm_format in ("xml", "kvn"):
            missing = [key for key in cls._OMM_METADATA if key not in fields]
            if missing:
                raise ValueError(
                    "Missing required OMM metadata: {}.".format(", ".join(missing))
                )
        else:
            # CelesTrak JSON and CSV may omit mandatory metadata; fill in defaults
            for key, value in cls._OMM_METADATA.items():
                fields.setdefault(key, value)
        for key, expected in cls._OMM_METADATA.items():
            value = str(fields[key]).strip().upper()
            if key == "MEAN_ELEMENT_THEORY":
                if value not in ("SGP4", "SGP/SGP4"):
                    raise ValueError(f"Unsupported OMM `{key}`: {fields[key]!r}.")
            elif value != expected:
                raise ValueError(f"Unsupported OMM `{key}`: {fields[key]!r}.")
            fields[key] = value

        if "MEAN_MOTION" not in fields and "SEMI_MAJOR_AXIS" in fields:
            raise ValueError(
                "`SEMI_MAJOR_AXIS` is not supported in place of `MEAN_MOTION`."
            )
        missing = sorted(cls._OMM_CORE_FIELDS - fields.keys())
        if missing:
            raise ValueError(
                "Missing required OMM fields: {}.".format(", ".join(missing))
            )
        # Supply safe defaults when fields are absent
        for key, value in cls._OMM_DEFAULTS.items():
            fields.setdefault(key, value)

        try:
            ephemeris_type = int(fields["EPHEMERIS_TYPE"])
        except (TypeError, ValueError):
            raise ValueError("`EPHEMERIS_TYPE` must be an integer.") from None
        if ephemeris_type not in (0, 2):
            raise ValueError(f"Unsupported OMM `EPHEMERIS_TYPE`: {ephemeris_type}.")
        fields["EPHEMERIS_TYPE"] = str(ephemeris_type)
        fields["EPOCH"] = cls._normalize_epoch(fields["EPOCH"])

        values = {}
        for key in (
            "MEAN_MOTION",
            "ECCENTRICITY",
            "INCLINATION",
            "RA_OF_ASC_NODE",
            "ARG_OF_PERICENTER",
            "MEAN_ANOMALY",
            "BSTAR",
            "MEAN_MOTION_DOT",
            "MEAN_MOTION_DDOT",
        ):
            try:
                values[key] = float(fields[key])
            except (TypeError, ValueError):
                raise ValueError(f"OMM field `{key}` must be numeric.") from None
            if not np.isfinite(values[key]):
                raise ValueError(f"OMM field `{key}` must be finite.")
        if not 0 <= values["ECCENTRICITY"] < 1:
            raise ValueError("OMM `ECCENTRICITY` must satisfy 0 <= ECCENTRICITY < 1.")
        if values["MEAN_MOTION"] <= 0:
            raise ValueError("OMM `MEAN_MOTION` must be greater than zero.")
        try:
            catalog_number = int(fields["NORAD_CAT_ID"])
        except (TypeError, ValueError):
            raise ValueError("`NORAD_CAT_ID` must be an integer.") from None
        if catalog_number < 0:
            raise ValueError("`NORAD_CAT_ID` must be non-negative.")
        return fields

    @staticmethod
    def _build_satrec_from_tle(
        line1, line2, gravity_constant, validation_gravity_constant, validate=True
    ):
        """Initialize a `Satrec` from TLE element lines."""
        from sgp4.api import Satrec

        if validate:
            from sgp4.io import twoline2rv

            twoline2rv(line1, line2, validation_gravity_constant)
        return Satrec.twoline2rv(line1, line2, gravity_constant)

    @staticmethod
    def _build_satrec_from_omm(fields, gravity_constant):
        """Initialize and validate a `Satrec` from normalized OMM fields."""
        from sgp4 import omm
        from sgp4.api import Satrec
        from sgp4.conveniences import check_satrec

        satrec = Satrec()
        omm.initialize(satrec, fields, gravity_constant)
        check_satrec(satrec)
        return satrec

    @property
    def epoch(self):
        """Element epoch as an Astropy `~astropy.time.Time`."""
        return Time(
            self._satrec.jdsatepoch, self._satrec.jdsatepochF, format="jd", scale="utc"
        )

    @property
    def catalog_number(self):
        """Integer catalog number, or zero if an OMM record did not supply one."""
        return int(self._satrec.satnum)

    @property
    def source_format(self):
        """Element source format, either ``'tle'`` or ``'omm'``."""
        return self._source_format

    def to_tle(self):
        """
        Return the orbital elements as a normalized two-line element set.

        Returns
        -------
        line1, line2 : tuple of str
            Canonically formatted TLE element lines including checksums.
        """
        from sgp4.exporter import export_tle
        return export_tle(self._satrec)

    def to_omm(self):
        """
        Return the orbital elements as a normalized OMM record.

        Returns
        -------
        fields : dict
            OMM fields using the values represented by the internal SGP4 record.
        """
        from sgp4.exporter import export_omm
        return export_omm(self._satrec, self.name)

    def get_teme(self, times):
        """
        Return the observer-independent geocentric “True Equator Mean Equinox” (TEME) state.

        Parameters
        ----------
        times : `~astropy.time.Time` or time-like
            Time(s) at which to propagate the satellite.

        Returns
        -------
        teme : `~astropy.coordinates.TEME`
            TEME position and velocity at the requested times.
        """
        if times is None:
            raise ValueError("`times` is required to evaluate an SGP4SatelliteTarget.")
        if not isinstance(times, Time):
            times = Time(times)
        utc = times.utc
        shape = utc.shape
        jd1 = np.atleast_1d(utc.jd1).ravel()
        jd2 = np.atleast_1d(utc.jd2).ravel()
        errors, positions, velocities = self._satrec.sgp4_array(jd1, jd2)
        errors = np.asarray(errors)
        failed = errors != 0
        if np.any(failed):
            from sgp4.api import SGP4_ERRORS

            if utc.isscalar:
                message = SGP4_ERRORS.get(
                    int(errors[0]), f"unknown error code {errors[0]}"
                )
                warnings.warn(
                    f"SGP4 propagation failed: {message}", SatellitePropagationWarning
                )
            else:
                details = []
                for code in np.unique(errors[failed]):
                    indices = ", ".join(str(i) for i in np.flatnonzero(errors == code))
                    message = SGP4_ERRORS.get(int(code), "unknown error")
                    details.append(f"code {code} ({message}): indices {indices}")
                warnings.warn(
                    f"SGP4 propagation failed for {failed.sum()} of {errors.size} times:\n"
                    + "\n".join(details),
                    SatellitePropagationWarning,
                )
            positions = np.asarray(positions).copy()
            velocities = np.asarray(velocities).copy()
            positions[failed] = np.nan
            velocities[failed] = np.nan

        positions = np.asarray(positions).reshape(shape + (3,))
        velocities = np.asarray(velocities).reshape(shape + (3,))
        position = CartesianRepresentation(np.moveaxis(positions, -1, 0) * u.km)
        velocity = CartesianDifferential(np.moveaxis(velocities, -1, 0) * u.km / u.s)
        return TEME(position.with_differentials(velocity), obstime=utc)

    def get_skycoord(self, times, observer=None):
        """
        Return the observer-corrected ICRS coordinate.

        Parameters
        ----------
        times : `~astropy.time.Time` or time-like
            Time(s) at which to evaluate the target.
        observer : `~astroplan.Observer`
            Observer from which to evaluate the target.

        Returns
        -------
        coord : `~astropy.coordinates.SkyCoord`
            Topocentric direction transformed to ICRS.
        """
        # https://docs.astropy.org/en/stable/coordinates/satellites.html
        if times is None:
            raise ValueError("`times` is required to evaluate an SGP4SatelliteTarget.")
        if observer is None:
            raise ValueError(
                "`observer` is required to evaluate an SGP4SatelliteTarget."
            )
        if not isinstance(times, Time):
            times = Time(times)
        times = times.utc

        # Intermediate topocentric ITRS step to avoid change in stellar aberration
        geocentric_itrs = self.get_teme(times).transform_to(ITRS(obstime=times))
        observer_itrs = observer.location.get_itrs(obstime=times)
        topocentric_position = (
            geocentric_itrs.cartesian.without_differentials() - observer_itrs.cartesian
        )
        topocentric_itrs = ITRS(
            topocentric_position, obstime=times, location=observer.location
        )
        topocentric_direction = topocentric_itrs.realize_frame(
            topocentric_itrs.represent_as(UnitSphericalRepresentation)
        )
        return SkyCoord(topocentric_direction).icrs

    def __repr__(self):
        return (
            f'<{self.__class__.__name__} "{self.name}" catalog #{self.catalog_number}>'
        )

    def __str__(self):
        return self.name


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
