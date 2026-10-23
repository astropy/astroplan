# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Standard library
from collections.abc import Mapping
from datetime import datetime
from functools import cached_property
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
    ITRS,
    SkyCoord,
    TEME,
    UnitSphericalRepresentation,
)

# Package
from .exceptions import SatellitePropagationWarning
from .target import Target


__all__ = ["SGP4SatelliteTarget"]


class SGP4SatelliteTarget(Target):
    """
    A satellite target propagated from TLE or OMM elements using SGP4.

    Notes
    -----
    GP element sets are generally useful only near their epoch. Scheduling
    fast-moving satellites requires a time resolution appropriate to the
    orbit; for a LEO target ``time_resolution=10*u.s`` can be appropriate.

    Satellite targets may rise and set multiple times within a day.
    Rise and set searches return a single event and may require increasing
    ``n_grid_points`` to resolve short passes. Meridian-transit searches
    locate altitude extrema, which do not necessarily coincide with meridian
    crossings for fast-moving satellites. Their results should not be
    interpreted as true meridian-transit times.

    The propagated TEME position is observer-independent. `get_skycoord`
    returns a topocentric ICRS direction for the provided observer without
    atmospheric refraction. No light-time correction is applied; the
    satellite position is evaluated at the requested time rather than
    at the time the observed light was emitted.
    `~astroplan.Observer.altaz` uses the observer's atmospheric parameters
    when transforming this direction.
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
            lines = [line.strip() for line in tle.splitlines() if line.strip()]
            if len(lines) not in (2, 3):
                raise ValueError(
                    f"Expected TLE string to contain 2 or 3 lines, got {len(lines)}."
                )
            if len(lines) == 3:
                name, lines = lines[0].removeprefix("0 ").strip(), lines[1:]
        else:
            try:
                lines = list(tle)
            except TypeError:
                raise TypeError(
                    "`tle` must be a two- or three-line string or a two-element sequence."
                ) from None
            if len(lines) != 2 or not all(isinstance(line, str) for line in lines):
                raise ValueError("TLE sequences must contain exactly two strings.")
            lines = [line.strip() for line in lines]
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
        # OMM may use YYYY-DDDTHH:MM:SS, while Astropy's yday format uses YYYY:DDD:HH:MM:SS
        match = re.fullmatch(r"(\d{4}-\d{3})(T.+)", value)
        if match:
            try:
                date = datetime.strptime(match.group(1), "%Y-%j")
                # Leap-year check
                if date.year != int(match.group(1)[:4]):
                    raise ValueError
            except ValueError:
                raise ValueError("Unsupported OMM `EPOCH`.") from None
            value = date.strftime("%Y-%m-%d") + match.group(2)

        try:
            return Time(value, scale="utc", precision=6).isot
        except ValueError:
            raise ValueError("Unsupported OMM `EPOCH`.") from None

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

    @cached_property
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

        age_days = np.abs((utc - self.epoch).to_value(u.day))
        if np.any((age_days > 14) & ~failed.reshape(shape)):
            warnings.warn(
                "SGP4 elements are more than 14 days from their epoch; "
                "propagated positions may be inaccurate.",
                SatellitePropagationWarning,
            )

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
                flat_times = utc.ravel()
                for code in np.unique(errors[failed]):
                    code_times = flat_times[errors == code]
                    first = code_times[np.argmin(code_times.jd)].isot
                    last = code_times[np.argmax(code_times.jd)].isot
                    message = SGP4_ERRORS.get(int(code), "unknown error")
                    details.append(
                        f"code {code} ({message}): {code_times.size} times "
                        f"between {first} and {last}"
                    )
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

        teme = self.get_teme(times)
        # Times where SGP4 propagation failed have NaN positions (get_teme has
        # already warned about them); don't let the transforms warn again.
        with np.errstate(invalid="ignore"):
            # Intermediate topocentric ITRS step to avoid change in stellar aberration
            geocentric_itrs = teme.realize_frame(
                teme.cartesian.without_differentials()
            ).transform_to(ITRS(obstime=times))
            observer_itrs = observer.location.get_itrs(obstime=times)
            topocentric_position = (
                geocentric_itrs.cartesian - observer_itrs.cartesian
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
