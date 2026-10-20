# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Standard library
import csv
import io
import json

# Third-party
import astropy.units as u
import pytest
from astropy.coordinates import AltAz, EarthLocation, GCRS, ICRS, ITRS, SkyCoord, TEME
from astropy.time import Time
import numpy as np

try:
    import sgp4  # noqa: F401

    HAS_SGP4 = True
except ImportError:
    HAS_SGP4 = False

# Package
from astroplan.exceptions import AstroplanWarning, SatellitePropagationWarning
from astroplan.target import FixedTarget, AltAzTarget, SGP4SatelliteTarget, get_skycoord
from astroplan.observer import Observer
from astroplan.utils import time_grid_from_range


ISS_TLE = (
    "1 25544U 98067A   23215.27256123  .00041610  00000-0  73103-3 0  9990",
    "2 25544  51.6403  95.2411 0000623 157.9606 345.0624 15.50085581409092",
)
ISS_OMM = {
    "OBJECT_NAME": "ISS (ZARYA)",
    "OBJECT_ID": "1998-067A",
    "CENTER_NAME": "EARTH",
    "REF_FRAME": "TEME",
    "TIME_SYSTEM": "UTC",
    "MEAN_ELEMENT_THEORY": "SGP4",
    "EPOCH": "2023-08-03T06:32:29.290272",
    "MEAN_MOTION": "15.50085581",
    "ECCENTRICITY": "0.0000623",
    "INCLINATION": "51.6403",
    "RA_OF_ASC_NODE": "95.2411",
    "ARG_OF_PERICENTER": "157.9606",
    "MEAN_ANOMALY": "345.0624",
    "EPHEMERIS_TYPE": "0",
    "CLASSIFICATION_TYPE": "U",
    "NORAD_CAT_ID": "25544",
    "ELEMENT_SET_NO": "999",
    "REV_AT_EPOCH": "40909",
    "BSTAR": "0.00073103",
    "MEAN_MOTION_DOT": "0.00041610",
    "MEAN_MOTION_DDOT": "0.0",
}
OMM_METADATA = {"CENTER_NAME", "REF_FRAME", "TIME_SYSTEM", "MEAN_ELEMENT_THEORY"}
ISS_OMM_GP = {key: value for key, value in ISS_OMM.items() if key not in OMM_METADATA}


def _omm_csv(records=(ISS_OMM_GP,)):
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=records[0])
    writer.writeheader()
    writer.writerows(records)
    return stream.getvalue()


def _omm_xml(records=(ISS_OMM,)):
    metadata = (
        "OBJECT_NAME",
        "OBJECT_ID",
        "CENTER_NAME",
        "REF_FRAME",
        "TIME_SYSTEM",
        "MEAN_ELEMENT_THEORY",
    )
    mean = (
        "EPOCH",
        "MEAN_MOTION",
        "ECCENTRICITY",
        "INCLINATION",
        "RA_OF_ASC_NODE",
        "ARG_OF_PERICENTER",
        "MEAN_ANOMALY",
    )
    tle = (
        "EPHEMERIS_TYPE",
        "CLASSIFICATION_TYPE",
        "NORAD_CAT_ID",
        "ELEMENT_SET_NO",
        "REV_AT_EPOCH",
        "BSTAR",
        "MEAN_MOTION_DOT",
        "MEAN_MOTION_DDOT",
    )

    def elements(record, keys):
        return "".join(f"<{key}>{record[key]}</{key}>" for key in keys if key in record)

    segments = [
        f"<segment><metadata>{elements(record, metadata)}</metadata>"
        f"<data><meanElements>{elements(record, mean)}</meanElements>"
        f"<tleParameters>{elements(record, tle)}</tleParameters></data></segment>"
        for record in records
    ]
    return "<ndm><body>{}</body></ndm>".format("".join(segments))


def _omm_kvn(fields=ISS_OMM):
    units = {
        "MEAN_MOTION": "rev/day",
        "INCLINATION": "deg",
        "RA_OF_ASC_NODE": "deg",
        "ARG_OF_PERICENTER": "deg",
        "MEAN_ANOMALY": "deg",
        "BSTAR": "1/ER",
        "MEAN_MOTION_DOT": "rev/day**2",
        "MEAN_MOTION_DDOT": "rev/day**3",
    }
    lines = ["CCSDS_OMM_VERS = 2.0", "COMMENT Offline test fixture"]
    lines += [
        f"{key} = {value}" + (f" [{units[key]}]" if key in units else "")
        for key, value in fields.items()
    ]
    return "\n".join(lines)


@pytest.mark.remote_data
def test_FixedTarget_from_name():
    """
    Check that resolving target names with the `SkyCoord.from_name` constructor
    to produce a `FixedTarget` accurately resolves the coordinates of Polaris.
    """

    # Resolve coordinates with SkyCoord.from_name classmethod
    polaris_from_name = FixedTarget.from_name('Polaris')
    polaris_from_name = FixedTarget.from_name('Polaris', name='Target 1')
    # Coordinates grabbed from SIMBAD
    polaris_from_SIMBAD = SkyCoord('02h31m49.09456s', '+89d15m50.7923s')

    # Make sure separation is small
    assert polaris_from_name.coord.separation(polaris_from_SIMBAD) < 1*u.arcsec


@pytest.mark.remote_data
def test_FixedTarget_ra_dec():
    """
    Confirm that FixedTarget.ra and FixedTarget.dec are the same as the
    right ascension and declination stored in the FixedTarget.coord variable -
    which is a SkyCoord
    """

    vega_coords = SkyCoord('18h36m56.33635s', '+38d47m01.2802s')
    vega = FixedTarget(vega_coords, name='Vega')
    assert vega.coord == vega_coords, 'Store coordinates directly'
    assert vega.coord.ra == vega_coords.ra == vega.ra, ('Retrieve RA from '
                                                        'SkyCoord')
    assert vega.coord.dec == vega_coords.dec == vega.dec, ('Retrieve Dec from '
                                                           'SkyCoord')


def test_AltAzTarget_quantity_validation():
    """
    `AltAzTarget` should require angle quantities for `alt` and `az`, and a valid
    EarthLocation.
    """
    location = EarthLocation.from_geodetic(lon=0 * u.deg, lat=0 * u.deg, height=0 * u.m)

    with pytest.raises(TypeError):
        AltAzTarget(alt=30, az=120, location=location)

    with pytest.raises(TypeError):
        AltAzTarget(alt=30 * u.deg, az=120 * u.deg, location="not a location")

    with pytest.raises(TypeError, match="temprature"):
        AltAzTarget(alt=30*u.deg, az=120*u.deg, location=location,
                    temprature=20*u.deg_C)


def test_AltAzTarget_get_skycoord_requires_times():
    """
    `AltAzTarget.get_skycoord` and `get_skycoord` should require `times` when
    evaluating time-dependent targets.
    """
    location = EarthLocation.from_geodetic(lon=0 * u.deg, lat=0 * u.deg, height=0 * u.m)
    target = AltAzTarget(alt=45 * u.deg, az=0 * u.deg, location=location)

    with pytest.raises(ValueError):
        target.get_skycoord(None)

    with pytest.raises(ValueError):
        get_skycoord([target])


def test_AltAzTarget_from_observer():
    """
    `AltAzTarget.from_observer` should construct an AltAzTarget that inherits the
    observer's location and atmospheric parameters, and evaluates identically to
    a manually-constructed AltAzTarget with the same inputs.
    """
    location = EarthLocation.from_geodetic(
        lon=10 * u.deg, lat=45 * u.deg, height=100 * u.m
    )
    observer = Observer(
        location=location,
        pressure=800 * u.hPa,
        temperature=10 * u.deg_C,
        relative_humidity=0.25,
        timezone="UTC",
        name="Test Observer",
    )

    target = AltAzTarget.from_observer(
        alt=30 * u.deg,
        az=120 * u.deg,
        observer=observer,
        name="AltAz via Observer",
        marker="test",
        obswl=2 * u.micron,
    )

    assert target.location == observer.location
    assert target.pressure == observer.pressure
    assert target.temperature == observer.temperature
    assert target.relative_humidity == observer.relative_humidity
    assert target.name == "AltAz via Observer"
    assert target.marker == "test"

    manual = AltAzTarget(
        alt=30 * u.deg,
        az=120 * u.deg,
        location=observer.location,
        pressure=observer.pressure,
        temperature=observer.temperature,
        relative_humidity=observer.relative_humidity,
        obswl=2 * u.micron,
    )

    times = Time(
        ["2026-02-05 00:00", "2026-02-05 06:00", "2026-02-05 12:00"]
    )

    coord_from_observer = target.get_skycoord(times)
    coord_manual = manual.get_skycoord(times)

    # They should be effectively identical
    assert coord_from_observer.separation(coord_manual).max() < 1e-6 * u.arcsec

    default_observer = Observer(location=location)
    default_target = AltAzTarget.from_observer(alt=30*u.deg, az=120*u.deg,
                                               observer=default_observer)
    assert default_target.pressure == 0*u.hPa


def test_AltAzTarget_warns_for_different_observer():
    location = EarthLocation.from_geodetic(10*u.deg, 45*u.deg, 0*u.m)
    other_location = EarthLocation.from_geodetic(20*u.deg, 45*u.deg, 0*u.m)
    target = AltAzTarget(alt=30*u.deg, az=120*u.deg, location=location)
    observer = Observer(location=other_location)
    time = Time("2026-02-05T00:00:00", scale="utc")

    with pytest.warns(AstroplanWarning, match="different location"):
        target.get_skycoord(time, observer=observer)


def test_AltAzTarget_get_skycoord_vector_times_shape_and_frame():
    """
    Evaluating an `AltAzTarget` at vector times should return an ICRS SkyCoord
    with shape matching `times.shape` and values that vary with time.
    """
    location = EarthLocation.from_geodetic(lon=0 * u.deg, lat=0 * u.deg, height=0 * u.m)
    target = AltAzTarget(
        alt=45 * u.deg,
        az=0 * u.deg,
        location=location,
        pressure=None,
        temperature=None,
        relative_humidity=None,
        obswl=None,
    )

    t0 = Time("2026-02-05 00:00")
    times = t0 + np.array([0, 3, 6]) * u.hour

    coord = target.get_skycoord(times)

    assert coord.is_equivalent_frame(ICRS())
    assert coord.shape == times.shape
    assert coord.size == times.size

    # Coordinate should change with time
    assert coord[0].separation(coord[-1]) > 50 * u.deg


def test_AltAzTarget_apparent_vs_vacuum_differ_when_pressure_nonzero():
    location = EarthLocation.from_geodetic(10*u.deg, 45*u.deg, 0*u.m)

    t = Time("2026-01-01T00:00:00", scale="utc")

    # Use a low altitude where refraction matters
    alt = 10*u.deg
    az = 90*u.deg

    apparent = AltAzTarget(alt=alt, az=az, location=location, name="apparent",
                           pressure=1*u.bar, temperature=10*u.deg_C, relative_humidity=0.5)
    vacuum = AltAzTarget(alt=alt, az=az, location=location, name="vacuum")

    ca = apparent.get_skycoord(t)
    cv = vacuum.get_skycoord(t)

    # Require a non-trivial difference
    assert ca.separation(cv) > 100*u.arcsec


@pytest.mark.skipif(not HAS_SGP4, reason="sgp4 is not installed")
def test_SGP4SatelliteTarget():
    tle_string = "\n".join(ISS_TLE)
    three_line_tle = "ISS (ZARYA)\n" + tle_string
    target = SGP4SatelliteTarget(tle=ISS_TLE, name="ISS")
    assert target.name == "ISS"
    assert target.is_time_dependent
    assert target.source_format == "tle"
    assert target.catalog_number == 25544
    assert str(target) == "ISS"
    assert repr(target) == '<SGP4SatelliteTarget "ISS" catalog #25544>'
    assert SGP4SatelliteTarget(tle=tle_string).catalog_number == 25544
    assert SGP4SatelliteTarget(tle=three_line_tle).name == "ISS (ZARYA)"
    assert SGP4SatelliteTarget(tle=three_line_tle, name="ISS").name == "ISS"
    SGP4SatelliteTarget(tle=ISS_TLE, validate=False)
    for model in ("wgs72", "wgs72old", "wgs84"):
        SGP4SatelliteTarget(tle=ISS_TLE, gravity_model=model)

    with pytest.raises(ValueError, match="exactly one"):
        SGP4SatelliteTarget()
    with pytest.raises(ValueError, match="exactly one"):
        SGP4SatelliteTarget(tle=ISS_TLE, omm=ISS_OMM)
    with pytest.raises(ValueError, match="2 or 3 lines"):
        SGP4SatelliteTarget(tle=ISS_TLE[0])
    with pytest.raises(ValueError, match="begin with"):
        SGP4SatelliteTarget(tle=(ISS_TLE[1], ISS_TLE[0]))
    with pytest.raises(ValueError):
        SGP4SatelliteTarget(tle=(ISS_TLE[0].replace("   23215", " 23215"), ISS_TLE[1]))
    with pytest.raises(ValueError, match="gravity_model"):
        SGP4SatelliteTarget(tle=ISS_TLE, gravity_model="other")

    json_text = json.dumps(ISS_OMM_GP)
    omm_targets = [
        SGP4SatelliteTarget(omm=ISS_OMM_GP),
        SGP4SatelliteTarget(omm=json_text),
        SGP4SatelliteTarget(omm=json.dumps([ISS_OMM_GP])),
        SGP4SatelliteTarget(omm=_omm_csv()),
        SGP4SatelliteTarget(omm="\n  " + _omm_csv()),
        SGP4SatelliteTarget(omm=_omm_xml()),
        SGP4SatelliteTarget(omm=_omm_kvn()),
        SGP4SatelliteTarget(omm=json_text, omm_format="json"),
        SGP4SatelliteTarget(omm=("\ufeff" + json_text).encode()),
        SGP4SatelliteTarget(omm=io.StringIO(json_text)),
    ]
    times = Time(
        ["2023-08-03T06:32:29.290272", "2023-08-03T06:42:29.290272"], scale="utc"
    )
    reference = target.get_teme(times)
    for omm_target in omm_targets:
        assert omm_target.name == "ISS (ZARYA)"
        assert omm_target.source_format == "omm"
        assert omm_target.catalog_number == 25544
        assert abs(omm_target.epoch - target.epoch) < 1 * u.us
        np.testing.assert_allclose(
            omm_target.get_teme(times).cartesian.xyz.to_value(u.km),
            reference.cartesian.xyz.to_value(u.km),
            atol=1e-6,
            rtol=0,
        )
        np.testing.assert_allclose(
            omm_target.get_teme(times)
            .cartesian.differentials["s"]
            .d_xyz.to_value(u.km / u.s),
            reference.cartesian.differentials["s"].d_xyz.to_value(u.km / u.s),
            atol=1e-9,
            rtol=0,
        )

    # CelesTrak XML may use null for an unavailable international designator
    xml = _omm_xml().replace(
        f"<OBJECT_ID>{ISS_OMM['OBJECT_ID']}</OBJECT_ID>", "<OBJECT_ID />"
    )
    target_without_object_id = SGP4SatelliteTarget(omm=xml)
    assert target_without_object_id.catalog_number == 25544
    assert target_without_object_id.name == "ISS (ZARYA)"

    xml = _omm_xml().replace(
        f"<OBJECT_NAME>{ISS_OMM['OBJECT_NAME']}</OBJECT_NAME>", "<OBJECT_NAME />"
    )
    assert SGP4SatelliteTarget(omm=xml).name == "Satellite 25544"

    for value in ("[]", json.dumps([ISS_OMM_GP, ISS_OMM_GP])):
        with pytest.raises(ValueError, match="exactly one OMM record"):
            SGP4SatelliteTarget(omm=value)
    with pytest.raises(ValueError, match="exactly one OMM record"):
        SGP4SatelliteTarget(omm=_omm_csv((ISS_OMM_GP, ISS_OMM_GP)))
    with pytest.raises(ValueError, match="exactly one OMM record"):
        SGP4SatelliteTarget(omm=_omm_xml((ISS_OMM, ISS_OMM)))

    for key, value in (
        ("CENTER_NAME", "MOON"),
        ("CENTER_NAME", "MARS"),
        ("REF_FRAME", "EME2000"),
        ("TIME_SYSTEM", "TAI"),
        ("MEAN_ELEMENT_THEORY", "SGP4-XP"),
        ("EPHEMERIS_TYPE", "4"),
        ("INCLINATION", "180"),
        ("INCLINATION", "nan"),
        ("RA_OF_ASC_NODE", "360"),
    ):
        fields = dict(ISS_OMM)
        fields[key] = value
        with pytest.raises(ValueError):
            SGP4SatelliteTarget(omm=fields)

    fields = dict(ISS_OMM_GP)
    fields["MEAN_ELEMENT_THEORY"] = "SGP/SGP4"
    SGP4SatelliteTarget(omm=fields)

    fields = dict(ISS_OMM)
    del fields["CENTER_NAME"]
    with pytest.raises(ValueError, match="metadata"):
        SGP4SatelliteTarget(omm=_omm_kvn(fields))
    with pytest.raises(ValueError, match="metadata"):
        SGP4SatelliteTarget(omm=_omm_xml((fields,)))

    fields = dict(ISS_OMM_GP)
    del fields["MEAN_MOTION"]
    fields["SEMI_MAJOR_AXIS"] = "6790"
    with pytest.raises(ValueError, match="SEMI_MAJOR_AXIS"):
        SGP4SatelliteTarget(omm=fields)

    fields = dict(ISS_OMM_GP)
    fields["MEAN_MOTION"] = "nan"
    with pytest.raises(ValueError, match="finite"):
        SGP4SatelliteTarget(omm=fields)

    fields = dict(ISS_OMM_GP)
    fields["epoch"] = fields["EPOCH"]
    with pytest.raises(ValueError, match="Duplicate"):
        SGP4SatelliteTarget(omm=fields)

    with pytest.raises(ValueError, match="Duplicate"):
        SGP4SatelliteTarget(omm=_omm_kvn() + "\nEPOCH = 2023-08-03T06:32:29.290272")
    with pytest.raises(ValueError, match="keyword"):
        SGP4SatelliteTarget(omm=_omm_kvn().replace("EPOCH =", "epoch ="))
    with pytest.raises(ValueError, match="KVN line"):
        SGP4SatelliteTarget(omm=_omm_kvn() + "\nINVALID")
    with pytest.raises(ValueError, match="Unsupported unit"):
        SGP4SatelliteTarget(
            omm=_omm_kvn().replace(
                "INCLINATION = 51.6403 [deg]", "INCLINATION = 51.6403 [rad]"
            )
        )

    scalar = target.get_teme(times[0])
    vector = target.get_teme(times)
    matrix_times = times[[0, 1, 1, 0]].reshape(2, 2)
    assert isinstance(scalar, TEME) and scalar.isscalar
    assert vector.shape == (2,)
    assert target.get_teme(matrix_times).shape == (2, 2)
    assert scalar.cartesian.xyz.unit == u.km
    assert scalar.cartesian.differentials["s"].d_xyz.unit == u.km / u.s
    assert abs(scalar.obstime - times[0]) < 1 * u.ns
    np.testing.assert_allclose(
        target.get_teme(times[0].tai).cartesian.xyz.to_value(u.km),
        scalar.cartesian.xyz.to_value(u.km),
        atol=1e-9,
        rtol=0,
    )

    observer = Observer(
        longitude=-155.476111 * u.deg, latitude=19.825555 * u.deg, elevation=4139 * u.m
    )
    other_observer = Observer(
        longitude=0 * u.deg, latitude=0 * u.deg, elevation=0 * u.m
    )
    with pytest.raises(ValueError, match="observer"):
        target.get_skycoord(times[0])
    assert (
        target.get_skycoord(times[0], observer).separation(
            target.get_skycoord(times[0], other_observer)
        )
        > 1 * u.arcsec
    )

    fixed = FixedTarget(SkyCoord(279.23458, 38.78369, unit="deg"), name="Vega")
    combined = get_skycoord([fixed, target], times=times, observer=observer)
    assert combined.is_equivalent_frame(ICRS())
    assert combined.shape == (2, 2)

    assert observer.altaz(
        times[0], target, grid_times_targets=True
    ).shape == (1, 1)
    assert observer.altaz(
        times[0], [fixed, target], grid_times_targets=True
    ).shape == (2, 1)
    assert observer.altaz(
        times, [fixed, target], grid_times_targets=True
    ).shape == (2, 2)

    with pytest.raises(ValueError, match="`times` is required"):
        target.get_teme(None)
    with pytest.raises(ValueError, match="`times` is required"):
        target.get_skycoord(None, observer=observer)


@pytest.mark.skipif(not HAS_SGP4, reason="sgp4 is not installed")
def test_SGP4SatelliteTarget_accuracy():
    target = SGP4SatelliteTarget(tle=ISS_TLE)
    observer = Observer(
        longitude=-155.476111 * u.deg, latitude=19.825555 * u.deg, elevation=4139 * u.m
    )

    # Below Horizon
    time = Time("2023-08-02 10:00", scale="utc")
    # '08h29m27.30648375s +07d31m31.61825139s'
    ra_dec = target.get_skycoord(time, observer=observer)

    # Comparison with the JPL Horizons System
    ra_dec_horizon_icrf = SkyCoord("08h29m27.029117s +07d31m28.35610s")
    # ICRF: Compensated for the down-leg light-time delay aberration
    assert ra_dec.separation(ra_dec_horizon_icrf) < 10 * u.arcsec  # 5.25″
    # Distance estimation: ~ 2 * tan(5.25/2/3600) * 11801.56 = 17 km

    ra_dec_horizon_ref_apparent = SkyCoord("08h30m54.567398s +08d05m32.72764s")
    # Refracted Apparent: In an equatorial coordinate system with all compensations
    assert ra_dec.separation(ra_dec_horizon_ref_apparent) > 2000 * u.arcsec  # 2418.21″

    ra_dec_horizon_icrf_ref_apparent = SkyCoord("08h29m37.373866s +08d10m14.78811s")
    # ICRF Refracted Apparent: In the ICRF reference frame with all compensations
    assert (
        ra_dec.separation(ra_dec_horizon_icrf_ref_apparent) > 2000 * u.arcsec
    )  # 2327.98″

    # Above Horizon
    time = Time("2023-08-02 07:20", scale="utc")
    # '11h19m49.85714655s +44d49m34.4600926s'
    ra_dec_ah = target.get_skycoord(time, observer=observer)
    ra_dec_ah_horizon_icrf = SkyCoord("11h19m49.660349s +44d49m34.65875s")
    assert ra_dec_ah.separation(ra_dec_ah_horizon_icrf) < 10 * u.arcsec  # 2.10″

    teme = target.get_teme(time)
    geocentric_itrs = teme.transform_to(ITRS(obstime=time))
    topocentric_position = (
        geocentric_itrs.cartesian.without_differentials()
        - observer.location.get_itrs(time).cartesian
    )
    manual = ITRS(
        topocentric_position, obstime=time, location=observer.location
    ).transform_to(AltAz(obstime=time, location=observer.location))
    assert manual.separation(observer.altaz(time, target)) < 1e-3 * u.arcsec
    # Direct transformation incorrectly changes stellar aberration
    direct = teme.transform_to(AltAz(obstime=time, location=observer.location))
    assert manual.alt > 0 * u.deg
    assert manual.separation(direct) > 20 * u.arcsec  # 62.70″


@pytest.mark.skipif(not HAS_SGP4, reason="sgp4 is not installed")
def test_SGP4SatelliteTarget_export():
    target = SGP4SatelliteTarget(tle=ISS_TLE)

    line1, line2 = target.to_tle()
    assert line1.startswith("1 25544")
    assert line2.startswith("2 25544")
    assert len(line1) == 69
    assert len(line2) == 69

    roundtrip_tle = SGP4SatelliteTarget(tle=(line1, line2))
    assert roundtrip_tle.catalog_number == target.catalog_number
    assert roundtrip_tle.epoch == target.epoch

    omm = target.to_omm()
    assert omm["OBJECT_NAME"] == target.name
    assert omm["NORAD_CAT_ID"] == target.catalog_number
    assert omm["CENTER_NAME"] == "EARTH"
    assert omm["REF_FRAME"] == "TEME"
    assert omm["TIME_SYSTEM"] == "UTC"
    assert omm["MEAN_ELEMENT_THEORY"] == "SGP4"

    roundtrip_omm = SGP4SatelliteTarget(omm=omm)
    assert roundtrip_omm.to_tle() == target.to_tle()

    # TLE Alpha-5 support
    omm = """[{
        "OBJECT_NAME": "STARLINK-38222",
        "OBJECT_ID": "2026-166A",
        "EPOCH": "2026-08-18T22:42:27.673632",
        "MEAN_MOTION": 15.77313337,
        "ECCENTRICITY": 0.00058888,
        "INCLINATION": 97.2853,
        "RA_OF_ASC_NODE": 247.7747,
        "ARG_OF_PERICENTER": 131.8444,
        "MEAN_ANOMALY": 228.3326,
        "EPHEMERIS_TYPE": 0,
        "CLASSIFICATION_TYPE": "U",
        "NORAD_CAT_ID": 100101,
        "ELEMENT_SET_NO": 999,
        "REV_AT_EPOCH": 448,
        "BSTAR": -0.00083032139,
        "MEAN_MOTION_DOT": -0.00134116,
        "MEAN_MOTION_DDOT": 0
    }]"""
    target = SGP4SatelliteTarget(omm=omm)
    assert target.catalog_number == 100101
    omm = target.to_omm()
    line1, line2 = target.to_tle()
    assert omm["NORAD_CAT_ID"] == 100101
    assert line1.startswith("1 A0101U")  # Alpha-5

    roundtrip_tle = SGP4SatelliteTarget(tle=(line1, line2))
    assert roundtrip_tle.catalog_number == target.catalog_number
    omm = roundtrip_tle.to_omm()
    line1, line2 = roundtrip_tle.to_tle()
    assert omm["NORAD_CAT_ID"] == 100101
    assert line1.startswith("1 A0101U")  # Alpha-5


@pytest.mark.skipif(not HAS_SGP4, reason="sgp4 is not installed")
def test_SGP4SatelliteTarget_propagation_warning():
    class FakeSatrec:
        def sgp4_array(self, jd, fraction):
            """Returns last entry as erroneous"""
            count = len(jd)
            errors = np.zeros(count, dtype=int)
            errors[-1] = 6
            return errors, np.ones((count, 3)), np.ones((count, 3))

    target = SGP4SatelliteTarget(tle=ISS_TLE)
    observer = Observer(
        longitude=-155.476111 * u.deg, latitude=19.825555 * u.deg, elevation=4139 * u.m
    )
    target._satrec = FakeSatrec()
    times = Time(["2023-08-03T06:32:29", "2023-08-03T06:33:29"], scale="utc")
    with pytest.warns(SatellitePropagationWarning, match="1 of 2 times"):
        teme = target.get_teme(times)
    assert np.all(np.isfinite(teme.cartesian.xyz[:, 0]))
    assert np.all(np.isnan(teme.cartesian.xyz[:, 1]))

    target = SGP4SatelliteTarget(tle=ISS_TLE)
    time_invalid = Time("2025-08-02 10:00", scale="utc")
    with pytest.warns(SatellitePropagationWarning, match="SGP4 propagation failed"):
        assert np.isnan(target.get_skycoord(time_invalid, observer=observer).ra)


@pytest.mark.remote_data
def test_get_skycoord():
    m31 = SkyCoord(10.6847083*u.deg, 41.26875*u.deg)
    m31_with_distance = SkyCoord(10.6847083*u.deg, 41.26875*u.deg, 780*u.kpc)
    subaru = Observer.at_site('subaru')
    time = Time("2016-01-22 12:00")
    pos, vel = subaru.location.get_gcrs_posvel(time)
    gcrs_frame = GCRS(obstime=Time("2016-01-22 12:00"), obsgeoloc=pos, obsgeovel=vel)
    m31_gcrs = m31.transform_to(gcrs_frame)
    m31_gcrs_with_distance = m31_with_distance.transform_to(gcrs_frame)

    coo = get_skycoord(m31)
    assert coo.is_equivalent_frame(ICRS())
    with pytest.raises(TypeError):
        len(coo)

    coo = get_skycoord([m31])
    assert coo.is_equivalent_frame(ICRS())
    assert len(coo) == 1

    coo = get_skycoord([m31, m31_gcrs])
    assert coo.is_equivalent_frame(ICRS())
    assert len(coo) == 2

    coo = get_skycoord([m31_with_distance, m31_gcrs_with_distance])
    assert coo.is_equivalent_frame(ICRS())
    assert len(coo) == 2

    coo = get_skycoord([m31, m31_gcrs, m31_gcrs_with_distance, m31_with_distance])
    assert coo.is_equivalent_frame(ICRS())
    assert len(coo) == 4

    coo = get_skycoord([m31_gcrs, m31_gcrs_with_distance])
    assert coo.is_equivalent_frame(m31_gcrs.frame)
    assert len(coo) == 2


def test_get_skycoord_broadcasts_fixed_targets_when_time_dependent_present():
    """
    When at least one target is time-dependent and `times` is provided,
    `get_skycoord` should broadcast fixed targets to match `times.shape` and
    stack along the target axis.
    """
    location = EarthLocation.from_geodetic(lon=0 * u.deg, lat=0 * u.deg, height=0 * u.m)
    altaz_target = AltAzTarget(alt=45 * u.deg, az=0 * u.deg, location=location)

    m31 = SkyCoord(10.6847083 * u.deg, 41.26875 * u.deg)

    t0 = Time("2026-02-05 00:00")
    times = t0 + np.arange(4) * u.hour

    coo = get_skycoord([m31, altaz_target], times=times)

    assert coo.is_equivalent_frame(ICRS())
    assert coo.shape == (2,) + times.shape

    # Fixed target should be repeated across time
    assert np.allclose(coo[0].ra.to_value(u.deg), m31.ra.to_value(u.deg))
    assert np.allclose(coo[0].dec.to_value(u.deg), m31.dec.to_value(u.deg))

    # Time-dependent target should vary across time
    assert coo[1][0].separation(coo[1][-1]) > 10 * u.deg


@pytest.mark.skipif(not HAS_SGP4, reason="sgp4 is not installed")
def test_get_skycoord_with_SGP4SatelliteTarget():
    skycoord_target = SkyCoord(10.6847083*u.deg, 41.26875*u.deg)
    fixed_target = FixedTarget(name="fixed1", coord=SkyCoord(279.23458, 38.78369, unit='deg'))

    observer = Observer(longitude=-155.476111*u.deg, latitude=19.825555*u.deg, elevation=4139*u.m)
    tle_target = SGP4SatelliteTarget(tle=ISS_TLE)
    time = Time("2023-08-02 10:00", scale="utc")
    times = time_grid_from_range(
        [time, time + 3.1 * u.hour], time_resolution=1 * u.hour
    )

    with pytest.raises(ValueError, match="observer"):
        get_skycoord(tle_target, time)

    tle_output = get_skycoord(tle_target, time, observer=observer)
    assert tle_output.size == 1

    tle_output = get_skycoord(tle_target, times, observer=observer)
    assert tle_output.shape == (4,)

    tle_output = get_skycoord([tle_target, tle_target], time, observer=observer)
    assert tle_output.shape == (2,)

    tle_output = get_skycoord([tle_target, tle_target], times, observer=observer)
    assert tle_output.shape == (2, 4)

    mixed_output = get_skycoord(
        [skycoord_target, fixed_target, tle_target], time, observer=observer
    )
    assert mixed_output.shape == (3,)

    mixed_output = get_skycoord(
        [skycoord_target, fixed_target, tle_target], times, observer=observer
    )
    assert mixed_output.shape == (3, 4)


def test_get_skycoord_does_not_broadcast_when_all_targets_are_fixed():
    """
    If all targets are fixed, providing `times` should not change the output
    shape (backwards-compatible behavior).
    """
    m31 = SkyCoord(10.6847083 * u.deg, 41.26875 * u.deg)
    m32 = SkyCoord(10.6747083 * u.deg, 40.26875 * u.deg)

    t0 = Time("2026-02-05 00:00")
    times = t0 + np.arange(3) * u.hour

    coo = get_skycoord([m31, m32], times=times)

    assert coo.is_equivalent_frame(ICRS())
    assert coo.shape == (2,)


def test_get_skycoord_mixed_distances_with_time_dependent_target_fills_unitspherical():
    """
    With a mixture of targets with distances and unit-spherical targets, and at
    least one time-dependent target present, `get_skycoord` should return a
    distance-bearing SkyCoord and fill large distances for unit-spherical entries.
    """
    location = EarthLocation.from_geodetic(lon=0 * u.deg, lat=0 * u.deg, height=0 * u.m)
    altaz_target = AltAzTarget(alt=45 * u.deg, az=0 * u.deg, location=location)

    m31 = SkyCoord(10.6847083 * u.deg, 41.26875 * u.deg)  # unit-spherical
    m31_with_distance = SkyCoord(10.6847083 * u.deg, 41.26875 * u.deg, 780 * u.kpc)

    t0 = Time("2026-02-05 00:00")
    times = t0 + np.arange(4) * u.hour

    coo = get_skycoord([m31, m31_with_distance, altaz_target], times=times)

    assert coo.is_equivalent_frame(ICRS())
    assert coo.shape == (3,) + times.shape
    assert coo.distance.shape == (3,) + times.shape

    # Filled distances for unit-spherical targets
    assert np.allclose(coo.distance[0].to_value(u.kpc), 100.0)
    assert np.allclose(coo.distance[2].to_value(u.kpc), 100.0)
    # Preserved distance for the distance-bearing target
    assert np.allclose(coo.distance[1].to_value(u.kpc), 780.0)


def test_get_skycoord_forwards_observer(observer_dependent_target):
    location = EarthLocation.from_geodetic(
        lon=10 * u.deg,
        lat=45 * u.deg,
        height=100 * u.m,
    )
    observer = Observer(location=location)
    target = observer_dependent_target
    times = Time(["2026-02-05T00:00:00", "2026-02-05T01:00:00"], scale="utc")

    with pytest.raises(ValueError, match="observer"):
        get_skycoord(target, times=times)

    coord = get_skycoord(target, times=times, observer=observer)

    assert coord.shape == times.shape
    assert np.allclose(coord.ra.to_value(u.deg), observer.location.lon.to_value(u.deg))
    assert np.allclose(coord.dec.to_value(u.deg), observer.location.lat.to_value(u.deg))

    combined = get_skycoord(
        [SkyCoord(ra=0 * u.deg, dec=0 * u.deg), target], times=times, observer=observer
    )

    assert combined.shape == (2,) + times.shape

    altaz = observer.altaz(times, target)

    assert altaz.shape == times.shape
