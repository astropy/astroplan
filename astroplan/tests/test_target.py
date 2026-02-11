# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Third-party
import numpy as np
import astropy.units as u
import pytest
from astropy.coordinates import SkyCoord, GCRS, ICRS, EarthLocation
from astropy.time import Time

# Package
from astroplan.target import FixedTarget, AltAzTarget, get_skycoord
from astroplan.observer import Observer


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
