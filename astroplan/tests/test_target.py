# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Third-party
import astropy.units as u
import pytest
from astropy.coordinates import SkyCoord, GCRS, ICRS, EarthLocation
from astropy.time import Time
import numpy as np
try:
    import skyfield  # noqa
    HAS_SKYFIELD = True
except ImportError:
    HAS_SKYFIELD = False

# Package
from astroplan.target import Target, FixedTarget, TLETarget, get_skycoord
from astroplan.observer import Observer
from astroplan.utils import time_grid_from_range


class ObserverDependentTarget(Target):
    """Minimal target that requires observer context when evaluated."""

    def __init__(self, name="observer-dependent"):
        self.name = name

    @property
    def is_time_dependent(self):
        return True

    def get_skycoord(self, times, observer=None):
        if observer is None:
            raise ValueError("`observer` is required.")

        ra = np.broadcast_to(
            observer.location.lon.to_value(u.deg),
            times.shape
        ) * u.deg
        dec = np.broadcast_to(
            observer.location.lat.to_value(u.deg),
            times.shape,
        ) * u.deg

        return SkyCoord(ra=ra, dec=dec)


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


@pytest.mark.remote_data
@pytest.mark.skipif(not HAS_SKYFIELD, reason="skyfield is not installed")
def test_TLETarget():
    tle_string = ("ISS (ZARYA)\n"
                  "1 25544U 98067A   23215.27256123  .00041610  00000-0  73103-3 0  9990\n"
                  "2 25544  51.6403  95.2411 0000623 157.9606 345.0624 15.50085581409092")
    line1 = "1 25544U 98067A   23215.27256123  .00041610  00000-0  73103-3 0  9990"
    line2 = "2 25544  51.6403  95.2411 0000623 157.9606 345.0624 15.50085581409092"
    subaru = Observer.at_site('subaru')  # (lon, lat, el)=(-155.476111 deg, 19.825555 deg, 4139.0 m)
    subaru_temp_pressure = Observer.at_site('subaru', temperature=-10*u.Celsius, pressure=0.9*u.bar)
    time = Time("2023-08-02 10:00", scale='utc')
    times = time_grid_from_range([time, time + 3.1*u.hour],
                                 time_resolution=1 * u.hour)

    tle_target1 = TLETarget(name="ISS (ZARYA)", line1=line1, line2=line2)
    tle_target2 = TLETarget.from_string(tle_string=tle_string)

    assert tle_target1.name == "ISS (ZARYA)"
    assert tle_target2.name == "ISS (ZARYA)"
    assert repr(tle_target1) == repr(tle_target2)
    assert str(tle_target1) == str(tle_target2)

    with pytest.raises(ValueError, match="observer"):
        tle_target1.get_skycoord(time)

    # Single time (Below Horizon)
    ra_dec1 = tle_target1.get_skycoord(time, observer=subaru)   # '08h29m26.00003243s +07d31m36.65950907s'
    ra_dec2 = tle_target2.get_skycoord(time, observer=subaru)
    assert ra_dec1.to_string('hmsdms') == ra_dec2.to_string('hmsdms')

    # Comparison with the JPL Horizons System
    ra_dec_horizon_icrf = SkyCoord("08h29m27.029117s +07d31m28.35610s")
    # ICRF: Compensated for the down-leg light-time delay aberration
    assert ra_dec1.separation(ra_dec_horizon_icrf) < 20*u.arcsec  # 17.41″
    # Distance estimation: ~ 2 * tan(17,41/2/3600) * 11801,56 = 57 km

    ra_dec_horizon_ref_apparent = SkyCoord("08h30m54.567398s +08d05m32.72764s")
    # Refracted Apparent: In an equatorial coordinate system with all compensations
    assert ra_dec1.separation(ra_dec_horizon_ref_apparent) > 2000*u.arcsec  # 2424.44″

    ra_dec_horizon_icrf_ref_apparent = SkyCoord("08h29m37.373866s +08d10m14.78811s")
    # ICRF Refracted Apparent: In the ICRF reference frame with all compensations
    assert ra_dec1.separation(ra_dec_horizon_icrf_ref_apparent) > 2000*u.arcsec  # 2324.28″

    # Skyfield appears to use no compensations. According to this, it's not even recommended to
    # compensate for light travel time. Compensating changes the difference to ra_dec_horizon_icrf
    # to 20.05 arcsec.
    # https://rhodesmill.org/skyfield/earth-satellites.html#avoid-calling-the-observe-method

    # Single time (Above Horizon)
    time_ah = Time("2023-08-02 07:20", scale='utc')
    ra_dec_ah = tle_target1.get_skycoord(time_ah, observer=subaru)  # '11h19m48.53631001s +44d49m45.22194611s'

    ra_dec_ah_horizon_icrf = SkyCoord("11h19m49.660349s +44d49m34.65875s")
    assert ra_dec_ah.separation(ra_dec_ah_horizon_icrf) < 20*u.arcsec  # 15.95″
    ra_dec_ah_horizon_ref_apparent = SkyCoord("11h21m34.102381s +44d43m40.06899s")
    assert ra_dec_ah.separation(ra_dec_ah_horizon_ref_apparent) > 1000*u.arcsec  # 1181.84″
    ra_dec_ah_horizon_icrf_ref_apparent = SkyCoord("11h20m16.627261s +44d51m24.25337s")
    assert ra_dec_ah.separation(ra_dec_ah_horizon_icrf_ref_apparent) > 300*u.arcsec  # 314.75″

    # Default is WGS72 for Skyfield. Coordinates with WGS84 gravity model that Horizon uses:
    # '11h19m48.28084569s +44d49m46.33649241s' - 18.75″
    # See 'Build a satellite with a specific gravity model' in Skyfield's Earth Satellites docu

    # There are many potential sources of inaccuracies, and it's not all super precise.
    # Should the accuracy be better than < 25*u.arcsec when compared to the JPL Horizons System?

    # A single orbit can be evaluated independently from different sites.
    greenwich = Observer(
        latitude=51.4779 * u.deg, longitude=0 * u.deg, elevation=46 * u.m
    )
    ra_dec_greenwich = tle_target1.get_skycoord(time_ah, observer=greenwich)    
    assert ra_dec_ah.separation(ra_dec_greenwich) > 1 * u.arcmin

    # Multiple times
    ra_dec1 = tle_target1.get_skycoord(times, observer=subaru)
    ra_dec2 = tle_target2.get_skycoord(times, observer=subaru)
    ra_dec_from_horizon = SkyCoord(["08h29m27.029117s +07d31m28.35610s",
                                    "06h25m46.672661s -54d32m16.77533s",
                                    "13h52m08.854291s +04d26m49.56432s",
                                    "09h20m04.872215s -00d51m21.17432s"])
    # 17.41″, 22.05″, 3.20″, 17.55″
    assert all(list(ra_dec1.separation(ra_dec_from_horizon) < 25*u.arcsec))
    assert ra_dec1.to_string('hmsdms') == ra_dec2.to_string('hmsdms')

    # TLE Check
    line1_invalid = "1 25544U 98067A .00041610  00000-0  73103-3 0  9990"
    with pytest.raises(ValueError):
        TLETarget(name="ISS (ZARYA)", line1=line1_invalid, line2=line2)
    TLETarget(name="ISS (ZARYA)", line1=line1_invalid, line2=line2, skip_tle_check=True)

    # from_string
    tle_string = ("1 25544U 98067A   23215.27256123  .00041610  00000-0  73103-3 0  9990\n"
                  "2 25544  51.6403  95.2411 0000623 157.9606 345.0624 15.50085581409092")
    tle_target3 = TLETarget.from_string(tle_string=tle_string, name="ISS")
    assert tle_target3.name == "ISS"

    tle_string = "1 25544U 98067A   23215.27256123  .00041610  00000-0  73103-3 0  9990"
    with pytest.raises(ValueError):
        TLETarget.from_string(tle_string=tle_string, observer=subaru)

    # AltAz (This is slow)
    altaz_observer = subaru.altaz(time_ah, tle_target1)
    altaz_skyfield = tle_target1.altaz(time_ah, observer=subaru)

    assert altaz_observer.separation(altaz_skyfield) < 21*u.arcsec

    # AltAz with atmospheric refraction (temperature and pressure)
    altaz_observer = subaru_temp_pressure.altaz(time_ah, tle_target1)
    altaz_skyfield = tle_target1.altaz(time_ah, observer=subaru_temp_pressure)

    assert altaz_observer.separation(altaz_skyfield) < 26*u.arcsec

    # AltAz with multiple times
    # subaru.altaz(times, tle_target1)
    altaz_multiple = tle_target1.altaz(times, observer=subaru)
    assert len(altaz_multiple.obstime) == len(altaz_multiple) == len(times)

    # Time too far in the future where elements stop making physical sense
    with pytest.warns():  # ErfaWarning: ERFA function "dtf2d" yielded 1 of "dubious year (Note 6)
        time_invalid = Time("2035-08-02 10:00", scale='utc')
        times_list = list(times)
        times_list[2] = Time("2035-08-02 10:00", scale='utc')
        times_invalid = Time(times_list)
    # InvalidTLEDataWarning and
    # ErfaWarning: ERFA function "utctai" yielded 1 of "dubious year (Note 3)"
    with pytest.warns():
        assert np.isnan(tle_target1.get_skycoord(time_invalid, observer=subaru).ra)
    # InvalidTLEDataWarning and
    # ErfaWarning: ERFA function "utctai" yielded 1 of "dubious year (Note 3)"
    with pytest.warns():
        assert np.isnan(tle_target1.get_skycoord(times_invalid, observer=subaru)[2].ra)


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


@pytest.mark.remote_data
@pytest.mark.skipif(not HAS_SKYFIELD, reason="skyfield is not installed")
def test_get_skycoord_with_TLETarget():
    skycoord_target = SkyCoord(10.6847083*u.deg, 41.26875*u.deg)
    fixed_target = FixedTarget(name="fixed1", coord=SkyCoord(279.23458, 38.78369, unit='deg'))
    line1 = "1 25544U 98067A   23215.27256123  .00041610  00000-0  73103-3 0  9990"
    line2 = "2 25544  51.6403  95.2411 0000623 157.9606 345.0624 15.50085581409092"
    subaru = Observer.at_site('subaru')
    tle_target = TLETarget(name="ISS (ZARYA)", line1=line1, line2=line2)

    time = Time("2023-08-02 10:00", scale="utc")
    times = time_grid_from_range(
        [time, time + 3.1 * u.hour], time_resolution=1 * u.hour
    )

    with pytest.raises(ValueError, match="observer"):
        get_skycoord(tle_target, time)

    tle_output = get_skycoord(tle_target, time, observer=subaru)
    assert tle_output.size == 1

    tle_output = get_skycoord(tle_target, times, observer=subaru)
    assert tle_output.shape == (4,)

    tle_output = get_skycoord([tle_target, tle_target], time, observer=subaru)
    assert tle_output.shape == (2,)

    tle_output = get_skycoord([tle_target, tle_target], times, observer=subaru)
    assert tle_output.shape == (2, 4)

    mixed_output = get_skycoord(
        [skycoord_target, fixed_target, tle_target], time, observer=subaru
    )
    assert mixed_output.shape == (3,)

    mixed_output = get_skycoord(
        [skycoord_target, fixed_target, tle_target], times, observer=subaru
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


def test_get_skycoord_forwards_observer():
    location = EarthLocation.from_geodetic(
        lon=10 * u.deg,
        lat=45 * u.deg,
        height=100 * u.m,
    )
    observer = Observer(location=location)
    target = ObserverDependentTarget()
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
    