# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

# Third-party
import astropy.units as u
from astropy.coordinates import SkyCoord, GCRS, ICRS
from astropy.time import Time
import numpy as np
try:
    import skyfield  # noqa
    HAS_SKYFIELD = True
except ImportError:
    HAS_SKYFIELD = False

# Package
from ..target import FixedTarget, TLETarget, get_skycoord
from ..observer import Observer
from ..utils import time_grid_from_range
from ..exceptions import InvalidTLEDataWarning


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


@pytest.mark.skipif('not HAS_SKYFIELD')
def test_TLETarget():
    tle_string = ("ISS (ZARYA)\n"
                  "1 25544U 98067A   23215.27256123  .00041610  00000-0  73103-3 0  9990\n"
                  "2 25544  51.6403  95.2411 0000623 157.9606 345.0624 15.50085581409092")
    line1="1 25544U 98067A   23215.27256123  .00041610  00000-0  73103-3 0  9990"
    line2="2 25544  51.6403  95.2411 0000623 157.9606 345.0624 15.50085581409092"
    subaru = Observer.at_site('subaru')
    time = Time("2023-08-02 10:00", scale='utc')
    times = time_grid_from_range([time, time + 3.1*u.hour],
                                 time_resolution=1 * u.hour)

    tle_target1 = TLETarget(name="ISS (ZARYA)", line1=line1, line2=line2, observer=subaru)
    tle_target2 = TLETarget.from_string(tle_string=tle_string, observer=subaru)

    assert tle_target1.name == "ISS (ZARYA)"
    assert tle_target2.name == "ISS (ZARYA)"
    assert repr(tle_target1) == repr(tle_target2)
    assert str(tle_target1) == str(tle_target2)

    # Single time
    ra_dec1 = tle_target1.coord(time)
    ra_dec2 = tle_target2.coord(time)
    ra_dec_from_horizon = SkyCoord("08h29m26.36s +07d31m24.0s")

    assert ra_dec1.separation(ra_dec_from_horizon) < 30*u.arcsec
    assert ra_dec1.to_string('hmsdms') == ra_dec2.to_string('hmsdms')

    # Multiple times
    ra_dec1 = tle_target1.coord(times)
    ra_dec2 = tle_target2.coord(times)
    ra_dec_from_horizon = SkyCoord(["08h29m26.36s +07d31m24.0s",
                                    "06h25m46.81s -54d32m20.2s",
                                    "13h52m09.08s +04d26m40.1s",
                                    "09h20m04.74s -00d51m19.2s"])

    assert all(list(ra_dec1.separation(ra_dec_from_horizon) < 30*u.arcsec))
    assert ra_dec1.to_string('hmsdms') == ra_dec2.to_string('hmsdms')

    # TLE Check
    line1_invalid = "1 25544U 98067A .00041610  00000-0  73103-3 0  9990"
    with pytest.raises(ValueError):
        TLETarget(name="ISS (ZARYA)", line1=line1_invalid, line2=line2, observer=subaru)
    TLETarget(name="ISS (ZARYA)", line1=line1_invalid, line2=line2, observer=subaru,
              skip_tle_check=True)

    # from_string
    tle_string = ("1 25544U 98067A   23215.27256123  .00041610  00000-0  73103-3 0  9990\n"
                  "2 25544  51.6403  95.2411 0000623 157.9606 345.0624 15.50085581409092")
    tle_target3 = TLETarget.from_string(tle_string=tle_string, observer=subaru, name="ISS")
    assert tle_target3.name == "ISS"

    tle_string = "1 25544U 98067A   23215.27256123  .00041610  00000-0  73103-3 0  9990"
    with pytest.raises(ValueError):
        TLETarget.from_string(tle_string=tle_string, observer=subaru)

    # AltAz (This is slow)
    altaz_observer = subaru.altaz(time, tle_target1)
    altaz_skyfield = tle_target1.altaz(time)

    assert altaz_observer.separation(altaz_skyfield) < 20*u.arcsec

    # Time too far in the future where elements stop making physical sense
    with pytest.warns():    # ErfaWarning: ERFA function "dtf2d" yielded 1 of "dubious year (Note 6)
        time_invalid = Time("2035-08-02 10:00", scale='utc')
        times_list = list(times)
        times_list[2] = Time("2035-08-02 10:00", scale='utc')
        times_invalid = Time(times_list)
    with pytest.warns(InvalidTLEDataWarning):
        assert np.isnan(tle_target1.coord(time_invalid).ra)
    with pytest.warns(InvalidTLEDataWarning):
        assert np.isnan(tle_target1.coord(times_invalid)[2].ra)


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


@pytest.mark.skipif('not HAS_SKYFIELD')
def test_get_skycoord_with_TLETarget():
    skycoord_targed = SkyCoord(10.6847083*u.deg, 41.26875*u.deg)
    fixed_target1 = FixedTarget(name="fixed1", coord=SkyCoord(279.23458, 38.78369, unit='deg'))
    line1 = "1 25544U 98067A   23215.27256123  .00041610  00000-0  73103-3 0  9990"
    line2 = "2 25544  51.6403  95.2411 0000623 157.9606 345.0624 15.50085581409092"
    subaru = Observer.at_site('subaru')
    tle_target1 = TLETarget(name="ISS (ZARYA)", line1=line1, line2=line2, observer=subaru)

    time = Time("2023-08-02 10:00", scale='utc')
    times = time_grid_from_range([time, time + 3.1*u.hour],
                                 time_resolution=1 * u.hour)

    tle_output = get_skycoord(tle_target1, time)
    assert tle_output.size == 1
    tle_output = get_skycoord(tle_target1, times)
    assert tle_output.shape == (4,)
    tle_output = get_skycoord([tle_target1, tle_target1], time)
    assert tle_output.shape == (2,)
    tle_output = get_skycoord([tle_target1, tle_target1], times)
    assert tle_output.shape == (2, 4)

    mixed_output = get_skycoord([skycoord_targed, fixed_target1, tle_target1], time)
    assert mixed_output.shape == (3,)
    mixed_output = get_skycoord([skycoord_targed, fixed_target1, tle_target1], times)
    assert mixed_output.shape == (3, 4)

    # backwards_compatible
    tle_output = get_skycoord(skycoord_targed, time, backwards_compatible=False)
    assert tle_output.size == 1

    tle_output = get_skycoord(skycoord_targed, times)
    assert tle_output.size == 1
    tle_output = get_skycoord(skycoord_targed, times, backwards_compatible=False)
    assert tle_output.shape == (4,)

    tle_output = get_skycoord([skycoord_targed, skycoord_targed], time, backwards_compatible=False)
    assert tle_output.shape == (2,)

    tle_output = get_skycoord([skycoord_targed, skycoord_targed], times)
    assert tle_output.shape == (2,)
    tle_output = get_skycoord([skycoord_targed, skycoord_targed], times, backwards_compatible=False)
    assert tle_output.shape == (2, 4)
