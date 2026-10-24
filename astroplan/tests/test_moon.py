# Licensed under a 3-clause BSD style license - see LICENSE.rst
import astropy.units as u
import numpy as np
from astropy.coordinates import EarthLocation
from astropy.time import Time
from numpy.testing import assert_allclose

from astroplan.moon import moon_illumination, moon_phase_angle
from astroplan.observer import Observer


def test_illumination():
    time = Time(['1990-01-01 00:00:00', '1990-03-01 06:00:00',
                 '1990-06-01 12:00:00', '1990-11-01 18:00:00'])
    location = EarthLocation.from_geodetic(-155*u.deg, 19*u.deg, 0*u.m)
    obs = Observer(location)
    # Get illumination via time
    illumination1 = obs.moon_illumination(time)

    # Run print_pyephem_illumination() for PyEphem's solution
    pyephem_illumination = [0.15475513880925418, 0.19484233284757257,
                            0.6170840254669668, 0.9780219372563843]

    assert_allclose(illumination1, pyephem_illumination, atol=0.05)


def print_pyephem_illumination():
    """
    To run, use:
    python -c "from astroplan.tests.test_moon import print_pyephem_illumination as f; f()"
    """
    time = Time(['1990-01-01 00:00:00', '1990-03-01 06:00:00',
                 '1990-06-01 12:00:00', '1990-11-01 18:00:00'])
    location = EarthLocation.from_geodetic(-155*u.deg, 19*u.deg, 0*u.m)

    import ephem
    moon = ephem.Moon()
    pe_obs = ephem.Observer()
    pe_obs.lat = location.lat.to_string(u.deg, sep=':')
    pe_obs.lon = location.lon.to_string(u.deg, sep=':')
    pe_obs.elevation = location.height.to(u.m).value
    illuminations = []
    for t in time:
        pe_obs.date = t.datetime
        moon.compute(pe_obs)
        illuminations.append(moon.moon_phase)
    print(illuminations)


def test_illumination_from_observer_cache():
    obs = Observer.at_site("APO")
    times = Time("2020-01-01 00:00") + np.linspace(0, 60, 50) * u.day

    # the observer path uses the cached topocentric Moon shifted back to the
    # geocentre; it must reproduce the direct geocentric computation
    illum = moon_illumination(times, observer=obs)
    assert_allclose(illum, moon_illumination(times), atol=1e-5)
    dphase = abs(moon_phase_angle(times, observer=obs) - moon_phase_angle(times))
    assert (dphase < 2 * u.arcsec).all()

    # only the Moon is cached (the Sun comes from get_sun), and a second call
    # through the Observer reuses that entry rather than adding to the cache
    assert len(obs._body_cache) == 1
    assert_allclose(obs.moon_illumination(times), illum)
    assert len(obs._body_cache) == 1

    # once the Sun is cached too, it is used instead of get_sun
    obs.get_body("sun", times)
    assert_allclose(obs.moon_illumination(times), illum, atol=1e-5)
    assert len(obs._body_cache) == 2


def test_moon_phase_default_time():
    # with no time given, moon_phase uses the current time (#663)
    obs = Observer.at_site("APO")
    now = Time.now()
    phase = obs.moon_phase()

    assert phase.isscalar
    # the phase changes by ~2e-6 rad per second
    assert_allclose(phase.to_value(u.rad), obs.moon_phase(now).to_value(u.rad),
                    atol=1e-4)


def test_moon_illumination_default_time():
    # with no time given, moon_illumination uses the current time, like moon_phase
    obs = Observer.at_site("APO")
    now = Time.now()
    illumination = obs.moon_illumination()

    assert np.ndim(illumination) == 0
    assert_allclose(illumination, obs.moon_illumination(now), atol=1e-4)
