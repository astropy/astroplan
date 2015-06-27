from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from astropy.coordinates import (EarthLocation, Latitude, Longitude, SkyCoord)
import astropy.units as u
from astropy.time import Time
from astropy.tests.helper import remote_data
import numpy as np
from numpy.testing import assert_allclose
import pytz
import datetime
import unittest

from ..core import FixedTarget, Observer


def test_Observer_constructor_location():
    '''
    Show that location defined by latitude/longitude/elevation is parsed
    identically to passing in an `~astropy.coordinates.EarthLocation` directly.
    '''

    lat = '+19:00:00'
    lon = '-155:00:00'
    elevation = 0.0 * u.m
    location = EarthLocation.from_geodetic(lon, lat, elevation)

    environment_kwargs = dict(pressure=1*u.bar, relative_humidity=0.1,
                              temperature=10*u.deg_C)

    obs1 = Observer(name='Observatory',
                    latitude=lat,
                    longitude=lon,
                    elevation=elevation,
                    **environment_kwargs)

    obs2 = Observer(name='Observatory',
                    location=location,
                    **environment_kwargs)

    assert obs1.location == obs2.location, ('using latitude/longitude/'
                                            'elevation keywords gave a '
                                            'different answer from passing in '
                                            'an EarthLocation directly')

@remote_data
def test_FixedTarget_from_name():
    '''
    Check that resolving target names with the `SkyCoord.from_name` constructor
    to produce a `FixedTarget` accurately resolves the coordinates of Polaris.
    '''

    # Resolve coordinates with SkyCoord.from_name classmethod
    polaris_from_name = FixedTarget.from_name('Polaris')
    # Coordinates grabbed from SIMBAD
    polaris_from_SIMBAD = SkyCoord('02h31m49.09456s', '+89d15m50.7923s')

    # Make sure separation is small
    assert polaris_from_name.coord.separation(polaris_from_SIMBAD) < 1*u.arcsec

def test_Observer_altaz():
    '''
    Check that the altitude/azimuth computed by `Observer.altaz` is similar
    to the result from PyEphem when pressure = 0 (no atmosphere) for Vega at
    2000-01-01 12:00:00 UTC.
    '''
    # Define the test case
    latitude = '00:00:00'
    longitude = '00:00:00'
    elevation = 0 # [m]
    pressure = 0  # no atmosphere
    time = Time('2000-01-01 12:00:00')
    vega_coords = SkyCoord('18h36m56.33635s', '+38d47m01.2802s')

    # Calculate altitude/azimuth with astroplan
    location = EarthLocation.from_geodetic(longitude, latitude, elevation*u.m)
    astroplan_obs = Observer(name='Observatory', location=location,
                             pressure=pressure*u.bar)
    astroplan_vega = FixedTarget(vega_coords)
    altaz = astroplan_obs.altaz(time, astroplan_vega)
    astroplan_altitude = altaz.alt
    astroplan_azimuth = altaz.az

    # Calculate altitude/azimuth with PyEphem, like so with get_pyephem_altaz
    '''
    >>> pyephem_altitude, pyephem_azimuth = get_pyephem_altaz(latitude,
                                                              longitude,
                                                              elevation,
                                                              time,
                                                              pressure,
                                                              vega_coords)
    >>> pyephem_altitude, pyephem_azimuth
    (<Latitude 51.198848716510874 deg>, <Longitude 358.4676707379987 deg>)
    '''
    pyephem_altitude = Latitude('51.198848716510874 deg')
    pyephem_azimuth = Longitude('358.4676707379987 deg')

    # Assert that altitudes/azimuths are within 30 arcsec - this is a wide
    # tolerance because the IERS tables used by astroplan may offset astroplan's
    # positions due to leap seconds.
    tolerance = (30*u.arcsec).to('deg').value
    assert_allclose(pyephem_altitude.value, astroplan_altitude.value,
                    atol=tolerance)
    assert_allclose(pyephem_azimuth.value, astroplan_azimuth.value,
                    atol=tolerance)

def get_pyephem_altaz(latitude, longitude, elevation, time, pressure,
                      target_coords):
    '''
    Run PyEphem to compute the altitude/azimuth of a target at specified time
    and observatory, for comparison with astroplan calucation tested in
    `test_Observer_altaz`.
    '''
    import ephem
    pyephem_obs = ephem.Observer()
    pyephem_obs.lat = latitude
    pyephem_obs.lon = longitude
    pyephem_obs.elevation = elevation
    pyephem_obs.date = time.datetime
    pyephem_obs.pressure = pressure
    pyephem_target = ephem.FixedBody()
    pyephem_target._ra = ephem.degrees(np.radians(target_coords.ra.value))
    pyephem_target._dec = ephem.degrees(np.radians(target_coords.dec.value))
    pyephem_target.compute(pyephem_obs)
    pyephem_altitude = Latitude(np.degrees(pyephem_target.alt)*u.degree)
    pyephem_azimuth = Longitude(np.degrees(pyephem_target.az)*u.degree)
    return pyephem_altitude, pyephem_azimuth

def test_Observer_timezone_parser():
    lat = '+19:00:00'
    lon = '-155:00:00'
    elevation = 0.0 * u.m
    location = EarthLocation.from_geodetic(lon, lat, elevation)

    obs1 = Observer(name='Observatory', location=location,
                    timezone=pytz.timezone('UTC'))
    obs2 = Observer(name='Observatory', location=location, timezone='UTC')
    obs3 = Observer(name='Observatory', location=location)

    assert obs1.timezone == obs2.timezone, ('Accept both strings to pass to '
                                            'the pytz.timezone() constructor '
                                            'and instances of pytz.timezone')

    assert obs2.timezone == obs3.timezone, ('Default timezone should be UTC')

def test_FixedTarget_ra_dec():
    '''
    Confirm that FixedTarget.ra and FixedTarget.dec are the same as the
    right ascension and declination stored in the FixedTarget.coord variable -
    which is a SkyCoord
    '''

    vega_coords = SkyCoord('18h36m56.33635s', '+38d47m01.2802s')
    vega = FixedTarget(vega_coords, name='Vega')
    assert vega.coord == vega_coords, 'Store coordinates directly'
    assert vega.coord.ra == vega_coords.ra == vega.ra, ('Retrieve RA from '
                                                        'SkyCoord')
    assert vega.coord.dec == vega_coords.dec == vega.dec, ('Retrieve Dec from '
                                                           'SkyCoord')

def test_sunrise_sunset_equator():
    '''
    Check that time of sunrise/set for an observer on the equator is
    consistent with PyEphem results (for no atmosphere/pressure=0)
    '''
    lat = '00:00:00'
    lon = '00:00:00'
    elevation = 0.0 * u.m
    pressure = 0
    location = EarthLocation.from_geodetic(lon, lat, elevation)
    time = Time('2000-01-01 12:00:00')
    obs = Observer(location=location, pressure=pressure)
    astroplan_next_sunrise = obs.sunrise(time, location, which='next').datetime
    astroplan_next_sunset = obs.sunset(time, location, which='next').datetime

    astroplan_prev_sunrise = obs.sunrise(time, location,
                                         which='previous').datetime
    astroplan_prev_sunset = obs.sunset(time, location,
                                       which='previous').datetime

    # Run get_pyephem_sunrise_sunset() to compute analogous
    # result from PyEphem:
    pyephem_next_sunrise = datetime.datetime(2000, 1, 2, 6, 3, 39, 150790)
    pyephem_next_sunset = datetime.datetime(2000, 1, 1, 18, 3, 23, 676686)
    pyephem_prev_sunrise = datetime.datetime(2000, 1, 1, 6, 3, 10, 720052)
    pyephem_prev_sunset = datetime.datetime(1999, 12, 31, 18, 2, 55, 100786)

    # Is there an equivalent to assert_allclose() for datetimes?
    #assert_allclose(pyephem_next_sunrise, pyephem_next_sunset,
    #                atol=datetime.timedelta(minutes=10))

    # Typical difference in this example between PyEphem and astroplan is <2 min
    threshold_minutes = 5
    assert (abs(pyephem_next_sunrise - astroplan_next_sunrise) <
            datetime.timedelta(minutes=threshold_minutes))
    assert (abs(pyephem_next_sunset - astroplan_next_sunset) <
            datetime.timedelta(minutes=threshold_minutes))
    assert (abs(pyephem_prev_sunrise - astroplan_prev_sunrise) <
            datetime.timedelta(minutes=threshold_minutes))
    assert (abs(pyephem_prev_sunset - astroplan_prev_sunset) <
            datetime.timedelta(minutes=threshold_minutes))

def get_pyephem_sunrise_sunset():
    '''
    Calculate next sunrise and sunset with PyEphem for an observer
    on the equator.
    '''
    lat = '00:00:00'
    lon = '00:00:00'
    elevation = 0.0 * u.m
    pressure = 0
    time = Time('2000-01-01 12:00:00')

    import ephem
    obs = ephem.Observer()
    obs.lat = lat
    obs.lon = lon
    obs.elevation = elevation
    obs.date = time.datetime
    obs.pressure = pressure
    next_sunrise = obs.next_rising(ephem.Sun(), use_center=True)
    next_sunset = obs.next_setting(ephem.Sun(), use_center=True)
    prev_sunrise = obs.previous_rising(ephem.Sun(), use_center=True)
    prev_sunset = obs.previous_setting(ephem.Sun(), use_center=True)

    return (repr(next_sunrise.datetime()),
            repr(next_sunset.datetime()),
            repr(prev_sunrise.datetime()),
            repr(prev_sunset.datetime()))

def test_vega_rise_set_equator():
    '''
    Check that time of rise/set of Vega for an observer on the equator is
    consistent with PyEphem results (for no atmosphere/pressure=0)
    '''
    lat = '00:00:00'
    lon = '00:00:00'
    elevation = 0.0 * u.m
    pressure = 0
    location = EarthLocation.from_geodetic(lon, lat, elevation)
    time = Time('2000-01-01 12:00:00')
    vega_ra, vega_dec = (279.23473479*u.degree, 38.78368896*u.degree)
    vega = SkyCoord(vega_ra, vega_dec)

    obs = Observer(location=location, pressure=pressure)
    astroplan_next_rise = obs.calc_rise(vega, time, location,
                                        which='next').datetime
    astroplan_next_set = obs.calc_set(vega, time, location,
                                      which='next').datetime

    astroplan_prev_rise = obs.calc_rise(vega, time, location,
                                           which='previous').datetime
    astroplan_prev_set = obs.calc_set(vega, time, location,
                                         which='previous').datetime

    # Run get_pyephem_vega_rise_set() to compute analogous
    # result from PyEphem:
    pyephem_next_rise, pyephem_next_set, pyephem_prev_rise, pyephem_prev_set = (
        datetime.datetime(2000, 1, 2, 5, 52, 8, 257401),
        datetime.datetime(2000, 1, 1, 17, 54, 6, 211705),
        datetime.datetime(2000, 1, 1, 5, 56, 4, 165852),
        datetime.datetime(1999, 12, 31, 17, 58, 2, 120088))

    # Typical difference in this example between PyEphem and astroplan is <2 min
    threshold_minutes = 5
    assert (abs(pyephem_next_rise - astroplan_next_rise) <
            datetime.timedelta(minutes=threshold_minutes))
    assert (abs(pyephem_next_set - astroplan_next_set) <
            datetime.timedelta(minutes=threshold_minutes))
    assert (abs(pyephem_prev_rise - astroplan_prev_rise) <
            datetime.timedelta(minutes=threshold_minutes))
    assert (abs(pyephem_prev_set - astroplan_prev_set) <
            datetime.timedelta(minutes=threshold_minutes))

def get_pyephem_vega_rise_set():
    '''
    Calculate next rise and set of Vega with PyEphem for an observer
    on the equator.
    '''
    lat = '00:00:00'
    lon = '00:00:00'
    elevation = 0.0 * u.m
    pressure = 0
    time = Time('2000-01-01 12:00:00')
    vega_ra, vega_dec = (279.23473479*u.degree, 38.78368896*u.degree)
    vega = SkyCoord(vega_ra, vega_dec)

    import ephem
    obs = ephem.Observer()
    obs.lat = lat
    obs.lon = lon
    obs.elevation = elevation
    obs.date = time.datetime
    obs.pressure = pressure
    target = ephem.FixedBody()
    target._ra = ephem.degrees(np.radians(vega.ra.value))
    target._dec = ephem.degrees(np.radians(vega.dec.value))
    target.compute(obs)
    next_rising = obs.next_rising(target).datetime()
    next_setting = obs.next_setting(target).datetime()
    prev_rising = obs.previous_rising(target).datetime()
    prev_setting = obs.previous_setting(target).datetime()

    return next_rising, next_setting, prev_rising, prev_setting

class TestRisingSetting(unittest.TestCase):
    def test_polaris_always_up_at_north_pole(self):
        with self.assertRaises(ValueError):
            lat = '90:00:00'
            lon = '00:00:00'
            elevation = 0.0 * u.m
            location = EarthLocation.from_geodetic(lon, lat, elevation)
            time = Time('2000-01-01 12:00:00')
            polaris = SkyCoord(37.95456067*u.degree, 89.26410897*u.degree)

            obs = Observer(location=location)
            astroplan_next_rise = obs.calc_rise(polaris, time, location,
                                                which='next').datetime
