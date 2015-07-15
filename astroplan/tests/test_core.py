from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from astropy.coordinates import (EarthLocation, Latitude, Longitude, SkyCoord)
import astropy.units as u
from astropy.time import Time
from astropy.tests.helper import remote_data
import numpy as np
from numpy.testing import assert_allclose
import pytz

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

def test_parallactic_angle():
    '''
    Compute parallactic angle for targets at hour angle = {3, 19} for
    at observer at IRTF using the online SpeX calculator and PyEphem
    '''
    # Set up position for IRTF
    lat = 19.826218*u.deg
    lon = -155.471999*u.deg
    elevation = 4160.0 * u.m
    location = EarthLocation.from_geodetic(lon, lat, elevation)
    time = Time('2015-01-01 00:00:00')
    LST = time.sidereal_time('mean', longitude=lon)
    desired_HA_1 = 3*u.hourangle
    desired_HA_2 = 19*u.hourangle # = -5*u.hourangle
    target1 = SkyCoord(LST - desired_HA_1, -30*u.degree)
    target2 = SkyCoord(LST - desired_HA_2, -30*u.degree)

    obs = Observer(location=location)
    q1 = obs.parallactic_angle(time, target1)
    q2 = obs.parallactic_angle(time, target2)

    # Get values from PyEphem for comparison from print_pyephem_parallactic_angle()
    pyephem_q1 = 46.54610060782033*u.deg
    pyephem_q2 = -65.51818282032019*u.deg

    assert_allclose(q1.to(u.degree).value, pyephem_q1, atol=1)
    assert_allclose(q2.to(u.degree).value, pyephem_q2, atol=1)

    # Get SpeX parallactic angle calculator values for comparison from
    # http://irtfweb.ifa.hawaii.edu/cgi-bin/spex/parangle.cgi to produce

    SpeX_q1 = 46.7237968 # deg
    SpeX_q2 = -65.428924 # deg

    assert_allclose(q1.to(u.degree).value, SpeX_q1, atol=0.1)
    assert_allclose(q2.to(u.degree).value, SpeX_q2, atol=0.1)

def print_pyephem_parallactic_angle():
    lat = 19.826218*u.deg
    lon = -155.471999*u.deg
    time = Time('2015-01-01 00:00:00')
    LST = time.sidereal_time('mean', longitude=lon)
    desired_HA_1 = 3*u.hourangle
    desired_HA_2 = 19*u.hourangle # = -5*u.hourangle

    import ephem
    obs = ephem.Observer()
    obs.lat = '19:49:34.3848'
    obs.lon = '-155:28:19.1964'
    obs.elevation = elevation.value
    obs.date = time.datetime
    pyephem_target1 = ephem.FixedBody()
    pyephem_target1._ra = ephem.degrees((LST - desired_HA_1).to(u.rad).value)
    pyephem_target1._dec = ephem.degrees((-30*u.deg).to(u.rad).value)
    pyephem_target1.compute(obs)
    pyephem_q1 = (float(pyephem_target1.parallactic_angle())*u.rad).to(u.deg)

    pyephem_target2 = ephem.FixedBody()
    pyephem_target2._ra = ephem.degrees((LST - desired_HA_2).to(u.rad).value)
    pyephem_target2._dec = ephem.degrees((-30*u.deg).to(u.rad).value)
    pyephem_target2.compute(obs)
    pyephem_q2 = (float(pyephem_target2.parallactic_angle())*u.rad).to(u.deg)
    print(pyephem_q1, pyephem_q2)