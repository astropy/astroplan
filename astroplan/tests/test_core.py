from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from astropy.coordinates import (EarthLocation, Latitude, Longitude, SkyCoord)
import astropy.units as u
from astropy.time import Time
from astropy.tests.helper import remote_data
import numpy as np
import ephem
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

    assert obs1.location == obs2.location, ('Locations not parsed equivalently')

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
    altaz = astroplan_obs.altaz(astroplan_vega, time)
    astroplan_altitude = altaz.alt
    astroplan_azimuth = altaz.az

    # Calculate altitude/azimuth with PyEphem
    pyephem_obs = ephem.Observer()
    pyephem_obs.lat = latitude
    pyephem_obs.lon = longitude
    pyephem_obs.elevation = elevation
    pyephem_obs.date = time.datetime
    pyephem_obs.pressure = pressure
    pyephem_vega = ephem.FixedBody()
    pyephem_vega._ra = ephem.degrees(np.radians(vega_coords.ra.value))
    pyephem_vega._dec = ephem.degrees(np.radians(vega_coords.dec.value))
    pyephem_vega.compute(pyephem_obs)
    pyephem_altitude = Latitude(np.degrees(pyephem_vega.alt)*u.degree)
    pyephem_azimuth = Longitude(np.degrees(pyephem_vega.az)*u.degree)

    # Assert that altitudes/azimuths are within 30 arcsec - this is a wide
    # tolerance because the IERS tables used by astroplan may offset astroplan's
    # positions due to leap seconds.
    tolerance = 30*u.arcsec
    assert (abs(pyephem_altitude - astroplan_altitude) < tolerance and
            abs(pyephem_azimuth - astroplan_azimuth) < tolerance)

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