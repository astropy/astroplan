from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from ..sites import get_site, add_site, new_site_info_to_json
from ..core import Observer
from astropy.tests.helper import assert_quantity_allclose
from astropy.coordinates import Latitude, Longitude, EarthLocation
import astropy.units as u
import pytest

def test_get_site():
    # Compare to the IRAF observatory list available at:
    # http://tdc-www.harvard.edu/iraf/rvsao/bcvcorr/obsdb.html
    keck = get_site('keck')
    lon, lat, el = keck.to_geodetic()
    assert_quantity_allclose(lon, -1*Longitude('155:28.7', unit=u.deg),
                             atol=0.001*u.deg)
    assert_quantity_allclose(lat, Latitude('19:49.7', unit=u.deg),
                             atol=0.001*u.deg)
    assert_quantity_allclose(el, 4160*u.m, atol=1*u.m)

    keck = get_site('ctio')
    lon, lat, el = keck.to_geodetic()
    assert_quantity_allclose(lon, -1*Longitude('70.815', unit=u.deg),
                             atol=0.001*u.deg)
    assert_quantity_allclose(lat, Latitude('-30.16527778', unit=u.deg),
                             atol=0.001*u.deg)
    assert_quantity_allclose(el, 2215*u.m, atol=1*u.m)

    subaru = get_site('subaru')
    lon, lat, el = subaru.to_geodetic()
    assert_quantity_allclose(lon, -1*Longitude("155d28m34"), atol=0.001*u.deg)

def test_add_site():
    # Test observatory can be added and retrieved
    new_site_name = 'University of Washington'
    new_site_location = EarthLocation(-122.3080*u.deg, 47.6550*u.deg, 0*u.m)
    add_site(new_site_name, new_site_location)
    retrieved_location = get_site(new_site_name)
    assert retrieved_location == new_site_location

def test_Observer_classmethod():
    site_name = 'kpno'
    pressure = 1*u.bar
    temperature = 0*u.deg_C
    kpno = Observer.at_site(site_name, pressure=pressure,
                            temperature=temperature)
    assert kpno.location == get_site(site_name)
    assert kpno.name == site_name
    assert kpno.temperature == temperature
    assert kpno.pressure == pressure

def test_bad_site():
    with pytest.raises(KeyError):
        get_site('nonexistent site')

def test_new_site_info_to_json():
    location = EarthLocation.from_geodetic("-155d28m34s", "+19d49m32s", 4139*u.m)
    short_name = "New telescope (subaru)"
    aliases = ["example new telescope with subaru's coordinates"]
    source = "the tests module"
    new_site_json = new_site_info_to_json(short_name, location, aliases, source)
    assert new_site_json == ('{\n    "New telescope (subaru)": {\n        '
                             '"aliases": [\n            "example new telescope '
                             'with subaru\'s coordinates"\n        ], \n'
                             '        "elevation_meters": 4139.000000000389,'
                             ' \n        "latitude": "19d49m32s", \n        '
                             '"longitude": "-155d28m34s", \n        "name": '
                             '"New telescope (subaru)", \n        "source": '
                             '"the tests module"\n    }\n}')

    with pytest.raises(ValueError):
        # This name already exists
        new_site_info_to_json("Keck", location, aliases, source)

class TestExceptions(unittest.TestCase):
    def test_bad_site(self):
        with self.assertRaises(KeyError):
            get_site('nonexistent site')
