from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from .. import sites
from ..sites import get_site, add_site, new_site_info_to_json
from ..observer import Observer
from astropy.tests.helper import assert_quantity_allclose
from astropy.coordinates import Latitude, Longitude, EarthLocation
import astropy.units as u
from astropy.tests.helper import pytest
import json

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
    del sites._site_names[sites._site_names.index(new_site_name)]
    del sites._site_db[new_site_name.lower()]

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
    lon_str = "-155d28m34s"
    lat_str = "+19d49m32s"
    elevation = 4139*u.m
    location = EarthLocation.from_geodetic(lon_str, lat_str, elevation)
    short_name = "New telescope (subaru)"
    aliases = ["example new telescope with subaru's coordinates"]
    source = "the tests module"
    new_site_json = new_site_info_to_json(short_name, location, aliases, source)
    new_site = json.loads(new_site_json)

    ns = new_site[short_name]
    assert_quantity_allclose(Longitude(lon_str),
                             Longitude(ns["longitude"]*u.Unit(ns["longitude_unit"])),
                             atol=0.001*u.deg)
    assert_quantity_allclose(Latitude(lat_str),
                             Latitude(ns["latitude"]*u.Unit(ns["latitude_unit"])),
                             atol=0.001*u.deg)
    assert_quantity_allclose(elevation,
                             new_site[short_name]["elevation"]*u.Unit(ns["elevation_unit"]),
                             atol=1*u.m)
    assert short_name == new_site[short_name]['name']
    assert aliases == new_site[short_name]['aliases']

    with pytest.raises(ValueError):
        # This name already exists
        new_site_info_to_json("Keck", location, aliases, source)
