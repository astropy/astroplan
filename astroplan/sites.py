# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Observatories accessible by the `sites` module originate from the IRAF
Observatory Database, and are stored in astroplan/data/observatories.json.
Longitudes are listed with positive to the West.

Additions and corrections to the observatory list can be submitted via Pull
Request to the [astroplan GitHub repository](https://github.com/astropy/astroplan).

"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# Standard library
from difflib import get_close_matches
import json

# Third-party
from astropy.utils.data import get_pkg_data_contents
from astropy.coordinates import EarthLocation
import astropy.units as u

__all__ = ['get_site', 'get_site_names', 'add_site']

# Observatory database and list of names:
_site_db = None
_site_names = []


def _load_sites():
    """
    Load observatory database from astroplan/data/observatories.json
    """
    global _site_db, _site_names
    _site_db = dict()
    db = json.loads(get_pkg_data_contents('data/observatories.json'))
    for site in db:
        location = EarthLocation.from_geodetic(db[site]['longitude'],
                                   db[site]['latitude'],
                                   db[site]['elevation_meters'])
        _site_names.append(db[site]['name'])
        _site_names.append(site)

        _site_db[site.lower()] = location
        _site_db[db[site]['name'].lower()] = location
        for alias in db[site]['aliases']:
            _site_db[alias.lower()] = location


def get_site(site_name):
    """
    Return an `~astropy.coordinates.EarthLocation` object for known observatory.

    Use `~astroplan.get_site_names` to see a list of available
    observatories.

    Parameters
    ----------
    site_name : str
        Name of the observatory.

    Returns
    -------
    `~astropy.coordinates.EarthLocation`
        The location of the observatory.
    """
    if _site_db is None:
        _load_sites()

    if site_name.lower() not in _site_db.keys():
        # If site name not found, find close matches and suggest them in error
        close_names = get_close_matches(site_name, _site_db.keys())
        close_names = sorted(close_names, key=lambda x: len(x))
        if len(close_names) > 0:
            errmsg = ("Site not in database. Use `astroplan.get_site_names()` "
                      "to see available sites. Did you mean: '{}'?".format(
                      "', '".join(close_names)))
        else:
            errmsg = 'Site not in database.'
        raise KeyError(errmsg)

    return _site_db[site_name.lower()]


def get_site_names(full_list=True):
    """
    Get list of names of observatories for use with `~astroplan.get_site`.

    Parameters
    ----------
    full_list : bool
        Show full list observatory names and aliases (True), or just the list
        of names (False)? Default to True.

    Returns
    -------
    List of observatory names (strings)
    """
    if _site_db is None:
        _load_sites()

    if full_list:
        return sorted(_site_db.keys())
    else:
        return sorted(_site_names)

def add_site(site_name, location):
    """
    Add a site to the list of available observatories.

    Parameters
    ----------
    site_name : string
        Name of the observatory

    location : `~astropy.coordinates.EarthLocation`
        Location of the observatory
    """
    if _site_db is None:
        _load_sites()

    if not isinstance(location, EarthLocation):
        raise ValueError('Location must be an EarthLocation.')

    if site_name.lower() not in _site_db.keys():
        _site_db[site_name.lower()] = location
        _site_names.append(site_name)
    else:
        raise KeyError('The site "{}" already exists at (longitude,latitude,'
                       'elevation)={}'.format(site_name,
                                     _site_db[site_name.lower()].to_geodetic()))


def new_site_info_to_json(short_name, location, aliases, source):
    """
    Generate JSON-formatted observatory information string.

    Use this function to prepare new observatory submissions for the astroplan
    observatory database. To submit the observatory for inclusion in the next
    version of astroplan, post a pull request on the astroplan GitHub
    repository [1]_. The json string should be put into the database in
    `astroplan/data/observatories.json`.

    .. [1] https://github.com/astroplanners/astroplan/pulls

    Parameters
    ----------
    short_name : str
        Short name of observatory

    location : `~astropy.coordinates.EarthLocation`
        Location of the observatory

    elevation : `~astropy.units.Quantity`
        Elevation of the observatory

    aliases : list of strings
        List of common abbreviations for the observatory name

    source : str
        Source of observatory information

    Return
    ------
    json_str : str
        String representation of the JSON-formatted observatory information
        for submissions to the astroplan observatories database via pull request

    Example
    -------
    Make a new site for a telescope at the University of Washington
    >>> from astropy.coordinates import EarthLocation
    >>> from astroplan.sites import new_site_info_to_json
    >>> import astropy.units as u
    >>> location = EarthLocation.from_geodetic(-122.307703*u.deg, 47.653706*u.deg, 10*u.m)
    >>> short_name = "My Fake Telescope"
    >>> aliases = ["UW Drumheller Fountain Telescope"]
    >>> source = "Brett Morris, personal communication"
    >>> print(new_site_info_to_json(short_name, location, aliases, source)) # doctest: +SKIP
    {
        "My Fake Telescope": {
            "aliases": [
                "UW Drumheller Fountain Telescope"
            ],
            "elevation_meters": 9.999999999832553,
            "latitude": "47d39m13.3416s",
            "longitude": "-122d18m27.7308s",
            "name": "My Fake Telescope",
            "source": "Brett Morris, personal communication"
        }
    }

    """
    if _site_db is None:
        _load_sites()

    # Gather all site names and aliases currently in database
    all_site_names_and_aliases = _site_db.keys()

    # Raise error if name/alias collision occurs
    if short_name.lower() in all_site_names_and_aliases:
        raise ValueError('short_name "{}" is already registered in the '
                         'observatory database.'.format(short_name))
    for alias in aliases:
        if alias.lower() in all_site_names_and_aliases:
            raise ValueError('alias {} is already registered in the '
                             'observatory database.'.format(alias))

    output_dict = {short_name :
                        dict(name=short_name,
                             longitude=location.longitude.to_string(unit=u.deg),
                             latitude=location.latitude.to_string(unit=u.deg),
                             elevation_meters=location.height.to(u.m).value,
                             aliases=aliases,
                             source=source)}
    return json.dumps(output_dict, indent=4, sort_keys=True)
