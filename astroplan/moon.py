# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This version of the `moon` module contains temporary solutions to calculating
the position of the moon using PyEphem, awaiting solutions to the astropy issue
[#3920](https://github.com/astropy/astropy/issues/3920).
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# Third-party
import numpy as np
from astropy.time import Time
from astropy.coordinates import (SkyCoord, get_sun, AltAz, Angle)
import astropy.units as u
from astropy.utils.data import download_file
from astropy.coordinates.angle_utilities import angular_separation

__all__ = ["get_moon"]

def get_spk_file():
    """
    Get the Satellite Planet Kernel (SPK) file `de430.bsp` from NASA JPL.

    Download the file from the JPL webpage once and subsequently access a
    cached copy. This file is ~120 MB, and covers years ~1550-2650 CE [1]_.

    .. [1] http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/aareadme_de430-de431.txt
    """
    de430_url = ('http://naif.jpl.nasa.gov/pub/naif/'
                 'generic_kernels/spk/planets/de430.bsp')
    return download_file(de430_url, cache=True, show_progress=True)

def get_moon(time, location, pressure=None):
    """
    Position of the Earth's moon.

    Currently uses PyEphem to calculate the position of the moon by default
    (requires PyEphem to be installed). Set ``use_pyephem`` to `False` to
    calculate the moon position with jplephem (requires jplephem to be
    installed).

    Parameters
    ----------
    time : `~astropy.time.Time` or see below
        Time of observation. This will be passed in as the first argument to
        the `~astropy.time.Time` initializer, so it can be anything that
        `~astropy.time.Time` will accept (including a `~astropy.time.Time`
        object).

    location : `~astropy.coordinates.EarthLocation`
        Location of the observer on Earth

    pressure : `None` or `~astropy.units.Quantity` (optional)

    Returns
    -------
    moon_sc : `~astropy.coordinates.SkyCoord`
        Position of the moon at ``time``
    """
    if not isinstance(time, Time):
        time = Time(time)

    try:
        import ephem
    except ImportError:
        raise ImportError("The get_moon function currently requires "
                          "PyEphem to compute the position of the moon.")

    moon = ephem.Moon()
    obs = ephem.Observer()
    obs.lat = location.latitude.to(u.degree).to_string(sep=':')
    obs.lon = location.longitude.to(u.degree).to_string(sep=':')
    obs.elevation = location.height.to(u.m).value
    if pressure is not None:
        obs.pressure = pressure.to(u.bar).value*1000.0

    if time.isscalar:
        obs.date = time.datetime
        moon.compute(obs)
        moon_alt = float(moon.alt)
        moon_az = float(moon.az)
        moon_distance = moon.earth_distance
    else:
        moon_alt = []
        moon_az = []
        moon_distance = []
        for t in time:
            obs.date = t.datetime
            moon.compute(obs)
            moon_alt.append(float(moon.alt))
            moon_az.append(float(moon.az))
            moon_distance.append(moon.earth_distance)
    return SkyCoord(alt=moon_alt*u.rad, az=moon_az*u.rad,
                    distance=moon_distance*u.AU,
                    frame=AltAz(location=location, obstime=time))

def moon_phase_angle(time, location):
    """
    Calculate lunar orbital phase [radians].

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time of observation

    location : `~astropy.coordinates.EarthLocation`
        Location of observer

    Returns
    -------
    i : float
        Phase angle of the moon [radians]
    """
    # TODO: cache these sun/moon SkyCoord objects
    sun = get_sun(time).transform_to(AltAz(location=location, obstime=time))
    moon = get_moon(time, location)
    # The line below should have worked, but needs a workaround.
    # TODO: once bug has been fixed, replace workaround with simpler version.
    # elongation = sun.separation(moon)
    elongation = Angle(angular_separation(moon.spherical.lon,
                                   moon.spherical.lat,
                                   sun.spherical.lon,
                                   sun.spherical.lat), u.deg)
    return np.arctan2(sun.distance*np.sin(elongation),
                      moon.distance - sun.distance*np.cos(elongation))

def moon_illumination(time, location):
    """
    Calculate fraction of the moon illuminated

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time of observation

    location : `~astropy.coordinates.EarthLocation`
        Location of observer

    Returns
    -------
    k : float
        Fraction of moon illuminated
    """
    i = moon_phase_angle(time, location)
    k = (1 + np.cos(i))/2.0
    return k.value
