# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This version of the `moon` module calculates lunar phase angle for a geocentric
"""

# Third-party
import numpy as np
from astropy.coordinates import get_sun, get_body, GCRS, SkyCoord

__all__ = ["moon_phase_angle", "moon_illumination"]


def _geocentric(coord):
    """
    Shift a topocentric GCRS coordinate (as returned by
    `~astropy.coordinates.get_body` with a ``location``) back to the geocentre.
    """
    return SkyCoord(GCRS(coord.cartesian + coord.frame.obsgeoloc,
                         obstime=coord.obstime))


def moon_phase_angle(time, ephemeris=None, observer=None):
    """
    Calculate lunar orbital phase in radians.

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time of observation

    ephemeris : str, optional
        Ephemeris to use.  If not given, use the one set with
        `~astropy.coordinates.solar_system_ephemeris` (which is
        set to 'builtin' by default).

    observer : `~astroplan.Observer`, optional
        If given, the Moon position (and the Sun's, if already cached) is
        taken from the observer's cache (`~astroplan.Observer.get_body`) and
        shifted to the geocentre, avoiding a recomputation. The phase angle
        is geocentric either way.

    Returns
    -------
    i : `~astropy.units.Quantity`
        Phase angle of the moon [radians]
    """
    if observer is None:
        sun = get_sun(time)
        moon = get_body("moon", time, ephemeris=ephemeris)
    else:
        moon = _geocentric(observer.get_body("moon", time, ephemeris=ephemeris))
        sun = observer.get_body("sun", time, compute=False)
        sun = get_sun(time) if sun is None else _geocentric(sun)
    elongation = sun.separation(moon)
    return np.arctan2(sun.distance*np.sin(elongation),
                      moon.distance - sun.distance*np.cos(elongation))


def moon_illumination(time, ephemeris=None, observer=None):
    """
    Calculate fraction of the moon illuminated.

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time of observation

    ephemeris : str, optional
        Ephemeris to use.  If not given, use the one set with
        `~astropy.coordinates.solar_system_ephemeris` (which is
        set to 'builtin' by default).

    observer : `~astroplan.Observer`, optional
        If given, use the observer's cached positions (see
        `moon_phase_angle`).

    Returns
    -------
    k : float
        Fraction of moon illuminated
    """
    i = moon_phase_angle(time, ephemeris=ephemeris, observer=observer)
    k = (1 + np.cos(i))/2.0
    return k.value
