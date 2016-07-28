# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This version of the `moon` module contains temporary solutions to calculating
the position of the moon using astropy.coordinates.get_moon, awaiting solutions
to the astropy issue [#5069](https://github.com/astropy/astropy/issues/5069).
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# Third-party
import numpy as np
from astropy.coordinates import (get_moon, get_sun, AltAz, Angle)

__all__ = ["moon_phase_angle", "moon_illumination"]


def moon_phase_angle(time, location, ephemeris=None):
    """
    Calculate lunar orbital phase [radians].

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time of observation

    location : `~astropy.coordinates.EarthLocation`
        Location of observer

    ephemeris : str, optional
        Ephemeris to use.  If not given, use the one set with
        ``astropy.coordinates.solar_system_ephemeris.set`` (which is
        set to 'builtin' by default).

    Returns
    -------
    i : float
        Phase angle of the moon [radians]
    """
    # TODO: cache these sun/moon SkyCoord objects

    # TODO: when astropy/astropy#5069 is resolved, replace this workaround which
    # handles scalar and non-scalar time inputs differently

    if time.isscalar:
        altaz_frame = AltAz(location=location, obstime=time)
        sun = get_sun(time).transform_to(altaz_frame)

        moon = get_moon(time, location=location, ephemeris=ephemeris).transform_to(altaz_frame)
        elongation = sun.separation(moon)
        return np.arctan2(sun.distance*np.sin(elongation),
                          moon.distance - sun.distance*np.cos(elongation))

    else:
        phase_angles = []
        for t in time:
            altaz_frame = AltAz(location=location, obstime=t)
            moon_coord = get_moon(t, location=location, ephemeris=ephemeris).transform_to(altaz_frame)
            sun_coord = get_sun(t).transform_to(altaz_frame)
            elongation = sun_coord.separation(moon_coord)
            phase_angle = np.arctan2(sun_coord.distance*np.sin(elongation),
                                     moon_coord.distance -
                                     sun_coord.distance*np.cos(elongation))
            phase_angles.append(phase_angle)
        return Angle(phase_angles)


def moon_illumination(time, location, ephemeris=None):
    """
    Calculate fraction of the moon illuminated

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time of observation

    location : `~astropy.coordinates.EarthLocation`
        Location of observer

    ephemeris : str, optional
        Ephemeris to use.  If not given, use the one set with
        ``astropy.coordinates.solar_system_ephemeris.set`` (which is
        set to 'builtin' by default).

    Returns
    -------
    k : float
        Fraction of moon illuminated
    """
    i = moon_phase_angle(time, location, ephemeris=ephemeris)
    k = (1 + np.cos(i))/2.0
    return k.value
