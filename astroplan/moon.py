# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This version of the `moon` module contains temporary solutions to calculating
the position of the moon using PyEphem, awaiting solutions to the astropy issue
[#3920](https://github.com/astropy/astropy/issues/3920).
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

try:
    import ephem
except ImportError:
    raise ImportError("The astroplan.moon module currently requires "
                      "PyEphem to compute the position of the moon.")

import numpy as np

def calc_moon_phase_angle(moon, sun):
    """
    Calculate lunar orbital phase [radians].

    Parameters
    ----------
    moon : `~astropy.coordinates.SkyCoord`
        Position of the moon

    sun : `~astropy.coordinates.SkyCoord`
        Position of the Sun

    Returns
    -------
    i : float
        Phase angle of the moon [radians]
    """
    elongation = sun.separation(moon)
    return np.arctan2(sun.distance*np.sin(elongation),
                      moon.distance - sun.distance*np.cos(elongation))

def calc_moon_illumination(moon, sun):
    """
    Calculate fraction of the moon illuminated

    Parameters
    ----------
    moon : `~astropy.coordinates.SkyCoord`
        Position of the moon

    sun : `~astropy.coordinates.SkyCoord`
        Position of the Sun

    Returns
    -------
    k : float
        Fraction of moon illuminated
    """
    i = calc_moon_phase_angle(moon, sun)
    k = (1 + np.cos(i))/2.0
    return k.value
