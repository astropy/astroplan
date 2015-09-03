# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from astropy.utils.exceptions import AstropyWarning

__all__ = ["TargetAlwaysUpWarning", "TargetNeverUpWarning",
           "OldEarthOrientationDataWarning", "PlotWarning",
           "PlotBelowHorizonWarning"]

class TargetAlwaysUpWarning(AstropyWarning):
    """Target is circumpolar"""
    pass

class TargetNeverUpWarning(AstropyWarning):
    """Target never rises above horizon"""
    pass

class OldEarthOrientationDataWarning(AstropyWarning):
    """Using old Earth rotation data from IERS"""
    pass

class PlotWarning(AstropyWarning):
    """Warnings dealing with the plotting aspects of astroplan"""
    pass

class PlotBelowHorizonWarning(PlotWarning):
    """Warning for when something is hidden on a plot because it's below the horizon"""
    pass
