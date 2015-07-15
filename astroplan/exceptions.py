# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from astropy.utils.exceptions import AstropyWarning

__all__ = ["TargetAlwaysUpWarning", "TargetNeverUpWarning"]

class TargetAlwaysUpWarning(AstropyWarning):
    """Target is circumpolar"""
    pass

class TargetNeverUpWarning(AstropyWarning):
    """Target never rises above horizon"""
    pass
