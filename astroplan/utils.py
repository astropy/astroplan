# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from astropy.time import Time
from astropy.coordinates import get_sun

__all__ = ['observe']


def observe():
    """Dummy function to check if tests / docs work.
    """
    time = Time.now()
    pos = get_sun(time)
    print('Current time: {}'.format(time))
    print('Sun location: {}'.format(pos))
    return time, pos
