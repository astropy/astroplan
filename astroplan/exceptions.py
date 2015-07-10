# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

class AlwaysUpError(ValueError):
    '''Target is circumpolar'''
    pass

class NeverUpError(ValueError):
    '''Target never rises above horizon'''
    pass