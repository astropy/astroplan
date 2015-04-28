# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from astropy.time import Time
from ..utils import observe


def test_observe():
    time, pos = observe()
    millenium = Time('2000-01-01T00:00:00')
    assert time > millenium
