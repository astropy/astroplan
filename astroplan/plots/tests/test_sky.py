# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from astropy.tests.helper import pytest

try:
    import matplotlib
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


# TODO: replace this with actual plot checks once these
# issues are resolved:
# https://github.com/astropy/astroplan/issues/65
# https://github.com/astropy/astroplan/issues/74
@pytest.mark.skipif('not HAS_MATPLOTLIB')
@pytest.mark.mpl_image_compare
def test_image_example():
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot([1, 2, 3])
    return fig
