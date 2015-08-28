# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import pytest
import matplotlib.pyplot as plt


# TODO: replace this with actual plot checks once these
# issues are resolved:
# https://github.com/astropy/astroplan/issues/65
# https://github.com/astropy/astroplan/issues/74
@pytest.mark.mpl_image_compare
def test_image_example():
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot([1, 2, 3])
    return fig
