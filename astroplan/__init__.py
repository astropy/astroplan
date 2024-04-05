# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
astroplan is an open source (BSD licensed) observation planning package for
astronomers that can help you plan for everything but the clouds.

It is an in-development `Astropy <http://www.astropy.org>`__
`affiliated package <http://www.astropy.org/affiliated/index.html>`__ that
seeks to make your life as an observational astronomer a little less
infuriating.

* Code: https://github.com/astropy/astroplan
* Docs: https://astroplan.readthedocs.io/
"""

try:
    from .version import version as __version__
except ImportError:
    __version__ = ''

from .utils import *
from .observer import *
from .target import *
from .exceptions import *
from .moon import *
from .constraints import *
from .scheduling import *
from .periodic import *

download_IERS_A()
