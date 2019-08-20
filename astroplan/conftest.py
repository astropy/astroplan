"""
This file contains functions that configure py.test like astropy but with
additions for astroplan.  Py.test looks for specially-named functions
(like  ``pytest_configure``) and uses those to configure itself.

Here, we want to keep the behavior of astropy while *adding* more for astroplan.
To do that, in the functions below, we first invoke the functions from astropy,
and then after that do things specific to astroplan.  But we also want astropy
functionality for any functions we have *not* overriden, so that's why the
``import *`` happens at the top.
"""
from astropy.tests.plugins.display import (TESTED_VERSIONS,
                                           PYTEST_HEADER_MODULES)

from astropy.tests.helper import enable_deprecations_as_exceptions

# We do this to pick up the test header report even when using LTS astropy
try:
    from astropy.tests.pytest_plugins import pytest_report_header  # noqa
except ImportError:
    pass

import os

# This is to figure out the affiliated package version, rather than
# using Astropy's
try:
    from .version import version
except ImportError:
    version = 'dev'

packagename = os.path.basename(os.path.dirname(__file__))
TESTED_VERSIONS[packagename] = version


# Comment out this line to avoid deprecation warnings being raised as
# exceptions
enable_deprecations_as_exceptions()

# Define list of packages for which to display version numbers in the test log
try:
    PYTEST_HEADER_MODULES['Astropy'] = 'astropy'
    PYTEST_HEADER_MODULES['pytz'] = 'pytz'
    PYTEST_HEADER_MODULES['pyephem'] = 'ephem'
    PYTEST_HEADER_MODULES['matplotlib'] = 'matplotlib'
    PYTEST_HEADER_MODULES['pytest-mpl'] = 'pytest_mpl'
    del PYTEST_HEADER_MODULES['h5py']
except KeyError:
    pass
