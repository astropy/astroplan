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
from astropy.tests.pytest_plugins import *

# also save a copy of the astropy hooks so we can use them below when overriding
from astropy.tests import pytest_plugins as astropy_pytest_plugins

import warnings
from .utils import _mock_remote_data, _unmock_remote_data
from .exceptions import AstroplanWarning

## Uncomment the following line to treat all DeprecationWarnings as
## exceptions
enable_deprecations_as_exceptions()

## Uncomment and customize the following lines to add/remove entries
## from the list of packages for which version numbers are displayed
## when running the tests
try:
    PYTEST_HEADER_MODULES['pyephem'] = 'pyephem'
    del PYTEST_HEADER_MODULES['h5py']
except NameError:  # needed to support Astropy < 1.0
    pass


def pytest_configure(config):
    if hasattr(astropy_pytest_plugins, 'pytest_configure'):
        # sure ought to be true right now, but always possible it will change in
        # future versions of astropy
        astropy_pytest_plugins.pytest_configure(config)

    #make sure astroplan warnings always appear so we can test when they show up
    warnings.simplefilter('always', category=AstroplanWarning)

    try:
        import matplotlib
        import nose  # needed for the matplotlib testing tools
        HAS_MATPLOTLIB_AND_NOSE = True
    except ImportError:
        HAS_MATPLOTLIB_AND_NOSE = False

    if HAS_MATPLOTLIB_AND_NOSE and config.pluginmanager.hasplugin('mpl'):
            config.option.mpl = True
            config.option.mpl_baseline_path = 'astroplan/plots/tests/baseline_images'

def pytest_runtest_setup(item):
    """
    This overrides the tests so that if they are marked ``remote_data`` they get
    run without any mocking of functions, but if they are not, then the mocking
    happens.  This means that functionality that uses mock data should have both
    a ``remote_data`` test *and* a separate one that is not ``remote_data``.
    """
    if hasattr(astropy_pytest_plugins, 'pytest_runtest_setup'):
        # sure ought to be true right now, but always possible it will change in
        # future versions of astropy
        astropy_pytest_plugins.pytest_runtest_setup(item)

    # Make appropriate substitutions to mock internet querying methods
    # within the tests
    if item.get_marker('remote_data'):
        # a no-op if the last one was remote-data
        _unmock_remote_data()
    else:
        # a no-op if the last one was not remote-data
        _mock_remote_data()
