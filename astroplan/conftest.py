# this contains imports plugins that configure py.test for astropy tests.
# by importing them here in conftest.py they are discoverable by py.test
# no matter how it is invoked within the source tree.

from astropy.tests.pytest_plugins import *

# also save a copy of the astropy hooks so we can use them below when overriding
from astropy.tests import pytest_plugins as astropy_pytest_plugins

import warnings
from .utils import _mock_remote_data, _unmock_remote_data
from .exceptions import OldEarthOrientationDataWarning

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

    warnings.simplefilter('always', category=OldEarthOrientationDataWarning)

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
