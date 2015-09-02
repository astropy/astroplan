# this contains imports plugins that configure py.test for astropy tests.
# by importing them here in conftest.py they are discoverable by py.test
# no matter how it is invoked within the source tree.

from astropy.tests.pytest_plugins import *

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

# Make appropriate substitutions to mock internet querying methods
# within the tests
from .utils import _mock_remote_data
_mock_remote_data()

def pytest_configure(config):
    try:
        import matplotlib
        HAS_MATPLOTLIB = True
    except ImportError:
        HAS_MATPLOTLIB = False

    if HAS_MATPLOTLIB and config.pluginmanager.hasplugin('mpl'):
            config.option.mpl = True
            config.option.mpl_baseline_path = 'astroplan/plots/tests/baseline_images'
