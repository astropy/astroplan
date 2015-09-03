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




def pytest_configure(config):
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
    # Make appropriate substitutions to mock internet querying methods
    # within the tests
    if item.get_marker('remote_data'):
        1/0
    else:
        2/0
    _mock_remote_data()

def pytest_runtest_teardown(item, nextitem):
    #
    _unmock_remote_data()


def _mock_remote_data():
    # Mock FixedTarget.from_name class method for tests without remote data
    from .target import FixedTarget

    FixedTarget._real_from_name = FixedTarget.from_name
    FixedTarget.from_name = FixedTarget._from_name_mock

def _unmock_remote_data():
    # undo _mock_remote_data
    from .target import FixedTarget

    if hasattr(FixedTarget, '_real_from_name'):
        FixedTarget.from_name = FixedTarget._real_from_name
        del FixedTarget._real_from_name
    #otherwise assume it's already correct
