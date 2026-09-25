# Use a non-interactive backend so tests and doctests that call plt.show()
# don't open windows (or block) when run outside of tox.
try:
    import matplotlib
    matplotlib.use('Agg')
except ImportError:
    pass

try:
    from pytest_astropy_header.display import PYTEST_HEADER_MODULES, TESTED_VERSIONS
except ImportError:  # In case this plugin is not installed
    PYTEST_HEADER_MODULES = {}
    TESTED_VERSIONS = {}

# This is to figure out the affiliated package version, rather than
# using Astropy's
try:
    from astroplan import __version__
except ImportError:
    __version__ = 'dev'

TESTED_VERSIONS["astroplan"] = __version__

# Define list of packages for which to display version numbers in the test log
PYTEST_HEADER_MODULES['astropy'] = 'astropy'
PYTEST_HEADER_MODULES['tzdata'] = 'tzdata'
PYTEST_HEADER_MODULES['astroquery'] = 'astroquery'
PYTEST_HEADER_MODULES['pyvo'] = 'pyvo'
if "h5py" in PYTEST_HEADER_MODULES:
    del PYTEST_HEADER_MODULES['h5py']
