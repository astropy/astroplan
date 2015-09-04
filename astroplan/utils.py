# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# Standard library
import warnings

# Third-party
import numpy as np
from astropy.utils.data import download_file, clear_download_cache
from astropy.utils import iers
from astropy.time import Time
import astropy.units as u
from astropy.utils.data import (_get_download_cache_locs, CacheMissingWarning,
                                _open_shelve)

# Package
from .exceptions import OldEarthOrientationDataWarning

__all__ = ["get_IERS_A_or_workaround", "download_IERS_A",
           "time_grid_from_range"]

IERS_A_WARNING = ("For best precision (on the order of arcseconds), you must "
                  "download an up-to-date IERS Bulletin A table. To do so, run:"
                  "\n\n"
                  ">>> from astroplan import download_IERS_A\n"
                  ">>> download_IERS_A()\n")

BACKUP_Time_get_delta_ut1_utc = Time._get_delta_ut1_utc


def _low_precision_utc_to_ut1(self, jd1, jd2):
    """
    When no IERS Bulletin A is available (no internet connection), use low
    precision time conversion by assuming UT1-UTC=0 always.
    This method mimics `~astropy.coordinates.builtin_frames.utils.get_dut1utc`
    """
    try:
        return self.delta_ut1_utc
    except IndexError:
        warnings.warn(IERS_A_WARNING, OldEarthOrientationDataWarning)
        return np.zeros(self.shape)


def get_IERS_A_or_workaround():
    """
    Get the cached IERS Bulletin A table if one exists. If one does not exist,
    monkey patch `~astropy.time.Time._get_delta_ut1_utc` so that
    `~astropy.time.Time` objects don't raise errors by computing UT1-UTC off
    the end of the IERS table.
    """
    if IERS_A_in_cache():
        iers.IERS.iers_table = _get_IERS_A_table()
    else:
        Time._get_delta_ut1_utc = _low_precision_utc_to_ut1


def IERS_A_in_cache():
    """
    Check if the IERS Bulletin A table is locally cached.
    """
    url_key = iers.IERS_A_URL
    # The below code which accesses ``urlmapfn`` is stolen from
    # astropy.utils.data.download_file()
    try:
        dldir, urlmapfn = _get_download_cache_locs()
    except (IOError, OSError) as e:
        msg = 'Remote data cache could not be accessed due to '
        estr = '' if len(e.args) < 1 else (': ' + str(e))
        warnings.warn(CacheMissingWarning(msg + e.__class__.__name__ + estr))
        return False
    with _open_shelve(urlmapfn, True) as url2hash:
        # TODO: try to figure out how to test this in the unicode case
        if str(url_key) in url2hash:
            return True
    return False


def _get_IERS_A_table(warn_update=14*u.day):
    """
    Grab the locally cached copy of the IERS Bulletin A table. Check to see
    if it's up to date, and warn the user if it is not.

    This will fail and raise OSError if the file is not in the cache.
    """
    if IERS_A_in_cache():
        path = download_file(iers.IERS_A_URL, cache=True, show_progress=True)
        table = iers.IERS_A.open(path)
        # Use polar motion flag to identify last observation before predictions
        index_of_last_observation = ''.join(table['PolPMFlag_A']).index('IP')
        time_of_last_observation = Time(table['MJD'][index_of_last_observation],
                                        format='mjd')
        time_since_last_update = Time.now() - time_of_last_observation

        # If the IERS bulletin is more than `warn_update` days old, warn user
        if warn_update < time_since_last_update:
            warnmsg = ("Your version of the IERS Bulletin A is {:.1f} days "
                       "old. ".format(time_since_last_update.to(u.day).value) +
                       IERS_A_WARNING)
            warnings.warn(warnmsg, OldEarthOrientationDataWarning)
        return table
    else:
        raise OSError("No IERS A table has been downloaded.")


def download_IERS_A(show_progress=True):
    """
    Download and cache the IERS Bulletin A table.

    If one is already cached, download a new one and overwrite the old. Store
    table in the astropy cache, and undo the monkey patching done by
    `~astroplan.get_IERS_A_or_workaround`.

    Parameters
    ----------
    show_progress : bool
        `True` shows a progress bar during the download.
    """
    if IERS_A_in_cache():
        clear_download_cache(iers.IERS_A_URL)

    local_iers_a_path = download_file(iers.IERS_A_URL, cache=True,
                                      show_progress=show_progress)
    # Undo monkey patch set up by get_IERS_A_or_workaround
    iers.IERS.iers_table = iers.IERS_A.open(local_iers_a_path)
    Time._get_delta_ut1_utc = BACKUP_Time_get_delta_ut1_utc


@u.quantity_input(time_resolution=u.hour)
def time_grid_from_range(time_range, time_resolution=0.5*u.hour):
    """
    Get linearly-spaced sequence of times.

    Parameters
    ----------
    time_range : `~astropy.time.Time` (length = 2)
        Lower and upper bounds on time sequence.

    time_resolution : `~astropy.units.quantity` (optional)
        Time-grid spacing

    Returns
    -------
    times : `~astropy.time.Time`
        Linearly-spaced sequence of times
    """
    return Time(np.arange(time_range[0].jd, time_range[1].jd,
                          time_resolution.to(u.day).value), format='jd')

def _mock_remote_data():
    """
    Mock FixedTarget.from_name class method for tests without remote data

    Actually called in conftest.py
    """
    from .target import FixedTarget

    if not hasattr(FixedTarget, '_real_from_name'):
        FixedTarget._real_from_name = FixedTarget.from_name
        FixedTarget.from_name = FixedTarget._from_name_mock
    #otherwise already mocked

def _unmock_remote_data():
    """
    undo _mock_remote_data
    """
    from .target import FixedTarget

    if hasattr(FixedTarget, '_real_from_name'):
        FixedTarget.from_name = FixedTarget._real_from_name
        del FixedTarget._real_from_name
    #otherwise assume it's already correct
