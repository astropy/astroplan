# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Standard library
import warnings

# Third-party
import numpy as np
from astropy.utils.iers import IERS_Auto
from astropy.time import Time
import astropy.units as u

# Package
from .exceptions import OldEarthOrientationDataWarning

__all__ = ["download_IERS_A",
           "time_grid_from_range", "_set_mpl_style_sheet",
           "stride_array"]

IERS_A_WARNING = ("For best precision (on the order of arcseconds), you must "
                  "download an up-to-date IERS Bulletin A table. To do so, run:"
                  "\n\n"
                  ">>> from astroplan import download_IERS_A\n"
                  ">>> download_IERS_A()\n")

# IF IERS table is unavailable we override the time deltas but need a way to
# restore them next time table is available.
BACKUP_Time_get_delta_ut1_utc = Time._get_delta_ut1_utc


def _low_precision_utc_to_ut1(self, jd1, jd2):
    """
    When no IERS Bulletin A is available (no internet connection), use low
    precision time conversion by assuming UT1-UTC=0 always.
    This method mimics `~astropy.coordinates.builtin_frames.utils.get_dut1utc`
    """
    try:
        if self.mjd*u.day not in IERS_Auto.open()['MJD']:
            warnings.warn(IERS_A_WARNING, OldEarthOrientationDataWarning)
        return self.delta_ut1_utc

    except (AttributeError, ValueError):
        warnings.warn(IERS_A_WARNING, OldEarthOrientationDataWarning)
        return np.zeros(self.shape)


def download_IERS_A(show_progress=True):
    """
    Download and cache the IERS Bulletin A table.

    If one is already cached, download a new one and overwrite the old. Store
    table in the astropy cache, and undo the monkey patching caused by earlier
    failure (if applicable).

    If one does not exist, monkey patch `~astropy.time.Time._get_delta_ut1_utc`
    so that `~astropy.time.Time` objects don't raise errors by computing UT1-UTC
    off the end of the IERS table.

    Parameters
    ----------
    show_progress : bool
        `True` shows a progress bar during the download.
    """
    # Let astropy handle all the details.
    try:
        IERS_Auto()
        # Undo monkey patch set up by exception below.
        if Time._get_delta_ut1_utc != BACKUP_Time_get_delta_ut1_utc:
            Time._get_delta_ut1_utc = BACKUP_Time_get_delta_ut1_utc
        return
    except Exception:
        warnings.warn(IERS_A_WARNING, OldEarthOrientationDataWarning)
        Time._get_delta_ut1_utc = _low_precision_utc_to_ut1


@u.quantity_input(time_resolution=u.hour)
def time_grid_from_range(time_range, time_resolution=0.5*u.hour):
    """
    Get linearly-spaced sequence of times.

    Parameters
    ----------
    time_range : `~astropy.time.Time` (length = 2)
        Lower and upper bounds on time sequence.

    time_resolution : `~astropy.units.Quantity` (optional)
        Time-grid spacing

    Returns
    -------
    times : `~astropy.time.Time`
        Linearly-spaced sequence of times
    """
    try:
        start_time, end_time = time_range
    except ValueError:
        raise ValueError("time_range should have a length of 2: lower and "
                         "upper bounds on the time sequence.")
    return Time(np.arange(start_time.jd, end_time.jd,
                          time_resolution.to(u.day).value), format='jd')


def _set_mpl_style_sheet(style_sheet):
    """
    Import matplotlib, set the style sheet to ``style_sheet`` using
    the most backward compatible import pattern.
    """
    import matplotlib
    matplotlib.rcdefaults()
    matplotlib.rcParams.update(style_sheet)


def stride_array(arr, window_width):
    """
    Computes all possible sequential subarrays of arr with length = window_width

    Parameters
    ----------
    arr : array-like (length = n)
        Linearly-spaced sequence

    window_width : int
        Number of elements in each new sub-array

    Returns
    -------
    strided_arr : array (shape = (n-window_width, window_width))
        Linearly-spaced sequence of times
    """
    as_strided = np.lib.stride_tricks.as_strided

    new_shape = (len(arr) - window_width + 1, window_width)

    strided_arr = as_strided(arr, new_shape, (arr.strides[0], arr.strides[0]))

    return strided_arr
