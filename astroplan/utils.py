# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# Standard library
import warnings

# Third-party
import numpy as np
from astropy.utils.iers import IERS_Auto
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import EarthLocation

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

    time_resolution : `~astropy.units.quantity` (optional)
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


def _mock_remote_data():
    """
    Apply mocks (i.e. monkey-patches) to avoid the need for internet access
    for certain things.

    This is currently called in `astroplan/conftest.py` when the tests are run
    and the `--remote-data` option isn't used.

    The way this setup works is that for functionality that usually requires
    internet access, but has mocks in place, it is possible to write the test
    without adding a `@remote_data` decorator, and `py.test` will do the right
    thing when running the tests:

    1. Access the internet and use the normal code if `--remote-data` is used
    2. Not access the internet and use the mock code if `--remote-data` is not used

    Both of these cases are tested on travis-ci.
    """
    from .target import FixedTarget
    from astropy.coordinates import EarthLocation

    if not hasattr(FixedTarget, '_real_from_name'):
        FixedTarget._real_from_name = FixedTarget.from_name
        FixedTarget.from_name = FixedTarget._from_name_mock

    if not hasattr(EarthLocation, '_real_of_site'):
        EarthLocation._real_of_site = EarthLocation.of_site
        EarthLocation.of_site = EarthLocation_mock.of_site_mock

    # otherwise already mocked


def _unmock_remote_data():
    """
    undo _mock_remote_data
    currently unused
    """
    from .target import FixedTarget

    if hasattr(FixedTarget, '_real_from_name'):
        FixedTarget.from_name = FixedTarget._real_from_name
        del FixedTarget._real_from_name

    if hasattr(EarthLocation, '_real_of_site'):
        EarthLocation.of_site = EarthLocation._real_of_site
        del EarthLocation._real_of_site
    # otherwise assume it's already correct


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


class EarthLocation_mock(EarthLocation):
    """
    Mock the EarthLocation class if no remote data for locations commonly
    used in the tests.
    """

    @classmethod
    def of_site_mock(cls, string):
        subaru = EarthLocation.from_geodetic(-155.4761111111111*u.deg,
                                             19.825555555555564*u.deg,
                                             4139*u.m)

        lco = EarthLocation.from_geodetic(-70.70166666666665*u.deg,
                                          -29.003333333333327*u.deg,
                                          2282*u.m)

        aao = EarthLocation.from_geodetic(149.06608611111113*u.deg,
                                          -31.277038888888896*u.deg,
                                          1164*u.m)

        vbo = EarthLocation.from_geodetic(78.8266*u.deg,
                                          12.576659999999999*u.deg,
                                          725*u.m)

        apo = EarthLocation.from_geodetic(-105.82*u.deg,
                                          32.78*u.deg,
                                          2798*u.m)

        keck = EarthLocation.from_geodetic(-155.47833333333332*u.deg,
                                           19.828333333333326*u.deg,
                                           4160*u.m)

        kpno = EarthLocation.from_geodetic(-111.6*u.deg,
                                           31.963333333333342*u.deg,
                                           2120*u.m)

        lapalma = EarthLocation.from_geodetic(-17.879999*u.deg,
                                              28.758333*u.deg,
                                              2327*u.m)

        observatories = dict(lco=lco, subaru=subaru, aao=aao, vbo=vbo, apo=apo,
                             keck=keck, kpno=kpno, lapalma=lapalma)

        return observatories[string.lower()]


def _open_shelve(shelffn, withclosing=False):
    """
    Opens a shelf file.  If ``withclosing`` is True, it will be opened with
    closing, allowing use like:

        with _open_shelve('somefile',True) as s:
            ...

    This workaround can be removed in favour of using shelve.open() directly
    once support for Python <3.4 is dropped.
    """
    import shelve
    import contextlib

    shelf = shelve.open(shelffn, protocol=2)

    if withclosing:
        return contextlib.closing(shelf)
    else:
        return shelf
