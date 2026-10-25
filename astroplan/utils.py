# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Third-party
import numpy as np
from astropy.time import Time
import astropy.units as u
from astropy.utils.decorators import deprecated

__all__ = ["download_IERS_A",
           "time_grid_from_range", "_set_mpl_style_sheet",
           "stride_array"]


@deprecated("0.11", message=(
    "download_IERS_A is deprecated and does nothing. astropy downloads and "
    "updates the IERS tables it needs itself; to update the tables bundled "
    "with astropy, upgrade the astropy-iers-data package. See the astroplan "
    "documentation page \"What are the IERS tables and how do I update "
    "them?\" for details."))
def download_IERS_A(show_progress=True):
    """
    Formerly downloaded and cached the IERS Bulletin A table.

    astropy now handles the IERS tables itself, so this function does
    nothing. See :ref:`iers` for how astropy gets and updates them.

    Parameters
    ----------
    show_progress : bool
        Ignored.
    """


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


def _make_cache_key(times, targets):
    """
    Make a unique key to reference this combination of ``times`` and ``targets``.

    Often, we wish to store expensive calculations for a combination of
    ``targets`` and ``times`` in a cache on an ``observer`` object. This
    routine will provide an appropriate, hashable, key to store these
    calculations in a dictionary.

    Parameters
    ----------
    times : `~astropy.time.Time`
        Array of times on which to test the constraint.
    targets : `~astropy.coordinates.SkyCoord`
        Target or list of targets.

    Returns
    -------
    cache_key : tuple
        A hashable tuple for use as a cache key
    """
    # make a tuple from times
    try:
        timekey = tuple(times.ravel().jd) + times.shape + (times.scale,)
    except BaseException:        # must be scalar
        timekey = (times.jd, times.scale)
    # make hashable thing from targets coords
    try:
        if hasattr(targets, 'frame'):
            # treat as a SkyCoord object. Accessing the longitude
            # attribute of the frame data should be unique and is
            # quicker than accessing the ra attribute.
            targkey = tuple(targets.frame.data.lon.value.ravel()) + targets.shape
        else:
            # assume targets is a string.
            targkey = (targets,)
    except BaseException:
        targkey = (targets.frame.data.lon,)
    return timekey + targkey
