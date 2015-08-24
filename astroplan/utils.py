# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from .exceptions import OldEarthRotationDataWarning
from astropy.time import Time
import astropy.units as u
from astropy.utils.data import get_pkg_data_filename, download_file
from astropy.utils import iers
import os

try:
    # Python 3
    from urllib.error import URLError
except ImportError:
    # Python 2
    from urllib2 import URLError as URLError

import warnings

__all__ = ["get_recent_IERS_A_table"]

def get_recent_IERS_A_table(warn_update=7*u.day):
    """
    Update the IERS Bulletin A table for accurate astropy calculations near
    the present time.

    The IERS A table records recent measurements and near-future predictions of
    the Earth's orientation. A warning will be raised if the tables used are
    more than ``warn_update`` old.
    """
    local_iers_a_path = get_pkg_data_filename('data/iers_bulletin_a.txt')
    table = iers.IERS_A.open(local_iers_a_path)

    # Use polar motion flag to identify last observation before predictions
    index_of_last_observation = ''.join(table['PolPMFlag_A']).index('IP')
    time_of_last_observation = Time(table['MJD'][index_of_last_observation],
                                    format='mjd')
    time_since_last_update = Time.now() - time_of_last_observation

    # If the IERS bulletin is more than 7 days old, try to download a fresh one.
    # Note that this doesn't ensure that the file that you download will be
    # newer, it just ensures that you've tried to download the latest copy.
    if warn_update < time_since_last_update:
        try:
            path_to_tmp_file = download_file(iers.IERS_A_URL, cache=False,
                                             show_progress=False)
            os.remove(local_iers_a_path)
            os.rename(path_to_tmp_file, local_iers_a_path)
            table = iers.IERS_A.open(local_iers_a_path)

        except URLError:
            warning_msg = ("Your local copy of the IERS Bulletin A table is "
                           "{:.1f} days old, and the internet can not be "
                           "accessed to download the latest copy. Earth "
                           "rotation calculations will not be as precise as "
                           "possible until you reconnect to the internet and "
                           "attempt to access the IERS Bulletin A table again. "
                           "").format(time_since_last_update.to(u.day).value)
            warnings.warn(warning_msg, OldEarthRotationDataWarning)
    return table
