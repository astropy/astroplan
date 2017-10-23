# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
import pytest

from astropy import units as u
from astropy.time import Time
from astropy.tests.helper import remote_data
from astropy.utils.data import clear_download_cache
from astropy.utils import iers

from ..utils import (download_IERS_A, IERS_A_in_cache,
                     get_IERS_A_or_workaround, BACKUP_Time_get_delta_ut1_utc,
                     stride_array, time_grid_from_range)

from ..exceptions import OldEarthOrientationDataWarning


@remote_data
def test_iers_download(monkeypatch, recwarn):
    # the monkeypatch is here to undo the changes that importing astroplan does
    # if the IERS A tables already exist
    if IERS_A_in_cache():
        clear_download_cache(iers.IERS_A_URL)

    monkeypatch.setattr(iers.IERS, 'iers_table', None)
    monkeypatch.setattr(iers.IERS_A, 'iers_table', None)
    monkeypatch.setattr(Time, '_get_delta_ut1_utc', BACKUP_Time_get_delta_ut1_utc)

    # now make sure the state is what it should be given the above changes
    get_IERS_A_or_workaround()

    # first make sure a future time gives a warning with IERS A missing
    nowplusoneyear = Time.now() + 1*u.year
    nowplusoneyear.ut1
    recwarn.pop(OldEarthOrientationDataWarning)

    download_IERS_A()

    # now test that it actually works post-IERS A download:
    nowplusoneyear.ut1


arr10 = np.arange(10)


def test_stride_array():
    stride10by5 = stride_array(arr10, 5)
    assert len(stride10by5) == 6


def test_stride_floats():
    arr_float = np.asarray(arr10, float)
    stride10by3 = stride_array(arr_float, 3)


def test_time_grid_from_range():
    times2 = ['2010-01-01 00:00:00', '2010-01-01 01:00:10']
    times3 = ['2010-01-01 00:00:00', '2010-01-01 01:00:00',
              '2010-01-01 03:00:00']
    tgrid = time_grid_from_range(Time(times2))
    assert len(tgrid) == 3
    assert np.all(tgrid.iso == Time(['2010-01-01 00:00:00',
                                     '2010-01-01 00:30:00',
                                     '2010-01-01 01:00:00']))

    with pytest.raises(ValueError):
        time_grid_from_range(Time(times3))
