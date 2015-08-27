# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from astropy import units as u
from astropy.time import Time
from astropy.tests.helper import remote_data
from astropy.utils.data import clear_download_cache
from astropy.utils import iers

from ..utils import (download_IERS_A, IERS_A_in_cache, get_IERS_A_or_workaround,
                     BACKUP_Time_get_delta_ut1_utc)
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

    #now test that it actually works post-IERS A download:
    nowplusoneyear.ut1
