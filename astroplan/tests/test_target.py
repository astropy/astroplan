# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# Third-party
import astropy.units as u
from astropy.tests.helper import remote_data
from astropy.coordinates import SkyCoord

# Package
from ..target import FixedTarget

@remote_data
def test_FixedTarget_from_name():
    """
    Check that resolving target names with the `SkyCoord.from_name` constructor
    to produce a `FixedTarget` accurately resolves the coordinates of Polaris.
    """

    # Resolve coordinates with SkyCoord.from_name classmethod
    polaris_from_name = FixedTarget.from_name('Polaris')
    polaris_from_name = FixedTarget.from_name('Polaris', name='Target 1')
    # Coordinates grabbed from SIMBAD
    polaris_from_SIMBAD = SkyCoord('02h31m49.09456s', '+89d15m50.7923s')

    # Make sure separation is small
    assert polaris_from_name.coord.separation(polaris_from_SIMBAD) < 1*u.arcsec

def test_FixedTarget_ra_dec():
    """
    Confirm that FixedTarget.ra and FixedTarget.dec are the same as the
    right ascension and declination stored in the FixedTarget.coord variable -
    which is a SkyCoord
    """

    vega_coords = SkyCoord('18h36m56.33635s', '+38d47m01.2802s')
    vega = FixedTarget(vega_coords, name='Vega')
    assert vega.coord == vega_coords, 'Store coordinates directly'
    assert vega.coord.ra == vega_coords.ra == vega.ra, ('Retrieve RA from '
                                                        'SkyCoord')
    assert vega.coord.dec == vega_coords.dec == vega.dec, ('Retrieve Dec from '
                                                           'SkyCoord')
