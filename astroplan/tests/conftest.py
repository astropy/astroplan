# Licensed under a 3-clause BSD style license - see LICENSE.rst

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import SkyCoord

from astroplan.target import Target


class ObserverDependentTarget(Target):
    """Minimal target that requires observer context when evaluated."""

    def __init__(self, name="observer-dependent"):
        self.name = name

    @property
    def is_time_dependent(self):
        return True

    def get_skycoord(self, times, observer=None):
        if observer is None:
            raise ValueError("`observer` is required.")

        ra = np.broadcast_to(
            observer.location.lon.to_value(u.deg),
            times.shape
        ) * u.deg
        dec = np.broadcast_to(
            observer.location.lat.to_value(u.deg),
            times.shape,
        ) * u.deg
        return SkyCoord(ra=ra, dec=dec)


@pytest.fixture
def observer_dependent_target():
    return ObserverDependentTarget()
