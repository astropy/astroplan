# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

__all__ = ["_mock_remote_data"]

def _mock_remote_data():
    # Mock FixedTarget.from_name class method for tests without remote data
    from astroplan.core import FixedTarget
    FixedTarget.from_name = FixedTarget._fixed_target_from_name_mock
