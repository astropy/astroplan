# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from astropy.time import Time
import astropy.units as u
from astropy.coordinates import SkyCoord

from ..observer import Observer
from ..target import FixedTarget
from ..constraints import (AirmassConstraint)
from ..scheduling import (ObservingBlock, PriorityScheduler, SequentialScheduler, Transitioner)

vega = FixedTarget(coord=SkyCoord(ra=279.23473479 * u.deg, dec=38.78368896 * u.deg),
                   name="Vega")
rigel = FixedTarget(coord=SkyCoord(ra=78.63446707 * u.deg, dec=8.20163837 * u.deg),
                    name="Rigel")
polaris = FixedTarget(coord=SkyCoord(ra=37.95456067 * u.deg,
                                     dec=89.26410897 * u.deg), name="Polaris")

mdm = Observer.at_site('mdm')
targets = [vega, rigel, polaris]


def test_priority_scheduler():
    constraints = [AirmassConstraint(3, boolean_constraint=False)]
    blocks = [ObservingBlock(t, 55*u.minute, i) for i, t in enumerate(targets)]
    start_time = Time('2016-02-06 00:00:00')
    end_time = start_time + 24*u.hour
    scheduler = PriorityScheduler(start_time, end_time,
                                  constraints=constraints, observer=mdm)
    schedule = scheduler(blocks)
    assert len(schedule.observing_blocks) > 0


def test_sequential_scheduler():
    constraints = [AirmassConstraint(2.5, boolean_constraint=False)]
    blocks = [ObservingBlock(t, 55 * u.minute, i) for i, t in enumerate(targets)]
    start_time = Time('2016-02-06 00:00:00')
    end_time = start_time + 24 * u.hour
    trans = Transitioner(slew_rate=1*u.deg/u.second)
    scheduler = SequentialScheduler(start_time, end_time,
                                    constraints=constraints, observer=mdm,
                                    transitioner=trans)
    schedule = scheduler(blocks)
    assert len(schedule.observing_blocks) > 0
