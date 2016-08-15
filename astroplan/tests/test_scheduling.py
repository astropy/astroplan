# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import SkyCoord

from ..utils import time_grid_from_range
from ..observer import Observer
from ..target import FixedTarget
from ..constraints import (AirmassConstraint, _get_altaz)
from ..scheduling import (ObservingBlock, PriorityScheduler, SequentialScheduler,
                          Transitioner, TransitionBlock, Schedule, Slot, Scorer)

vega = FixedTarget(coord=SkyCoord(ra=279.23473479 * u.deg, dec=38.78368896 * u.deg),
                   name="Vega")
rigel = FixedTarget(coord=SkyCoord(ra=78.63446707 * u.deg, dec=8.20163837 * u.deg),
                    name="Rigel")
polaris = FixedTarget(coord=SkyCoord(ra=37.95456067 * u.deg,
                                     dec=89.26410897 * u.deg), name="Polaris")

apo = Observer.at_site('apo')
targets = [vega, rigel, polaris]


def test_transitioner():
    blocks = [ObservingBlock(t, 55 * u.minute, i) for i, t in enumerate(targets)]
    slew_rate = 1 * u.deg / u.second
    trans = Transitioner(slew_rate=slew_rate)
    start_time = Time('2016-02-06 00:00:00')
    transition = trans(blocks[0], blocks[1], start_time, apo)
    aaz = _get_altaz(Time([start_time]), apo,
                     [blocks[0].target, blocks[1].target])['altaz']
    sep = aaz[0].separation(aaz[1])[0]
    assert isinstance(transition, TransitionBlock)
    assert transition.duration == sep/slew_rate
transitioner = Transitioner(slew_rate=1 * u.deg / u.second)


def test_priority_scheduler():
    constraints = [AirmassConstraint(3, boolean_constraint=False)]
    blocks = [ObservingBlock(t, 55*u.minute, i) for i, t in enumerate(targets)]
    start_time = Time('2016-02-06 00:00:00')
    end_time = start_time + 24*u.hour
    scheduler = PriorityScheduler(transitioner=transitioner,
                                  constraints=constraints, observer=apo)
    schedule = Schedule(start_time, end_time)
    scheduler(blocks, schedule)
    assert len(schedule.observing_blocks) == 3
    assert all(np.abs(block.end_time - block.start_time - block.duration) <
               1*u.second for block in schedule.scheduled_blocks)


def test_sequential_scheduler():
    constraints = [AirmassConstraint(2.5, boolean_constraint=False)]
    blocks = [ObservingBlock(t, 55 * u.minute, i) for i, t in enumerate(targets)]
    start_time = Time('2016-02-06 00:00:00')
    end_time = start_time + 24 * u.hour
    scheduler = SequentialScheduler(constraints=constraints, observer=apo,
                                    transitioner=transitioner)
    schedule = Schedule(start_time, end_time)
    scheduler(blocks, schedule)
    assert len(schedule.observing_blocks) > 0
    assert all(np.abs(block.end_time - block.start_time - block.duration) <
               1*u.second for block in schedule.scheduled_blocks)


def test_slot():
    start_time = Time('2016-02-06 00:00:00')
    end_time = start_time + 24 * u.hour
    slot = Slot(start_time, end_time)
    slots = slot.split_slot(start_time, start_time+1*u.hour)
    assert len(slots) == 2
    assert slots[0].end == slots[1].start


def test_schedule():
    start_time = Time('2016-02-06 00:00:00')
    end_time = start_time + 24 * u.hour
    schedule = Schedule(start_time, end_time)
    assert schedule.slots[0].start == start_time
    assert schedule.slots[0].end == end_time
    assert schedule.slots[0].duration == 24*u.hour
    schedule.new_slots(0, start_time, end_time)
    assert len(schedule.slots) == 1
    new_slots = schedule.new_slots(0, start_time+1*u.hour, start_time+4*u.hour)
    assert np.abs(new_slots[0].duration - 1*u.hour) < 1*u.second
    assert np.abs(new_slots[1].duration - 3*u.hour) < 1*u.second
    assert np.abs(new_slots[2].duration - 20*u.hour) < 1*u.second


def test_scorer():
    constraint = AirmassConstraint(max=4)
    times = time_grid_from_range(Time(['2016-02-06 00:00', '2016-02-06 08:00']),
                                 time_resolution=20*u.minute)
    c = constraint(apo, [vega, rigel], times)
    block = ObservingBlock(vega, 1*u.hour, 0, constraints=[constraint])
    block2 = ObservingBlock(rigel, 1*u.hour, 0, constraints=[constraint])
    scorer = Scorer.from_start_end([block, block2], apo, Time('2016-02-06 00:00'),
                                   Time('2016-02-06 08:00'))
    scores = scorer.create_score_array(time_resolution=20*u.minute)
    assert np.array_equal(c, scores)

    constraint2 = AirmassConstraint(max=2, boolean_constraint=False)
    c2 = constraint2(apo, [vega, rigel], times)
    block = ObservingBlock(vega, 1*u.hour, 0, constraints=[constraint])
    block2 = ObservingBlock(rigel, 1*u.hour, 0, constraints=[constraint2])
    # vega's score should be = c[0], rigel's should be =  c2[1]
    scorer = Scorer.from_start_end([block, block2], apo, Time('2016-02-06 00:00'),
                                   Time('2016-02-06 08:00'))
    scores = scorer.create_score_array(time_resolution=20 * u.minute)
    assert np.array_equal(c[0], scores[0])
    assert np.array_equal(c2[1], scores[1])

    block = ObservingBlock(vega, 1*u.hour, 0)
    block2 = ObservingBlock(rigel, 1*u.hour, 0)
    scorer = Scorer.from_start_end([block, block2], apo, Time('2016-02-06 00:00'),
                                   Time('2016-02-06 08:00'), [constraint2])
    scores = scorer.create_score_array(time_resolution=20 * u.minute)
    # the ``global_constraint``: constraint2 should have applied to the blocks
    assert np.array_equal(c2, scores)
