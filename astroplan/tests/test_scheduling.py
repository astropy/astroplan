# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import SkyCoord

from ..observer import Observer
from ..target import FixedTarget
from ..constraints import (AirmassConstraint, AtNightConstraint, _get_altaz)
from ..scheduling import (ObservingBlock, PriorityScheduler, SequentialScheduler,
                          Transitioner, TransitionBlock, Schedule, Slot)

vega = FixedTarget(coord=SkyCoord(ra=279.23473479 * u.deg, dec=38.78368896 * u.deg),
                   name="Vega")
rigel = FixedTarget(coord=SkyCoord(ra=78.63446707 * u.deg, dec=8.20163837 * u.deg),
                    name="Rigel")
polaris = FixedTarget(coord=SkyCoord(ra=37.95456067 * u.deg,
                                     dec=89.26410897 * u.deg), name="Polaris")

apo = Observer.at_site('apo')
targets = [vega, rigel, polaris]
default_time = Time('2016-02-06 00:00:00')
only_at_night = [AtNightConstraint()]


def test_observing_block():
    block = ObservingBlock(rigel, 1*u.minute, 0, configuration={'filter': 'b'})
    assert(block.configuration['filter'] == 'b')
    assert(block.target == rigel)
    times_per_exposure = [1*u.minute, 4*u.minute, 15*u.minute, 5*u.minute]
    numbers_of_exposures = [100, 4, 3, 12]
    readout_time = 0.5*u.minute
    for index in range(len(times_per_exposure)):
        block = ObservingBlock.from_exposures(vega, 0, times_per_exposure[index],
                                              numbers_of_exposures[index], readout_time)
        assert(block.duration == numbers_of_exposures[index] *
               (times_per_exposure[index] + readout_time))


def test_transitioner():
    blocks = [ObservingBlock(t, 55 * u.minute, i) for i, t in enumerate(targets)]
    slew_rate = 1 * u.deg / u.second
    trans = Transitioner(slew_rate=slew_rate)
    start_time = default_time
    transition = trans(blocks[0], blocks[1], start_time, apo)
    aaz = _get_altaz(Time([start_time]), apo,
                     [blocks[0].target, blocks[1].target])['altaz']
    sep = aaz[0].separation(aaz[1])[0]
    assert isinstance(transition, TransitionBlock)
    assert transition.duration == sep/slew_rate
    blocks = [ObservingBlock(vega, 10*u.minute, 0, configuration={'filter': 'v'}),
              ObservingBlock(vega, 10*u.minute, 0, configuration={'filter': 'i'}),
              ObservingBlock(rigel, 10*u.minute, 0, configuration={'filter': 'i'})]
    trans = Transitioner(slew_rate, instrument_reconfig_times={'filter': {('v', 'i'): 2*u.minute,
                                                                          'default': 5*u.minute}})
    transition1 = trans(blocks[0], blocks[1], start_time, apo)
    transition2 = trans(blocks[0], blocks[2], start_time, apo)
    transition3 = trans(blocks[1], blocks[0], start_time, apo)
    assert np.abs(transition1.duration - 2*u.minute) < 1*u.second
    assert np.abs(transition2.duration - 2*u.minute - transition.duration) < 1*u.second
    # to test the default transition
    assert np.abs(transition3.duration - 5*u.minute) < 1*u.second

transitioner = Transitioner(slew_rate=1 * u.deg / u.second)


def test_priority_scheduler():
    constraints = [AirmassConstraint(3, boolean_constraint=False)]
    blocks = [ObservingBlock(t, 55*u.minute, i) for i, t in enumerate(targets)]
    start_time = default_time
    end_time = start_time + 24*u.hour
    scheduler = PriorityScheduler(start_time, end_time, transitioner=transitioner,
                                  constraints=constraints, observer=apo)
    schedule = scheduler(blocks)
    assert len(schedule.observing_blocks) == 3
    assert all(np.abs(block.end_time - block.start_time - block.duration) <
               1*u.second for block in schedule.scheduled_blocks)
    assert None is None


def test_sequential_scheduler():
    constraints = [AirmassConstraint(2.5, boolean_constraint=False)]
    blocks = [ObservingBlock(t, 55 * u.minute, i) for i, t in enumerate(targets)]
    start_time = default_time
    end_time = start_time + 24 * u.hour
    scheduler = SequentialScheduler(start_time, end_time,
                                    constraints=constraints, observer=apo,
                                    transitioner=transitioner)
    schedule = scheduler(blocks)
    assert len(schedule.observing_blocks) > 0
    assert all(np.abs(block.end_time - block.start_time - block.duration) <
               1*u.second for block in schedule.scheduled_blocks)


def test_scheduling_target_down():
    lco = Observer.at_site('lco')
    block = [ObservingBlock(FixedTarget.from_name('polaris'), 1 * u.min, 0)]
    start_time = default_time
    end_time = start_time + 5*u.day
    scheduler1 = SequentialScheduler(start_time, end_time, only_at_night, lco,
                                     transitioner)
    schedule1 = scheduler1(block)
    assert len(schedule1.observing_blocks) == 0
    scheduler2 = PriorityScheduler(start_time, end_time, only_at_night, lco,
                                   transitioner)
    schedule2 = scheduler2(block)
    assert len(schedule2.observing_blocks) == 0


def test_scheduling_during_day():
    block = [ObservingBlock(FixedTarget.from_name('polaris'), 1 * u.min, 0)]
    day = default_time
    start_time = apo.midnight(day) + 10*u.hour
    end_time = start_time + 6*u.hour
    scheduler1 = SequentialScheduler(start_time, end_time, only_at_night, apo,
                                     transitioner)
    schedule1 = scheduler1(block)
    assert len(schedule1.observing_blocks) == 0
    scheduler2 = PriorityScheduler(start_time, end_time, only_at_night, apo,
                                   transitioner)
    schedule2 = scheduler2(block)
    assert len(schedule2.observing_blocks) == 0


def test_slot():
    start_time = default_time
    end_time = start_time + 24 * u.hour
    slot = Slot(start_time, end_time)
    slots = slot.split_slot(start_time, start_time+1*u.hour)
    assert len(slots) == 2
    assert slots[0].end == slots[1].start


def test_schedule():
    start_time = default_time
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

