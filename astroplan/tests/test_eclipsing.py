from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from astropy.time import Time
import astropy.units as u

from ..eclipsing import PeriodicEvent, EclipsingSystem


def test_phase():
    epoch = Time('2016-01-01')
    period = 3*u.day
    duration = 1*u.hour
    pe = PeriodicEvent(epoch=epoch, period=period, duration=duration,
                       name='test event')

    assert pe.phase(Time('2016-01-02 12:00')) == 0.5
    assert pe.phase(Time('2016-01-04 00:00')) == 0.0


def test_primary_secondary_eclipse():
    epoch = Time('2016-01-01')
    period = 3*u.day
    duration = 1*u.hour
    es = EclipsingSystem(epoch=epoch, period=period, duration=duration,
                       name='test event')

    assert es.next_primary_eclipse(epoch) == Time('2016-01-04 00:00')
    assert es.next_secondary_eclipse(epoch) == Time('2016-01-02 12:00')


def test_next_eclipse():
    epoch = Time('2016-01-01')
    period = 3*u.day
    duration = 1*u.hour
    es = EclipsingSystem(epoch=epoch, period=period, duration=duration,
                       name='test event')

    assert es.next_primary_eclipse(epoch) == Time('2016-01-04 00:00')
    assert es.next_secondary_eclipse(epoch) == Time('2016-01-02 12:00')


def test_out_of_eclipse():
    epoch = Time('2016-01-01')
    period = 3*u.day
    duration = 1*u.hour
    es = EclipsingSystem(epoch=epoch, period=period, duration=duration,
                         name='test event')

    test_times = Time(['2016-01-02 12:00', '2016-01-03 00:00',
                       '2016-01-04 00:00', '2016-01-05 00:00'])
    assert np.all(es.out_of_eclipse(test_times) ==
                  np.array([False, True, False, True]))


def test_next_eclipse():
    epoch = Time('2016-01-01')
    period = 3*u.day
    duration = 1*u.hour
    es = EclipsingSystem(epoch=epoch, period=period, duration=duration,
                         name='test event')

    test_time = epoch + period + 0.1*u.min
    next_primary = es.next_primary_eclipse(test_time)
    next_secondary = es.next_secondary_eclipse(test_time)

    assert next_primary == Time('2016-01-04')
    assert next_secondary == Time('2016-01-05 12:00')
