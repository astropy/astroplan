# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import astropy.units as u
import numpy as np

__all__ = ['Occultation']


class PeriodicEvent(object):
    """
    A periodic event defined by an epoch and period.
    """
    @u.quantity_input(period=u.day, duration=u.day)
    def __init__(self, epoch, period, duration=None, name=None):
        self.epoch = epoch
        self.period = period
        self.name = name
        self.duration = duration

    def phase(self, time):
        return ((time - self.epoch).to(u.day).value %
                 self.period.to(u.day).value) / self.period.to(u.day).value


class EclipsingSystem(PeriodicEvent):
    """
    Define some parameters for an eclipsing system, either an eclipsing binary
    or a transiting exoplanet.
    """
    def primary_eclipse(self, time):
        phases = self.phase(time)
        return ((phases < float(self.duration/self.period)/2) |
                (phases > 1 - float(self.duration/self.period)/2))

    def secondary_eclipse(self, time, secondary_eclipse_phase=0.5):
        phases = self.phase(time)
        return ((phases < secondary_eclipse_phase + float(self.duration/self.period)/2) |
                (phases > secondary_eclipse_phase - float(self.duration/self.period)/2))

    def out_of_eclipse(self, time):
        return np.logical_not(self.primary_eclipse(time) | self.secondary_eclipse(time))

    def next_primary_eclipse(self, time, n_eclipses=1):
        eclipse_times = ((1-self.phase(time)) * self.period + time +
                         np.arange(n_eclipses) * self.period)
        return eclipse_times if len(eclipse_times) > 1 else eclipse_times[0]

    def next_secondary_eclipse(self, time, n_eclipses=1):
        phase = self.phase(time)
        if phase >= 0.5:
            next_eclipse_phase = 1.5
        else:
            next_eclipse_phase = 0.5
        eclipse_times = ((next_eclipse_phase - phase) * self.period + time +
                         np.arange(n_eclipses) * self.period)
        return eclipse_times if len(eclipse_times) > 1 else eclipse_times[0]
