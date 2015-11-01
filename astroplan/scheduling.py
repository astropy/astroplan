# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tools for scheduling observations.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from abc import ABCMeta, abstractmethod

from astropy import units as u

__all__ = ['ObservingBlock']

class ObservingBlock(object):
    def __init__(self, target, duration, configuration={}):
        self.target = target
        self.duration = duration
        self.configuration = configuration
        self.start_time = self.end_time = None

    def __repr__(self):
        orig_repr = object.__repr__(self)
        if self.start_time is None or self.end_time is None:
            return orig_repr.replace('object at', 'object (unscheduled) at')
        else:
            s = 'object ({0} to {1}) at'.format(self.start_time, self.end_time)
            return orig_repr.replace('object at', s)

    @classmethod
    def from_exposures(cls, target, timeperexp, nexp, readouttime=0*u.second,
                            configuration={}):
        duration = nexp*(timeperexp + readouttime)
        ob = cls(target, duration, configuration)
        ob.timeperexp = timeperexp
        ob.nexp = nexp
        ob.readouttime = readouttime
        return ob

class TransitionBlock(ObservingBlock):
    """
    An object that represents "dead time" between observations, while the
    telescope is slewing, instrument is reconfiguring, etc.

    Parameters
    ----------
    components : dict 
        A dictionary mapping the reason for an observation's dead tune to 
        `Quantity`s with time units
    start_time : Quantity with time units

    """
    @u.quantity_input(start_time=u.second)
    def __init__(self, components, start_time=None):
        self.start_time = start_time
        self.components = components

    @property
    def end_time(self):
        return self.start_time + self.duration

    @property
    def components(self):
        return self._components
    @components.setter
    def components(self, val):
        duration = 0*u.second
        for t in components.values():
            duration += t

        self._components = val
        self.duration = duration
    
       
class Scheduler(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def __call__(self, target):
        """
        Does the actual business of scheduling.

        Parameters
        ----------
        blocks : sequence of `ObservingBlock`s
            The blocks to schedule

        Returns
        -------
        schedule : list
            A list of `ObservingBlock`s and `TransitionBlock`s with populated
            `start_time` and `end_time` attributes
        """

class SequentialScheduler(object):
    """
    A scheduler that does "stupid simple sequential scheduling".  That is, it
    simply looks at all the blocks, picks the best one, schedules it, and then
    moves on.
    """
    def __init__(self, constraints):
        self.constraints = constraints