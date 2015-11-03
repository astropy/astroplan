# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tools for scheduling observations.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import copy
from abc import ABCMeta, abstractmethod

import numpy as np

from astropy import units as u

__all__ = ['ObservingBlock', 'TransitionBlock',
           'Scheduler', 'SequentialScheduler']

class ObservingBlock(object):
    @u.quantity_input(duration=u.second)
    def __init__(self, target, duration, configuration={}, constraints=None):
        self.target = target
        self.duration = duration
        self.configuration = configuration
        self.constraints = constraints
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

class TransitionBlock(object):
    """
    An object that represents "dead time" between observations, while the
    telescope is slewing, instrument is reconfiguring, etc.

    Parameters
    ----------
    components : dict
        A dictionary mapping the reason for an observation's dead time to
        `Quantity`s with time units
    start_time : Quantity with time units

    """
    def __init__(self, components, start_time=None):
        self.start_time = start_time
        self.components = components

    def __repr__(self):
        orig_repr = object.__repr__(self)
        if self.start_time is None or self.end_time is None:
            return orig_repr.replace('object at', 'object ({0}, unscheduled) at'.format(self.components))
        else:
            s = 'object ({0}, {1} to {2}) at'.format(self.components, self.start_time, self.end_time)
            return orig_repr.replace('object at', s)

    @property
    def end_time(self):
        return self.start_time + self.duration

    @property
    def components(self):
        return self._components
    @components.setter
    def components(self, val):
        duration = 0*u.second
        for t in val.values():
            duration += t

        self._components = val
        self.duration = duration


class Scheduler(object):
    __metaclass__ = ABCMeta

    def __call__(self, blocks):
        """
        Schedule a set of `ObservingBlock`s

        Parameters
        ----------
        blocks : iterable of `ObservingBlock`s
            The blocks to schedule.  Note that these blocks will *not*
            be modified - new ones will be created and returned.

        Returns
        -------
        schedule : list
            A list of `ObservingBlock`s and `TransitionBlock`s with populated
            `start_time` and `end_time` attributes
        """
        #these are *shallow* copies
        copied_blocks = [copy.copy(block) for block in blocks]
        new_blocks, already_sorted = self._make_schedule(copied_blocks)
        if not already_sorted:
            block_time_map = {block.start_time : block for block in new_blocks}
            new_blocks = [block_time_map[time] for time in sorted(block_time_map)]
        return new_blocks

    @abstractmethod
    def _make_schedule(self, blocks):
        """
        Does the actual business of scheduling. The `blocks` passed in should
        have their `start_time` and `end_time` modified to reflect the schedule.
        any necessary `TransitionBlock` should also be added.  Then the full set
        of blocks should be returned as a list of blocks, along with a boolean
        indicating whether or not they have been put in order already.

        Parameters
        ----------
        blocks : list of `ObservingBlock`s
            Can be modified as it is already copied by `__call__`

        Returns
        -------
        new_blocks : list of blocks
            The blocks from ``blocks``, as well as any necessary
            `TransitionBlock`s
        already_sorted : bool
            If True, the ``new_blocks`` come out pre-sorted, otherwise they need
            to be sorted.
        """
        raise NotImplementedError
        return new_blocks, already_sorted

class SequentialScheduler(Scheduler):
    """
    A scheduler that does "stupid simple sequential scheduling".  That is, it
    simply looks at all the blocks, picks the best one, schedules it, and then
    moves on.

    Parameters
    ----------
    start_time : `~astropy.time.Time`
        the start of the observation scheduling window.
    end_time : `~astropy.time.Time`
        the end of the observation scheduling window.
    constraints : sequence of `Constraint`s
        The constraints to apply to *every* observing block.  Note that
        constraints for specific blocks can go on each block individually.
    observer : `astroplan.Observer`
        The observer/site to do the scheduling for.
    configuration_transitions : dict
        TBD
    gap_time : `Quantity` with time units
        The minimal spacing to try over a gap where nothing can be scheduled.

    """
    @u.quantity_input(gap_time=u.second)
    def __init__(self, start_time, end_time, constraints, observer,
                       configuration_transitions={}, gap_time=30*u.min):
        self.constraints = constraints
        self.start_time = start_time
        self.end_time = end_time
        self.observer = observer
        self.configuration_transitions = configuration_transitions
        self.gap_time = gap_time

    @classmethod
    @u.quantity_input(duration=u.second)
    def from_timespan(cls, center_time, duration, **kwargs):
        """
        Create a new instance of this class given a time and
        """
        start_time = center_time - duration/2.
        end_time = center_time + duration/2.
        return cls(start_time, end_time, **kwargs)

    def _make_schedule(self, blocks):
        for b in blocks:
            if b.constraints is None:
                b._all_constraints = self.constraints
            else:
                b._all_constraints = self.constraints + b.constraints
            b._duration_offsets = u.Quantity([0*u.second, b.duration/2, b.duration])


        new_blocks = []
        current_time = self.start_time
        while (len(blocks) > 0) and (current_time < self.end_time):

            # first compute the value of all the constraints for each block
            # given the current starting time
            block_constraint_results = []
            for b in blocks:
                times = current_time + b._duration_offsets

                constraint_res = []
                for constraint in b._all_constraints:
                    constraint_res.append(constraint(self.observer, [b.target], times))
                # take the product over all the constraints *and* times
                block_constraint_results.append(np.prod(constraint_res))

            # now identify the block that's the best
            bestblock_idx = np.argmax(block_constraint_results)

            if block_constraint_results[bestblock_idx] == 0.:
                # if even the best is unobservable, we need a gap
                new_blocks.append(TransitionBlock({'nothing_observable': self.gap_time}, current_time))
                current_time += self.gap_time
            else:
                # If there's a best one that's observable, assign it times
                newb = blocks.pop(bestblock_idx)
                newb.start_time = current_time
                current_time += self.gap_time
                newb.end_time = current_time
                newb.constraints_value = block_constraint_results[bestblock_idx]

                new_blocks.append(newb)

        return new_blocks, True







