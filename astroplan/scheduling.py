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

from .utils import time_grid_from_range, stride_array

__all__ = ['ObservingBlock', 'TransitionBlock', 'Schedule', 'Slot', 'Scheduler',
           'SequentialScheduler', 'PriorityScheduler', 'Transitioner',
           'ScheduleScheduler']


class ObservingBlock(object):
    """
    An observation to be scheduled, consisting of a target and associated
    constraints on observations.
    """
    @u.quantity_input(duration=u.second)
    def __init__(self, target, duration, priority, configuration={}, constraints=None):
        """
        Parameters
        ----------
        target: `~astroplan.FixedTarget'
            Target to observe

        duration : `~astropy.units.Quantity`
            exposure time

        priority: integer or float
            priority of this object in the target list. 1 is highest priority,
            no maximum

        configuration : dict
            Configuration metadata

        constraints : list of `~astroplan.constraints.Constraint` objects
            The constraints to apply to this particular observing block.  Note
            that constraints applicable to the entire list should go into the
            scheduler.

        """
        self.target = target
        self.duration = duration
        self.priority = priority
        self.configuration = configuration
        self.constraints = constraints
        self.start_time = self.end_time = None

    def __repr__(self):
        orig_repr = object.__repr__(self)
        if self.start_time is None or self.end_time is None:
            return orig_repr.replace('object at',
                                     '({0}, unscheduled) at'
                                     .format(self.target.name))
        else:
            s = '({0}, {1} to {2}) at'.format(self.target.name, self.start_time,
                                              self.end_time)
            return orig_repr.replace('object at', s)

    @classmethod
    def from_exposures(cls, target, priority, time_per_exposure,
                       number_exposures, readout_time=0 * u.second,
                       configuration={}):
        duration = number_exposures * (time_per_exposure + readout_time)
        ob = cls(target, duration, priority, configuration)
        ob.time_per_exposure = time_per_exposure
        ob.number_exposures = number_exposures
        ob.readout_time = readout_time
        return ob


class TransitionBlock(object):
    """
    Parameterizes the "dead time", e.g. between observations, while the
    telescope is slewing, instrument is reconfiguring, etc.
    """
    def __init__(self, components, start_time=None):
        """
        Parameters
        ----------
        components : dict
            A dictionary mapping the reason for an observation's dead time to
            `~astropy.units.Quantity` objects with time units

        start_time : `~astropy.units.Quantity`
            Start time of observation
        """
        self.start_time = start_time
        self.components = components
        self._components = None
        self.duration = None

    def __repr__(self):
        orig_repr = object.__repr__(self)
        comp_info = ', '.join(['{0}: {1}'.format(c, t)
                               for c, t in self.components.items()])
        if self.start_time is None or self.end_time is None:
            return orig_repr.replace('object at', ' ({0}, unscheduled) at'.format(comp_info))
        else:
            s = '({0}, {1} to {2}) at'.format(comp_info, self.start_time, self.end_time)
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
        
    @classmethod
    @u.quantity_input(duration=u.second)
    def from_duration(cls, duration):
        # for testing how to put transitions between observations during
        # scheduling without considering the complexities of duration
        tb = TransitionBlock({None: 0*u.second})
        tb.duration = duration
        return tb


class Schedule(object):
    """
    An object that represents a astronomical schedule
    
    Parameters:
    -----------
    start_time : `~astropy.time.Time`
        The starting time of the schedule; the start of your 
        observing window
    end_time : `~astropy.time.Time`
        The ending time of the schedule; the end of your
        observing window
    constraints : sequence of `Constraint`s
        these are constraints that apply to the entire schedule
    """
    # as currently written, there should be no consecutive unoccupied slots
    # this should change to allow for more flexibility (e.g. dark slots, grey slots)
    
    def __init__(self, start_time, end_time, constraints=None):
        self.slots = [Slot(start_time, end_time)]
        self.constraints = constraints
        self.slew_duration = 4*u.min
        # TODO: replace/overwrite slew_duration with Transitioner calls

    def __repr__(self):
        # once TransitionBlocks are defined, this should be scheduled_blocks
        for block in self.observing_blocks:
            if hasattr(block, 'target'):
                try:
                    print(block.target.name, 'starting @',
                          block.start_time.iso, 'lasting', block.duration)
                except:
                    print(block.name, 'starting @', block.start_time.iso,
                          'lasting', block.duration)
            else:
                print('trans starting @', block.start_time.iso, 'lasting',
                      block.duration, ':', block.components)
        return 'done'

    def apply_constraints(self):
        # this needs to be able to handle being passed constraints
        # that are targeted and non-targeted, use the non-targeted
        # and place targeted (e.g. MoonSep) somewhere they can be used
        raise NotImplementedError
    
    @property
    def observing_blocks(self):
        return [slot.OB for slot in self.slots if slot.OB]

    @property
    def scheduled_blocks(self):
        blocks = []
        for slot in self.slots:
            if slot.OB:
                blocks.append(slot.OB)
            elif slot.TB:
                blocks.append(slot.TB)
        return blocks

    @property
    def open_slots(self):
        return [slot for slot in self.slots if not slot.occupied]

    def new_slots(self, slot_index, start_time, end_time):
        # this is intended to be used such that there aren't consecutive unoccupied slots
        new_slots = self.slots[slot_index].split_slot(start_time, end_time)
        return new_slots

    def insert_slot(self, slot_index, start_time, block):
        if block.duration + self.slew_duration > self.slots[slot_index].duration:
            raise ValueError('constraint application failed and did not leave enough time')
        if isinstance(block, ObservingBlock):
            tb = TransitionBlock.from_duration(self.slew_duration)
            if self.slots[slot_index].end-(start_time+block.duration) < tb.duration:
                self.insert_slot(slot_index, self.slots[slot_index].end-tb.duration, tb)
                start_time = self.slots[slot_index].end - block.duration
            else:
                self.insert_slot(slot_index, start_time+block.duration, tb)
            # TODO: make it shift observing/transition blocks to fill small amounts of open space
            block.end_time = start_time+block.duration
        earlier_slots = self.slots[:slot_index]
        later_slots = self.slots[slot_index+1:]
        end_time = start_time+block.duration
        block.start_time = start_time
        new_slots = self.new_slots(slot_index, start_time, end_time)
        for new_slot in new_slots:
            if new_slot.middle:
                new_slot.occupied = True
                if isinstance(block, ObservingBlock):
                    new_slot.OB = block
                elif isinstance(block, TransitionBlock):
                    new_slot.TB = block
        self.slots = earlier_slots + new_slots + later_slots
        return earlier_slots + new_slots + later_slots
    
    def remove_slot(self, slot_index, start_time, end_time):
        earlier_slots = self.slots[:slot_index]
        later_slots = self.slots[slot_index+1:]
        new_slots = self.new_slots(slot_index, start_time, end_time)
        for new_slot in new_slots:
            if new_slot.middle:
                new_slots.remove(new_slot)
        return earlier_slots + new_slots + later_slots
    
    @classmethod
    def from_constraints(cls, start_time, end_time, constraints):
        sch = cls(start_time, end_time, constraints=constraints)
        sch.apply_constraints()
        return sch


class Slot(object):
    """
    A time slot within the schedule
    
    Parameters:
    -----------
    start_time : `~astropy.time.Time`
        The starting time of the slot
    end_time : `~astropy.time.Time`
        The ending time of the slot    
    """
    
    def __init__(self, start_time, end_time):
        self.start = start_time
        self.end = end_time
        self.duration = end_time-start_time
        self.occupied = False
        self.middle = False
        self.OB = False
        self.TB = False
        
    def split_slot(self, early_time, later_time):
        # check if the new slot would overwrite occupied/other slots
        if self.occupied:
            raise ValueError('slot is already occupied')
            
        new_slot = Slot(early_time, later_time)
        new_slot.middle = True
        early_slot = Slot(self.start, early_time)
        late_slot = Slot(later_time, self.end)
        
        if early_time > self.start and later_time < self.end:
            return [early_slot, new_slot, late_slot]
        elif early_time > self.start:
            return [early_slot, new_slot]
        elif later_time < self.end:
            return [new_slot, late_slot]
        else:
            return [new_slot]

           
class Scheduler(object):
    """
    Schedule a set of `~astroplan.scheduling.ObservingBlock` objects
    """

    __metaclass__ = ABCMeta

    def __call__(self, blocks):
        """
        Parameters
        ----------
        blocks : list of `~astroplan.scheduling.ObservingBlock` objects
            The observing blocks to schedule.  Note that the input
            `~astroplan.scheduling.ObservingBlock` objects will *not* be
            modified - new ones will be created and returned.

        Returns
        -------
        schedule : list
            A list of `~astroplan.scheduling.ObservingBlock` objects and
            `~astroplan.scheduling.TransitionBlock` objects with populated
            ``start_time`` and ``end_time`` attributes
        """
        # these are *shallow* copies
        copied_blocks = [copy.copy(block) for block in blocks]
        new_blocks, already_sorted = self._make_schedule(copied_blocks)
        if not already_sorted:
            block_time_map = {block.start_time.datetime: block for block in new_blocks}
            new_blocks = [block_time_map[time] for time in sorted(block_time_map)]
        return new_blocks

    @abstractmethod
    def _make_schedule(self, blocks):
        """
        Does the actual business of scheduling. The ``blocks`` passed in should
        have their ``start_time` and `end_time`` modified to reflect the
        schedule. Any necessary `~astroplan.scheduling.TransitionBlock` should
        also be added.  Then the full set of blocks should be returned as a list
        of blocks, along with a boolean indicating whether or not they have been
        put in order already.

        Parameters
        ----------
        blocks : list of `~astroplan.scheduling.ObservingBlock` objects
            Can be modified as it is already copied by ``__call__``

        Returns
        -------
        new_blocks : list of blocks
            The blocks from ``blocks``, as well as any necessary
            `~astroplan.scheduling.TransitionBlock` objects
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
    constraints : sequence of `~astroplan.constraints.Constraint` objects
        The constraints to apply to *every* observing block.  Note that
        constraints for specific blocks can go on each block individually.
    observer : `~astroplan.Observer`
        The observer/site to do the scheduling for.
    transitioner : `~astroplan.scheduling.Transitioner` or None
        The object to use for computing transition times between blocks
    gap_time : `~astropy.units.Quantity` with time units
        The minimal spacing to try over a gap where nothing can be scheduled.

    """
    @u.quantity_input(gap_time=u.second)
    def __init__(self, start_time, end_time, constraints, observer,
                 transitioner=None, gap_time=30*u.min):
        self.constraints = constraints
        self.start_time = start_time
        self.end_time = end_time
        self.observer = observer
        self.transitioner = transitioner
        self.gap_time = gap_time

    @classmethod
    @u.quantity_input(duration=u.second)
    def from_timespan(cls, center_time, duration, **kwargs):
        """
        Create a new instance of this class given a center time and duration.

        Parameters
        ----------
        center_time : `~astropy.time.Time`
            Mid-point of time-span to schedule.

        duration : `~astropy.units.Quantity` or `~astropy.time.TimeDelta`
            Duration of time-span to schedule
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
            b._duration_offsets = u.Quantity([0*u.second, b.duration/2,
                                              b.duration])

        new_blocks = []
        current_time = self.start_time
        while (len(blocks) > 0) and (current_time < self.end_time):

            # first compute the value of all the constraints for each block
            # given the current starting time
            block_transitions = []
            block_constraint_results = []
            for b in blocks:
                # first figure out the transition
                if len(new_blocks) > 0:
                    trans = self.transitioner(new_blocks[-1], b, current_time,
                                              self.observer)
                else:
                    trans = None
                block_transitions.append(trans)
                transition_time = 0*u.second if trans is None else trans.duration

                times = current_time + transition_time + b._duration_offsets

                constraint_res = []
                for constraint in b._all_constraints:
                    constraint_res.append(constraint(self.observer, [b.target],
                                                     times))
                # take the product over all the constraints *and* times
                block_constraint_results.append(np.prod(constraint_res))

            # now identify the block that's the best
            bestblock_idx = np.argmax(block_constraint_results)

            if block_constraint_results[bestblock_idx] == 0.:
                # if even the best is unobservable, we need a gap
                new_blocks.append(TransitionBlock({'nothing_observable': self.gap_time},
                                                  current_time))
                current_time += self.gap_time
            else:
                # If there's a best one that's observable, first get its transition
                trans = block_transitions.pop(bestblock_idx)
                if trans is not None:
                    new_blocks.append(trans)
                    current_time += trans.duration

                # now assign the block itself times and add it to the schedule
                newb = blocks.pop(bestblock_idx)
                newb.start_time = current_time
                current_time += self.gap_time
                newb.end_time = current_time
                newb.constraints_value = block_constraint_results[bestblock_idx]

                new_blocks.append(newb)

        return new_blocks, True


class ScheduleScheduler(object):
    # temporary as I try to figure out how to deal with TransitionBlocks
    # for testing modifications to the PriorityScheduler

    __metaclass__ = ABCMeta

    def __call__(self, blocks):
        copied_blocks = [copy.copy(block) for block in blocks]
        schedule = self._make_schedule(copied_blocks)
        return schedule

    @abstractmethod
    def _make_schedule(self, blocks):
        """
        Does the actual business of scheduling. The ``blocks`` passed in should
        have their ``start_time` and `end_time`` modified to reflect the
        schedule. Any necessary `~astroplan.scheduling.TransitionBlock` should
        also be added.  Then the full set of blocks should be returned as a list
        of blocks, along with a boolean indicating whether or not they have been
        put in order already.

        Parameters
        ----------
        blocks : list of `~astroplan.scheduling.ObservingBlock` objects
            Can be modified as it is already copied by ``__call__``

        Returns
        -------
        new_blocks : list of blocks
            The blocks from ``blocks``, as well as any necessary
            `~astroplan.scheduling.TransitionBlock` objects
        already_sorted : bool
            If True, the ``new_blocks`` come out pre-sorted, otherwise they need
            to be sorted.
        """
        raise NotImplementedError
        return schedule


class PriorityScheduler(ScheduleScheduler):
    """
    A scheduler that optimizes a prioritized list.  That is, it
    finds the best time for each ObservingBlock, in order of priority.
    """

    @u.quantity_input(gap_time=u.second)
    def __init__(self, start_time, end_time, constraints, observer,
                 transitioner=None, gap_time=30 * u.min, slew_time=5 * u.min):
        """
        Parameters
        ----------
        start_time : `~astropy.time.Time`
            the start of the observation scheduling window.
        end_time : `~astropy.time.Time`
            the end of the observation scheduling window.
        constraints : sequence of `~astroplan.constraints.Constraint`
            The constraints to apply to *every* observing block.  Note that
            constraints for specific blocks can go on each block individually.
        observer : `~astroplan.Observer`
            The observer/site to do the scheduling for.
        transitioner : `~astroplan.scheduling.Transitioner` or None
            The object to use for computing transition times between blocks.
            Not currently used in this Scheduler.
        gap_time : `~astropy.units.Quantity` with time units
            The minimal spacing to try over a gap where nothing can be scheduled.
        slew_time : `~astropy.units.Quantity` with time units
            The time required between observations.
            Used instead of transitioner (for now)
        """
        self.constraints = constraints
        self.start_time = start_time
        self.end_time = end_time
        self.observer = observer
        self.transitioner = transitioner
        self.gap_time = gap_time
        self.slew_time = slew_time
        # make a schedule object, when apply_constraints works, add constraints
        self.schedule = Schedule(self.start_time, self.end_time,
                                 # constraints=self.constraints
                                 )
        self.schedule.slew_duration = self.slew_time

    @classmethod
    @u.quantity_input(duration=u.second)
    def from_timespan(cls, center_time, duration, **kwargs):
        """
        Create a new instance of this class given a center time and duration.

        Parameters
        ----------
        center_time : `~astropy.time.Time`
            Mid-point of time-span to schedule.

        duration : `~astropy.units.Quantity` or `~astropy.time.TimeDelta`
            Duration of time-span to schedule
        """
        start_time = center_time - duration / 2.
        end_time = center_time + duration / 2.
        return cls(start_time, end_time, **kwargs)

    def _make_schedule(self, blocks):
        # Combine individual constraints with global constraints, and
        # retrieve priorities from each block to define scheduling order
        _all_times = []
        _block_priorities = np.zeros(len(blocks))
        for i, b in enumerate(blocks):
            if b.constraints is None:
                b._all_constraints = self.constraints
            else:
                b._all_constraints = self.constraints + b.constraints
            b._duration_offsets = u.Quantity([0 * u.second, b.duration / 2, b.duration])
            _block_priorities[i] = b.priority
            _all_times.append(b.duration)

        # Define a master schedule
        # Generate grid of time slots, and a mask for previous observations

        # Find the minimum required time step
        # TODO: a common factorization of all times is probably better long-term
        time_resolution = min(min(_all_times), self.slew_time)
        times = time_grid_from_range([self.start_time, self.end_time],
                                     time_resolution=time_resolution)
        is_open_time = np.ones(len(times), bool)

        # Sort the list of blocks by priority
        sorted_indices = np.argsort(_block_priorities)

        unscheduled_blocks = []
        # Compute the optimal observation time in priority order
        for i in sorted_indices:
            b = blocks[i]

            # Compute possible observing times by combining object constraints
            # with the master schedule mask
            constraint_scores = np.zeros(len(times))
            for constraint in b._all_constraints:
                applied_constraint = constraint(self.observer, [b.target],
                                                times=times)
                applied_score = np.asarray(applied_constraint[0], np.float32)
                constraint_scores = constraint_scores + applied_score

            # Add up the applied constraints to prioritize the best blocks
            # And then remove any times that are already scheduled
            constraint_scores[is_open_time == False] = 0

            # Select the most optimal time
            _is_scheduled = False
            total_duration = b.duration + self.slew_time
            if np.all(constraint_scores == 0):
                # No further calculation if no times meet the constraints
                _is_scheduled = False

            else:
                # calculate the number of time slots needed for this exposure
                _stride_by = np.int(np.ceil(total_duration / time_resolution))

                # Stride the score arrays by that number
                _strided_scores = stride_array(constraint_scores, _stride_by)

                # Collapse the sub-arrays
                # (run them through scorekeeper again? Just add them?
                # If there's a zero anywhere in there, def. have to skip)
                good = np.all(_strided_scores > 1e-5, axis=1)
                sum_scores = np.zeros(len(_strided_scores))
                sum_scores[good] = np.sum(_strided_scores[good], axis=1)

                # If an optimal block is available, _is_scheduled=True
                best_time_idx = np.argmax(sum_scores)
                new_start_time = times[best_time_idx]
                _is_scheduled = True

                # And remove it from the master time list
                is_open_time[best_time_idx:best_time_idx + _stride_by] = False

            if _is_scheduled is False:
                print("could not schedule", b.target.name)
                unscheduled_blocks.append(b)
                continue
#                best_time_idx = np.argmax(constraint_scores)
#                new_start_time = times[best_time_idx]
            else:
                # now assign the block itself times and add it to the schedule
                newb = b
                for j, slot in enumerate(self.schedule.slots):
                    if slot.start <= new_start_time and slot.end >= new_start_time + total_duration:
                        slot_index = j
                newb.constraints = b._all_constraints
                self.schedule.insert_slot(slot_index, new_start_time, newb)

        return self.schedule


class Transitioner(object):
    """
    A class that defines how to compute transition times from one block to
    another.
    """
    u.quantity_input(slew_rate=u.deg/u.second)

    def __init__(self, slew_rate=None, instrument_reconfig_times=None):
        """
        Parameters
        ----------
        slew_rate : `~astropy.units.Quantity` with angle/time units
            The slew rate of the telescope
        instrument_reconfig_times : dict of dicts or None
            If not None, gives a mapping from property names to another
            dictionary. The second dictionary maps 2-tuples of states to the
            time it takes to transition between those states (as an
            `~astropy.units.Quantity`).
        """
        self.slew_rate = slew_rate
        self.instrument_reconfig_times = instrument_reconfig_times

    def __call__(self, oldblock, newblock, start_time, observer):
        """
        Determines the amount of time needed to transition from one observing
        block to another.  This uses the parameters defined in
        ``self.instrument_reconfig_times``.

        Parameters
        ----------
        oldblock : `~astroplan.scheduling.ObservingBlock` or None
            The initial configuration/target
        newblock : `~astroplan.scheduling.ObservingBlock` or None
            The new configuration/target to transition to
        start_time : `~astropy.time.Time`
            The time the transition should start
        observer : `astroplan.Observer`
            The observer at the time

        Returns
        -------
        transition : `~astroplan.scheduling.TransitionBlock` or None
            A transition to get from ``oldblock`` to ``newblock`` or `None` if
            no transition is necessary
        """
        components = {}
        if self.slew_rate is not None:
            # use the constraints cache for now, but should move that machinery
            # to observer
            from .constraints import _get_altaz
            from astropy.time import Time

            aaz = _get_altaz(Time([start_time]), observer,
                             [oldblock.target, newblock.target])['altaz']
            # TODO: make this [0] unnecessary by fixing _get_altaz to behave well in scalar-time case
            sep = aaz[0].separation(aaz[1])[0]

            components['slew_time'] = sep / self.slew_rate
        if self.instrument_reconfig_times is not None:
            components.update(self.compute_instrument_transitions(oldblock, newblock))

        if components:
            return TransitionBlock(components, start_time)
        else:
            return None

        def compute_instrument_transitions(self, oldblock, newblock):
            components = {}
            for conf_name, old_conf in oldblock.configuration.items():
                if conf_name in newblock:
                    conf_times = self.instrument_reconfig_times.get(conf_name,
                                                                    None)
                    if conf_times is not None:
                        new_conf = newblock[conf_name]
                        ctime = conf_times.get((old_conf, new_conf), None)
                        if ctime is not None:
                            s = '{0}:{1} to {2}'.format(conf_name, old_conf,
                                                        new_conf)
                            components[s] = ctime
            return components
