.. _scheduling_tutorial:

.. doctest-skip-all

***************************
Scheduling an Observing Run
***************************

.. note::

    Some terms used have been defined in :ref:`Terminology <terminology>`.

Contents
========

* :ref:`scheduling-defining_targets`

* :ref:`scheduling-creating_constraints_and_observing_blocks`

* :ref:`scheduling-creating_a_transitioner`

* :ref:`scheduling-scheduling`

* :ref:`scheduling-user_defined_schedulers`

.. _scheduling-defining_targets:

Defining Targets
================

We want to observe Deneb and M13 in the B, V and R filters. We are scheduled
for the first half-night on July 6 2016 and want to know the order we should
schedule the targets in.

First we define our `~astroplan.Observer` object (where we are observing from):

.. code-block:: python

    >>> from astroplan import Observer

    >>> apo = Observer.at_site('apo')

Now we want to define our list of targets (`~astroplan.FixedTarget` objects),
any object that is in SIMBAD can be called by an identifier.

.. code-block:: python

    >>> from astroplan import FixedTarget

    >>> Deneb = FixedTarget.from_name('Deneb')
    >>> M13 = FixedTarget.from_name('M13')

    >>> Deneb
    <FixedTarget "Deneb" at SkyCoord (ICRS): (ra, dec) in deg (310.35797975, 45.28033881)>

    >>> M13
    <FixedTarget "M13" at SkyCoord (ICRS): (ra, dec) in deg (250.423475, 36.4613194)>

We also need to define bounds within which our blocks will be scheduled
using `~astropy.time.Time` objects. Our bounds will be from the noon
before our observation, to the noon after (19:00 UTC). Later we will
account for only being able to use the first half of the night.

.. code-block:: python

    >>> from astropy.time import Time

    >>> noon_before = Time('2016-07-06 19:00')
    >>> noon_after = Time('2016-07-07 19:00')

:ref:`Return to Top <scheduling_tutorial>`

.. _scheduling-creating_constraints_and_observing_blocks:

Creating Constraints and Observing Blocks
=========================================

An in-depth tutorial on creating and using constraints can be found in
the :ref:`constraint tutorial <constraints>`.

Constraints, when evaluated, take targets and times, and give scores that
indicate how well the combination of target and time fulfill the constraint.
We want to make sure that our targets will be high in the sky while observed
and that they will be observed during the night. We don't want any object to
be observed at an airmass greater than 3, but we prefer a better airmass.
Usually constraints scores are boolean, but with ``boolean_constraint = False``
the constraint will output floats instead, indicated when it is closer to ideal.

.. code-block:: python

    >>> from astroplan.constraints import AtNightConstraint, AirmassConstraint

    >>> global_constraints = [AirmassConstraint(max = 3, boolean_constraint = False),
    ...                       AtNightConstraint.twilight_civil()]

Now that we have constraints that we will apply to every target, we need to
create an   `~astroplan.ObservingBlock` for each target+configuration
combination. An observing block needs a target, a duration, and a priority;
configuration information can also be given (i.e. filter, instrument, etc.).
For each filter we want 16 exposures per target (100 seconds for M13 and 60
seconds for Deneb) and the instrument has a read-out time of 20 seconds.
The half night goes from 7PM local time to 1AM local time, in UTC this will
be from 2AM to 8AM, so we use `~astroplan.constraints.TimeConstraint`.

.. code-block:: python

    >>> from astroplan import ObservingBlock
    >>> from astroplan.constraints import TimeConstraint
    >>> from astropy import units as u

    >>> rot = 20 * u.second
    >>> blocks = []

    >>> first_half_night = TimeConstraint(Time('2016-07-07 02:00'), Time('2016-07-07 08:00'))
    >>> for priority, bandpass in enumerate(['B', 'G', 'R']):
    ...     # We want each filter to have separate priority (so that target
    ...     # and reference are both scheduled)
    ...     blocks.append(ObservingBlock.from_exposures(Deneb, priority, 60*u.second, 16, rot,
    ...                                                 configuration = {'filter': bandpass},
    ...                                                 constraints = [first_half_night]))
    ...     blocks.append(ObservingBlock.from_exposures(M13, priority, 100*u.second, 16, rot,
    ...                                                 configuration = {'filter': bandpass},
    ...                                                 constraints = [first_half_night]))

.. _scheduling-creating_a_transitioner:

Creating a Transitioner
=======================

Now that we have observing blocks, we need to define how the telescope
transitions between them. The first parameter needed is the slew_rate
of the telescope (degrees/second) and the second is a dictionary that
tells how long it takes to transition between two configurations. You
can also give a default duration if you aren't able to give one for
each pair of configurations.

.. code-block:: python

    >>> from astroplan.scheduling import Transitioner

    >>> transitioner = Transitioner(.8*u.deg/u.second,
    ...                             {'filter':{('B','G'): 10*u.second,
    ...                                        ('G','R'): 10*u.second,
    ...                                        'default': 30*u.second}})

The transitioner now knows that it takes 10 seconds to go from 'B' to 'G',
or from 'G' to 'R' but has to use the default transition time of 30 seconds
for any other transition between filters. Non-transitions, like 'g' to 'g',
will not take any time though.

.. _scheduling-scheduling:

Scheduling
==========

Now all we have left is to initialize the scheduler, input our list
of blocks and the schedule to put them in. There are currently two
schedulers to chose from in astroplan.

The first is a sequential scheduler. It starts at the start_time and
scores each block (constraints and target) at that time and then
schedules it, it then moves to where the first observing block stops
and repeats the scoring and scheduling on the remaining blocks.

.. code-block:: python

    >>> from astroplan.scheduling import SequentialScheduler
    >>> from astroplan.scheduling import Schedule

    >>> seq_scheduler = SequentialScheduler(constraints = global_constraints,
    ...                                     observer = apo,
    ...                                     transitioner = transitioner)
    >>> sequential_schedule = Schedule(noon_before, noon_after)

    >>> seq_scheduler(blocks, sequential_schedule)

The second is a priority scheduler. It sorts the blocks by their
priority (multiple blocks with the same priority will stay in the
order they were in), then schedules them one-by-one at the best
time for that block (highest score).

.. code-block:: python

    >>> from astroplan.scheduling import PriorityScheduler

    >>> prior_scheduler = PriorityScheduler(constraints = global_constraints,
    ...                                     observer = apo,
    ...                                     transitioner = transitioner)
    >>> priority_schedule = Schedule(noon_before, noon_after)

    >>> prior_scheduler(blocks, priority_schedule)

Now that you have a schedule there are a few ways of viewing it.
One way is to have it print a table where you can show, or hide,
unused time and transitions with ``show_transitions`` and
``show_unused`` (default is showing transitions and not unused).

.. code-block:: python

    >>> priority_schedule.to_table()
         target         start time (UTC)         end time (UTC)     duration (minutes)        ra            dec         configuration
         str15               str23                   str23               float64            str32          str32            object
    --------------- ----------------------- ----------------------- ------------------ --------------- -------------- -----------------
                M13 2016-07-07 03:49:20.019 2016-07-07 04:21:20.019               32.0   250d25m24.51s 36d27m40.7498s   {'filter': 'R'}
    TransitionBlock 2016-07-07 04:21:20.019 2016-07-07 04:22:00.019     0.666666666667                                ['filter:R to B']
                M13 2016-07-07 04:25:20.021 2016-07-07 04:57:20.021               32.0   250d25m24.51s 36d27m40.7498s   {'filter': 'B'}
    TransitionBlock 2016-07-07 04:57:20.021 2016-07-07 04:57:40.021     0.333333333333                                ['filter:B to G']
                M13 2016-07-07 04:57:40.021 2016-07-07 05:29:40.021               32.0   250d25m24.51s 36d27m40.7498s   {'filter': 'G'}
    TransitionBlock 2016-07-07 05:29:40.021 2016-07-07 05:31:00.021      1.33333333333                                ['filter:G to R']
              Deneb 2016-07-07 06:44:00.026 2016-07-07 07:05:20.026      21.3333333333 310d21m28.7271s 45d16m49.2197s   {'filter': 'R'}
    TransitionBlock 2016-07-07 07:05:20.026 2016-07-07 07:06:00.026     0.666666666667                                ['filter:R to G']
              Deneb 2016-07-07 07:09:20.027 2016-07-07 07:30:40.027      21.3333333333 310d21m28.7271s 45d16m49.2197s   {'filter': 'G'}
    TransitionBlock 2016-07-07 07:30:40.027 2016-07-07 07:31:20.027     0.666666666667                                ['filter:G to B']
              Deneb 2016-07-07 07:34:40.028 2016-07-07 07:56:00.028      21.3333333333 310d21m28.7271s 45d16m49.2197s   {'filter': 'B'}

The other way is to plot the schedule against the airmass of the
targets.

.. code-block:: python

    >>> from astroplan.plots import plot_schedule_airmass
    >>> import matplotlib.pyplot as plt

    >>> plt.figure(figsize = (14,6))
    >>> plot_schedule_airmass(priority_schedule)
    >>> plt.legend(loc = "upper right")
    >>> plt.show()

.. plot::

    # first import everything we will need for the scheduling
    import astropy.units as u
    from astropy.time import Time
    from astroplan import (Observer, FixedTarget, ObservingBlock, Transitioner, PriorityScheduler,
                           Schedule)
    from astroplan.constraints import AtNightConstraint, AirmassConstraint, TimeConstraint
    from astroplan.plots import plot_schedule_airmass
    import matplotlib.pyplot as plt

    # Now we define the targets, observer, start time, and end time of the schedule.
    Deneb = FixedTarget.from_name('Deneb')
    M13 = FixedTarget.from_name('M13')

    noon_before = Time('2016-07-06 19:00')
    noon_after = Time('2016-07-07 19:00')
    apo = Observer.at_site('apo')

    # Then define the constraints (global and specific) and make a list of the
    # observing blocks that you want scheduled
    global_constraints = [AirmassConstraint(max = 3, boolean_constraint = False),
                          AtNightConstraint.twilight_civil()]
    rot = 20 * u.second
    blocks = []
    first_half_night = TimeConstraint(Time('2016-07-07 02:00'), Time('2016-07-07 08:00'))
    for priority, bandpass in enumerate(['B', 'G', 'R']):
        # We want each filter to have separate priority (so that target
        # and reference are both scheduled)
        blocks.append(ObservingBlock.from_exposures(Deneb, priority, 60*u.second, 16, rot,
                                                    configuration = {'filter': bandpass},
                                                    constraints = [first_half_night]))
        blocks.append(ObservingBlock.from_exposures(M13, priority, 100*u.second, 16, rot,
                                                    configuration = {'filter': bandpass},
                                                    constraints = [first_half_night]))

    # Define how the telescope transitions between the configurations defined in the
    # observing blocks (target, filter, instrument, etc.).
    transitioner = Transitioner(.8*u.deg/u.second,
                                {'filter':{('B','G'): 10*u.second,
                                           ('G','R'): 10*u.second,
                                           'default': 30*u.second}})

    # Initialize the scheduler
    prior_scheduler = PriorityScheduler(constraints = global_constraints,
                                        observer = apo, transitioner = transitioner)
    # Create a schedule for the scheduler to insert the blocks into, and run the scheduler
    priority_schedule = Schedule(noon_before, noon_after)
    prior_scheduler(blocks, priority_schedule)

    # To get a plot of the airmass vs where the blocks were scheduled
    plt.figure(figsize = (14,6))
    plot_schedule_airmass(priority_schedule)
    plt.tight_layout()
    plt.legend(loc="upper right")
    plt.show()

There is a lot of unfilled space in our schedule currently. We can
fill that time with more observations of our targets by calling our
scheduler using the same blocks and the schedule we already added to.

.. code-block:: python

    >>> prior_schedule(blocks, priority_schedule)
    >>> plt.figure(figsize = (14,6))
    >>> plot_schedule_airmass(priority_schedule)
    >>> plt.legend(loc = "upper right")
    >>> plt.show()

.. plot::

    # first import everything we will need for the scheduling
    import astropy.units as u
    from astropy.time import Time
    from astroplan import (Observer, FixedTarget, ObservingBlock, Transitioner, PriorityScheduler,
                           Schedule)
    from astroplan.constraints import AtNightConstraint, AirmassConstraint, TimeConstraint
    from astroplan.plots import plot_schedule_airmass
    import matplotlib.pyplot as plt

    # Now we define the targets, observer, start time, and end time of the schedule.
    Deneb = FixedTarget.from_name('Deneb')
    M13 = FixedTarget.from_name('M13')

    noon_before = Time('2016-07-06 19:00')
    noon_after = Time('2016-07-07 19:00')
    apo = Observer.at_site('apo')

    # Then define the constraints (global and specific) and make a list of the
    # observing blocks that you want scheduled
    global_constraints = [AirmassConstraint(max = 3, boolean_constraint = False),
                          AtNightConstraint.twilight_civil()]
    rot = 20 * u.second
    blocks = []
    first_half_night = TimeConstraint(Time('2016-07-07 02:00'), Time('2016-07-07 08:00'))
    for priority, bandpass in enumerate(['B', 'G', 'R']):
        # We want each filter to have separate priority (so that target
        # and reference are both scheduled)
        blocks.append(ObservingBlock.from_exposures(Deneb, priority, 60*u.second, 16, rot,
                                                    configuration = {'filter': bandpass},
                                                    constraints = [first_half_night]))
        blocks.append(ObservingBlock.from_exposures(M13, priority, 100*u.second, 16, rot,
                                                    configuration = {'filter': bandpass},
                                                    constraints = [first_half_night]))

    # Define how the telescope transitions between the configurations defined in the
    # observing blocks (target, filter, instrument, etc.).
    transitioner = Transitioner(.8*u.deg/u.second,
                                {'filter':{('B','G'): 10*u.second,
                                           ('G','R'): 10*u.second,
                                           'default': 30*u.second}})

    # Initialize the scheduler
    prior_scheduler = PriorityScheduler(constraints = global_constraints,
                                        observer = apo, transitioner = transitioner)
    # Create a schedule for the scheduler to insert the blocks into, and run the scheduler
    priority_schedule = Schedule(noon_before, noon_after)
    prior_scheduler(blocks, priority_schedule)
    # add more of the blocks into the schedule
    prior_scheduler(blocks, priority_schedule)

    # To get a plot of the airmass vs where the blocks were scheduled
    plt.figure(figsize = (14,6))
    plot_schedule_airmass(priority_schedule)
    plt.tight_layout()
    plt.legend(loc="upper right")
    plt.show()

We want to check if there is any way that we could observe Alpha
Centauri A as well during our time slot. So we create a new block
for it with priority over the others, add it to our list of blocks
and run the priority scheduler again.

.. code-block:: python

    >>> alf_cent = FixedTarget.from_name('Alpha Centauri A')
    >>> blocks.append(ObservingBlock(alf_cent, 20*u.minute, -1))
    >>> schedule = Schedule(start_time, end_time)
    >>> prior_scheduler(blocks, schedule)

    >>> plt.figure(figsize = (14,6))
    >>> plot_schedule_airmass(priority_schedule)
    >>> plt.legend(loc = "upper right")
    >>> plt.show()

.. plot::

    # first import everything we will need for the scheduling
    import astropy.units as u
    from astropy.time import Time
    from astroplan import (Observer, FixedTarget, ObservingBlock, Transitioner, PriorityScheduler,
                           Schedule)
    from astroplan.constraints import AtNightConstraint, AirmassConstraint, TimeConstraint
    from astroplan.plots import plot_schedule_airmass
    import matplotlib.pyplot as plt

    # Now we define the targets, observer, start time, and end time of the schedule.
    Deneb = FixedTarget.from_name('Deneb')
    M13 = FixedTarget.from_name('M13')

    noon_before = Time('2016-07-06 19:00')
    noon_after = Time('2016-07-07 19:00')
    apo = Observer.at_site('apo')

    # Then define the constraints (global and specific) and make a list of the
    # observing blocks that you want scheduled
    global_constraints = [AirmassConstraint(max = 3, boolean_constraint = False),
                          AtNightConstraint.twilight_civil()]
    rot = 20 * u.second
    blocks = []
    first_half_night = TimeConstraint(Time('2016-07-07 02:00'), Time('2016-07-07 08:00'))
    for priority, bandpass in enumerate(['B', 'G', 'R']):
        # We want each filter to have separate priority (so that target
        # and reference are both scheduled)
        blocks.append(ObservingBlock.from_exposures(Deneb, priority, 60*u.second, 16, rot,
                                                    configuration = {'filter': bandpass},
                                                    constraints = [first_half_night]))
        blocks.append(ObservingBlock.from_exposures(M13, priority, 100*u.second, 16, rot,
                                                    configuration = {'filter': bandpass},
                                                    constraints = [first_half_night]))
    # add the new target's block
    alf_cent = FixedTarget.from_name('Alpha Centauri A')
    blocks.append(ObservingBlock(alf_cent, 20*u.minute, -1))

    # Define how the telescope transitions between the configurations defined in the
    # observing blocks (target, filter, instrument, etc.).
    transitioner = Transitioner(.8*u.deg/u.second,
                                {'filter':{('B','G'): 10*u.second,
                                           ('G','R'): 10*u.second,
                                           'default': 30*u.second}})

    # Initialize the scheduler
    prior_scheduler = PriorityScheduler(constraints = global_constraints,
                                        observer = apo, transitioner = transitioner)
    # Create a schedule for the scheduler to insert the blocks into, and run the scheduler
    priority_schedule = Schedule(noon_before, noon_after)
    prior_scheduler(blocks, priority_schedule)

    # To get a plot of the airmass vs where the blocks were scheduled
    plt.figure(figsize = (14,6))
    plot_schedule_airmass(priority_schedule)
    plt.tight_layout()
    plt.legend(loc="upper right")
    plt.show()

Nothing new shows up because Alpha Centauri isn't visible from APO.

.. _scheduling-user_defined_schedulers:

User-Defined Schedulers
=======================

There are many ways that targets can be scheduled with, only two of which
are currently implemented. This example will walk through the steps for
creating your own scheduler that will be compatible with the tools of
the ``scheduling`` module.

As you may have noticed above, the schedulers are assembled by making a
call to the initializer of the class (e.g. `~astroplan.scheduling.PriorityScheduler`).
Each of the schedulers is subclassed from the abstract `astroplan.scheduling.Scheduler`
class, and our custom scheduler needs to be as well.

A scheduler needs to be able to schedule observing blocks where they have a non-zero
score (i.e. they satisfy all of their constraints). For our scheduler, we will make
one that schedules ``ObservingBlocks`` at the first unoccupied place they have a score
greater than zero: a ``SimpleScheduler``. We need to include two methods, ``__init__``
and ``_make_schedule`` for it to work:

* The ``__init__`` is already defined by the super class, and accepts global constraints,
  the `~astroplan.Observer`, the `~astroplan.scheduling.Transitioner`, a ``gap_time``,
  and a ``time_resolution`` for spacing during the creation of the schedule.

* It also needs a ``_make_schedule`` to do the heavy lifting. This takes a list of
  `~astroplan.scheduling.ObservingBlock` objects and a `~astroplan.scheduling.Schedule`
  object to input them into. This method needs to be able to check whether a
  block can be scheduled in a given spot, and be able to insert it into the
  schedule once a suitable spot has been found.

Here's the ``SimpleScheduler`` implementation::

    from astroplan.scheduling import Scheduler, Scorer
    from astroplan.utils import time_grid_from_range
    from astroplan.constraints import AltitudeConstraint
    from astropy import units as u

    import numpy as np

    class SimpleScheduler(Scheduler):
        """
        schedule blocks randomly
        """
        def __init__(self, *args, **kwargs):
            super(SimpleScheduler, self).__init__(*args, **kwargs)

        def _make_schedule(self, blocks):
            # gather all the constraints on each block into a single attribute
            for b in blocks:
                if b.constraints is None:
                    b._all_constraints = self.constraints
                else:
                    b._all_constraints = self.constraints + b.constraints
                # to make sure the scheduler has some constraint to work off of
                # and to prevent scheduling of targets below the horizon
                if b._all_constraints is None:
                    b._all_constraints = [AltitudeConstraint(min=0*u.deg)]
                    b.constraints = [AltitudeConstraint(min=0*u.deg)]
                elif not any(isinstance(c, AltitudeConstraint) for c in b._all_constraints):
                    b._all_constraints.append(AltitudeConstraint(min=0*u.deg))
                b.observer = self.observer

            # before we can schedule, we need to know where blocks meet the constraints
            scorer = Scorer(blocks,self.observer, self.schedule, global_constraints=self.constraints)
            score_array = scorer.create_score_array(self.time_resolution)
            # now we have an array with the scores for all of the blocks at intervals of time_resolution

            # we want to start from the beginning and start scheduling
            current_time = self.schedule.start_time
            while current_time < self.schedule.end_time:
                scheduled = False
                i=0
                while i < len(blocks) and scheduled is False:
                    block = blocks[i]
                    # the schedule starts with only 1 slot
                    if len(self.schedule.slots) > 1:
                        # make a transition between the last scheduled block and this one
                        transition = self.transitioner(schedule.observing_blocks[-1], block,
                                                       current_time, self.observer)
                        test_time = current_time + transition.duration
                    else:
                        test_time = current_time
                    # how far from the start is the time we are testing
                    start_idx = int((test_time - self.schedule.start_time)/self.time_resolution)
                    duration_idx = int(block.duration/self.time_resolution)
                    # if any score during the block's duration would be 0, reject it
                    if any(score_array[i][start_idx:start_idx+duration_idx] == 0):
                        i +=1
                    # if all of the scores are >0, accept and schedule it
                    else:
                        if len(self.schedule.slots) >1:
                            self.schedule.insert_slot(current_time, transition)
                        self.schedule.insert_slot(test_time, block)
                        # advance the time, break this while loop and remove the block from the list
                        current_time = test_time + block.duration
                        scheduled = True
                        blocks.remove(block)
                # if every block failed, progress the time
                if i == len(blocks):
                    current_time += self.gap_time
            return schedule

Then to use our new scheduler, we just need to call it how we did
up above::

    >>> from astroplan.constraints import AtNightConstraint
    >>> from astroplan.scheduling import Schedule, ObservingBlock
    >>> from astroplan import FixedTarget, Observer, Transitioner
    >>> from astropy.time import Time

    >>> apo = Observer.at_site('apo')
    >>> Deneb = FixedTarget.from_name('Deneb')
    >>> M13 = FixedTarget.from_name('M13')
    >>> blocks = [ObservingBlock(Deneb, 20*u.minute, 0)]
    >>> blocks.append(ObservingBlock(M13, 20*u.minute, 0))

    >>> # this will produce 0 second transitions
    >>> transitioner = Transitioner()

    >>> schedule = Schedule(Time('2016-07-06 19:00'), Time('2016-07-07 19:00'))

    >>> scheduler = SimpleScheduler(observer = apo, transitioner = transitioner,
    ...                                 constraints = [])
    >>> scheduler(blocks, schedule)

    >>> import matplotlib.pyplot as plt
    >>> from astroplan.plots import plot_schedule_airmass
    >>> plot_schedule_airmass(schedule)
    >>> plt.legend()
    >>> plt.show()

.. plot::

    from astroplan.scheduling import Scheduler, Scorer
    from astroplan.utils import time_grid_from_range
    from astroplan.constraints import AltitudeConstraint
    from astropy import units as u

    import numpy as np

    class SimpleScheduler(Scheduler):
        """
        schedule blocks randomly
        """
        def __init__(self, *args, **kwargs):
            super(SimpleScheduler, self).__init__(*args, **kwargs)

        def _make_schedule(self, blocks):
            # gather all the constraints on each block into a single attribute
            for b in blocks:
                if b.constraints is None:
                    b._all_constraints = self.constraints
                else:
                    b._all_constraints = self.constraints + b.constraints
                # to make sure the scheduler has some constraint to work off of
                # and to prevent scheduling of targets below the horizon
                if b._all_constraints is None:
                    b._all_constraints = [AltitudeConstraint(min=0*u.deg)]
                    b.constraints = [AltitudeConstraint(min=0*u.deg)]
                elif not any(isinstance(c, AltitudeConstraint) for c in b._all_constraints):
                    b._all_constraints.append(AltitudeConstraint(min=0*u.deg))
                b.observer = self.observer

            # before we can schedule, we need to know where blocks meet the constraints
            scorer = Scorer(blocks,self.observer, self.schedule, global_constraints=self.constraints)
            score_array = scorer.create_score_array(self.time_resolution)
            # now we have an array with the scores for all of the blocks at intervals of time_resolution

            # we want to start from the beginning and start scheduling
            current_time = self.schedule.start_time
            while current_time < self.schedule.end_time:
                scheduled = False
                i=0
                while i < len(blocks) and scheduled is False:
                    block = blocks[i]
                    # the schedule starts with only 1 slot
                    if len(self.schedule.slots) > 1:
                        # make a transition between the last scheduled block and this one
                        transition = self.transitioner(schedule.observing_blocks[-1], block,
                                                       current_time, self.observer)
                        test_time = current_time + transition.duration
                    else:
                        test_time = current_time
                    # how far from the start is the time we are testing
                    start_idx = int((test_time - self.schedule.start_time)/self.time_resolution)
                    duration_idx = int(block.duration/self.time_resolution)
                    # if any score during the block's duration would be 0, reject it
                    if any(score_array[i][start_idx:start_idx+duration_idx] == 0):
                        i +=1
                    # if all of the scores are >0, accept and schedule it
                    else:
                        if len(self.schedule.slots) >1:
                            self.schedule.insert_slot(current_time, transition)
                        self.schedule.insert_slot(test_time, block)
                        # advance the time, break this while loop and remove the block from the list
                        current_time = test_time + block.duration
                        scheduled = True
                        blocks.remove(block)
                # if every block failed, progress the time
                if i == len(blocks):
                    current_time += self.gap_time
            return schedule
    from astroplan.constraints import AtNightConstraint
    from astroplan.scheduling import Schedule, ObservingBlock
    from astroplan import FixedTarget, Observer, Transitioner
    from astropy.time import Time

    apo = Observer.at_site('apo')
    Deneb = FixedTarget.from_name('Deneb')
    M13 = FixedTarget.from_name('M13')
    blocks = [ObservingBlock(Deneb, 20*u.minute, 0)]
    blocks.append(ObservingBlock(M13, 20*u.minute, 0))

    transitioner = Transitioner()
    global_constraints = [AtNightConstraint.twilight_civil()]

    schedule = Schedule(Time('2016-07-06 19:00'), Time('2016-07-07 19:00'))
    scheduler = SimpleScheduler(observer = apo, transitioner = transitioner,
                                    constraints = [])
    scheduler(blocks, schedule)

    import matplotlib.pyplot as plt
    from astroplan.plots import plot_schedule_airmass
    plot_schedule_airmass(schedule)
    plt.tight_layout()
    plt.legend()
    plt.show()
