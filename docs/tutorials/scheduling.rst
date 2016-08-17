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

    >>> sequential_schedule.to_table()
         target         start time (UTC)         end time (UTC)     ...        ra            dec
         str15               str23                   str23          ...      str32          str32
    --------------- ----------------------- ----------------------- ... --------------- --------------
                M13 2016-07-07 02:45:00.000 2016-07-07 03:17:00.000 ...   250d25m24.51s 36d27m40.7498s
    TransitionBlock 2016-07-07 03:17:00.000 2016-07-07 03:17:30.000 ...
                M13 2016-07-07 03:17:30.000 2016-07-07 03:49:30.000 ...   250d25m24.51s 36d27m40.7498s
    TransitionBlock 2016-07-07 03:49:30.000 2016-07-07 03:50:00.000 ...
                M13 2016-07-07 03:50:00.000 2016-07-07 04:22:00.000 ...   250d25m24.51s 36d27m40.7498s
    TransitionBlock 2016-07-07 04:22:00.000 2016-07-07 04:23:26.384 ...
              Deneb 2016-07-07 04:23:26.384 2016-07-07 04:44:46.384 ... 310d21m28.7271s 45d16m49.2197s
    TransitionBlock 2016-07-07 04:44:46.384 2016-07-07 04:45:16.384 ...
              Deneb 2016-07-07 04:45:16.384 2016-07-07 05:06:36.384 ... 310d21m28.7271s 45d16m49.2197s
    TransitionBlock 2016-07-07 05:06:36.384 2016-07-07 05:07:06.384 ...
              Deneb 2016-07-07 05:07:06.384 2016-07-07 05:28:26.384 ... 310d21m28.7271s 45d16m49.2197s

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
    plt.legend(loc="upper right")
    plt.show()

Nothing new shows up because Alpha Centauri isn't visible from APO.