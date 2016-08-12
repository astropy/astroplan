.. _scheduling_tutorial:

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

    >>> Deneb = FixedTarget.from_name('deneb')
    >>> M13 = FixedTarget.from_name('m13')

    >>> Deneb
    <FixedTarget "deneb" at SkyCoord (ICRS): (ra, dec) in deg (310.35797975, 45.28033881)>

    >>> M13
    <FixedTarget "m13" at SkyCoord (ICRS): (ra, dec) in deg (250.423475, 36.4613194)>

We also need to define when we will be observing the targets, `~astropy.time.Time`
objects for the start and end of our observing window. The half night goes from 7PM
local time to 1AM local time, in UTC this will be from 2AM to 8AM.

.. code-block:: python

    >>> from astropy.time import Time

    >>> start_time = Time('2016-06-07 02:00')
    >>> end_time = Time('2016-06-07 08:00')

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
    ...                       AtNightConstraint()]

Now that we have constraints that we will apply to every target, we need to
create an   `~astroplan.ObservingBlock` for each target+configuration
combination. An observing block needs a target, a duration, and a priority;
configuration information can also be given (i.e. filter, instrument, etc.).
For each filter we want 17 exposures per target (100 seconds for M13 and 60
seconds for Deneb) and the instrument has a read-out time of 20 seconds.
We also need to consider that the moon can be a problem for observations
in the B filter, so we want the moon to be down or dim for the B filter.

.. code-block:: python

    >>> from astroplan import ObservingBlock
    >>> from astroplan.constraints import MoonIlluminationConstraint
    >>> from astropy import units as u

    >>> rot = 20 * u.second
    >>> blocks = []

    >>> for filter in ['B', 'G', 'R']:
    ...     if filter == 'B':
    ...         constraints = [MoonIlluminationConstraint(max = 0.25)]
    ...     else:
    ...         constraints = None
    ...     # M13 is the science target, so I will give it priority=0, and deneb priority=1
    ...     blocks.append(ObservingBlock.from_exposures(Deneb, 1, 60*u.second, 17, rot,
    ...                                                 configuration = {'filter': filter},
    ...                                                 constraints = constraints))
    ...     blocks.append(ObservingBlock.from_exposures(M13, 0, 100*u.second, 17, rot,
    ...                                                 configuration = {'filter': filter},
    ...                                                 constraints = constraints))

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

Now all we have left is to initialize the scheduler, and run it on our
list of blocks. There are currently two schedulers to chose from in
astroplan.

The first is a sequential scheduler. It starts at the start_time and
scores each block (constraints and target) at that time and then
schedules it, it then moves to where the first observing block stops
and repeats the scoring and scheduling on the remaining blocks.

.. code-block:: python

    >>> from astroplan.scheduling import SequentialScheduler

    >>> seq_scheduler = SequentialScheduler(start_time, end_time,
    ...                                     constraints = global_constraints,
    ...                                     observer = apo,
    ...                                     transitioner = transitioner)

    >>> sequential_schedule = seq_scheduler(blocks) # doctest: +SKIP

The second is a priority scheduler. It sorts the blocks by their
priority (multiple blocks with the same priority will stay in the
order they were in), then schedules them one-by-one at the best
time for that block (highest score).

.. code-block:: python

    >>> from astroplan.scheduling import PriorityScheduler

    >>> prior_scheduler = PriorityScheduler(start_time, end_time,
    ...                                     constraints = global_constraints,
    ...                                     observer = apo,
    ...                                     transitioner = transitioner)

    >>> priority_schedule = prior_scheduler(blocks) # doctest: +SKIP

Now that you have a schedule there are a few ways of viewing it.
One way is to have it print a table where you can show, or hide,
unused time and transitions with ``show_transitions`` and
``show_unused`` (default is showing transitions and not unused).

.. code-block:: python

    >>> sequential_schedule.to_table() # doctest: +SKIP

The other way is to plot the schedule against the airmass of the
targets.

.. code-block:: python

    >>> from astroplan.plots import plot_schedule_airmass
    >>> import matplotlib.pyplot as plt

    >>> plt.figure(figsize = (14,6)) # doctest: +SKIP
    >>> plot_schedule_airmass(priority_schedule) # doctest: +SKIP
    >>> plt.legend(loc = 1) # doctest: +SKIP
    >>> plt.show() # doctest: +SKIP

.. plot::

    import astropy.units as u
    from astropy.time import Time
    from astroplan import (Observer, FixedTarget, ObservingBlock, Transitioner, PriorityScheduler)
    from astroplan.constraints import AtNightConstraint, AirmassConstraint, MoonIlluminationConstraint
    from astroplan.plots import plot_schedule_airmass
    import matplotlib.pyplot as plt

    Deneb = FixedTarget.from_name('deneb')
    M13 = FixedTarget.from_name('m13')

    start_time = Time('2016-06-07 02:00')
    end_time = Time('2016-06-07 08:00')
    apo = Observer.at_site('apo')

    global_constraints = [AirmassConstraint(max = 3, boolean_constraint = False),
                          AtNightConstraint()]
    rot = 20 * u.second
    blocks = []
    for filter in ['B', 'G', 'R']:
        if filter == 'B':
            constraints = [MoonIlluminationConstraint(max = 0.25)]
        else:
            constraints = None
        blocks.append(ObservingBlock.from_exposures(Deneb, 1, 60*u.second, 17, rot,
                                                    configuration = {'filter': filter},
                                                    constraints = constraints))
        blocks.append(ObservingBlock.from_exposures(M13, 0, 100*u.second, 17, rot,
                                                    configuration = {'filter': filter},
                                                    constraints = constraints))


    transitioner = Transitioner(.8*u.deg/u.second,
                                {'filter':{('B','G'): 10*u.second,
                                           ('G','R'): 10*u.second,
                                           'default': 30*u.second}})

    prior_scheduler = PriorityScheduler(start_time, end_time, constraints = global_constraints,
                                        observer = apo, transitioner = transitioner)
    priority_schedule = prior_scheduler(blocks)

    plt.figure(figsize = (14,6))
    plot_schedule_airmass(priority_schedule)
    plt.legend(loc=1)
    plt.show()