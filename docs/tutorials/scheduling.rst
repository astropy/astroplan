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

Say we want to observe Vega, Arcturus, Deneb, and Polaris, as well as a few
other stars in Ursa Minor. For this we will be using the Apache Point Observatory.

First we define our `~astroplan.Observer` object (where we are observing from):

.. code-block:: python

    >>> from astroplan import Observer

    >>> apo = Observer.at_site('apo')

Now we want to define our list of targets (`~astroplan.FixedTarget` objects),
Most stars can be called by an identifier, but some targets may require coordinates.

.. code-block:: python

    >>> from astroplan import FixedTarget

    >>> targets = [FixedTarget.from_name('vega'),
    >>>            FixedTarget.from_name('Arcturus'),
    >>>            FixedTarget.from_name('Deneb'),
    >>>            FixedTarget.from_name('Polaris'),
    >>>            FixedTarget.from_name('eta umi'),
    >>>            FixedTarget.from_name('NLTT 42620'),
    >>>            FixedTarget.from_name('eps umi')
    >>>            ]

    >>> from astropy.coordinates import SkyCoord

    >>> coords = SkyCoord('05h28m51.016s', '89d26m59.58s')
    >>> interesting_target = FixedTarget(coords, name = 'UMi target')
    >>> targets.append(interesting_target)

    >>> targets
    [<FixedTarget "vega" at SkyCoord (ICRS): (ra, dec) in deg (279.23473479, 38.78368896)>,
     <FixedTarget "Arcturus" at SkyCoord (ICRS): (ra, dec) in deg (213.9153003, 19.18240916)>,
     <FixedTarget "Deneb" at SkyCoord (ICRS): (ra, dec) in deg (310.35797975, 45.28033881)>,
     <FixedTarget "Polaris" at SkyCoord (ICRS): (ra, dec) in deg (37.95456067, 89.26410897)>,
     <FixedTarget "eta umi" at SkyCoord (ICRS): (ra, dec) in deg (244.37619567, 75.75533014)>,
     <FixedTarget "NLTT 42620" at SkyCoord (ICRS): (ra, dec) in deg (244.587292, 75.718917)>,
     <FixedTarget "eps umi" at SkyCoord (ICRS): (ra, dec) in deg (251.49267367, 82.03725646)>,
     <FixedTarget "UMi target" at SkyCoord (ICRS): (ra, dec) in deg (82.21256667, 89.44988333)>]

We also need to define when we will be observing the targets, `~astropy.time.Time`
objects for the start and end of our observing window. Here we will take an entire
night from 8PM local time to 7AM local time, in UTC this will be from 3AM to 2PM

.. code-block:: python

    >>> from astropy.time import Time

    >>> start_time = Time('2015-06-16 03:00')
    >>> end_time = Time('2015-06-16 14:00')

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
be observed at an airmass greater than 2.5, but we prefer a better airmass.
Usually constraints scores are boolean, but with ``boolean_constraint = False``
the constraint will output floats instead, indicated when it is closer to ideal.

.. code-block:: python

    >>> from astroplan.constraints import AtNightConstraint, AirmassConstraint

    >>> global_constraints = [AirmassConstraint(max = 2.5, boolean_constraint = False),
    >>>                       AtNightConstraint()]

Now that we have constraints that we will apply to every target, we need to
create an   `~astroplan.ObservingBlock` for each target. An observing block needs
a target, a duration, and a priority; configuration information can also be
given (i.e. filter, instrument, etc.). We want all of the targets in the 'g'
filter, and also want 'r' and 'i' for our UMi target. We also want to make sure
that our observations of the UMi target occur between 10PM and 3AM local time
(5AM and 10 AM UTC).For each target we want 7 exposures (with length depending
on the target) and the instrument has a read-out time of 1 minute.

.. code-block:: python

    >>> from astroplan import ObservingBlock
    >>> from astroplan.constraints import TimeConstraint
    >>> from astropy import units as u

    >>> rot = 1 * u.minute
    >>> blocks = []
    >>> # first we will make the blocks for our UMi target
    >>> time_constraint = TimeConstraint(start_time + 2*u.hour, end_time - 4*u.hour)
    >>> for filter in ['g', 'r', 'i']:
    >>>     blocks.append(ObservingBlock.from_exposures(targets[-1], 0, 8*u.minute, 7, rot,
    >>>                                                 configuration = {'filter': filter},
    >>>                                                 constraints = [time_constraint]))
    >>> for target in targets[4:7]:
    >>>     blocks.append(ObservingBlock.from_exposures(target, 1, 4*u.minute, 7, rot,
    >>>                                                 configuration = {'filter': 'g'}))
    >>> for target in targets[:4]:
    >>>     blocks.append(ObservingBlock.from_exposures(target, 2, 2*u.minute, 7, rot,
    >>>                                                 configuration = {'filter': 'g'}))

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
    >>>                             {'filter':{('g','i'): 10*u.second,
    >>>                                        'default': 30*u.second}})

The transitioner knows that it takes 10 seconds to go from 'g' to 'i'
but has to use the default transition time of 30 seconds for any other
transition between filters. Non-transitions, like 'g' to 'g', will not
take any time though.

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
    >>>                                     constraints = global_constraints,
    >>>                                     observer = apo,
    >>>                                     transitioner = transitioner)

    >>> sequential_schedule = seq_scheduler(blocks)

The second is a priority scheduler. It sorts the blocks by their
priority (multiple blocks with the same priority will stay in the
order they were in), then schedules them one-by-one at the best
time for that block (highest score).

    >>> from astroplan.scheduling import PriorityScheduler

    >>> prior_scheduler = PriorityScheduler(start_time, end_time,
    >>>                                     constraints = global_constraints,
    >>>                                     observer = apo,
    >>>                                     transitioner = transitioner)

    >>> priority_schedule = prior_scheduler(blocks)

Now that you have a schedule there are a few ways of viewing it.
One way is to have it print a table where you can show, or hide,
unused time and transitions with ``show_transitions`` and
``show_unused`` (default is showing transitions and not unused).

.. code-block:: python

    >>> sequential_schedule.to_table()

The other way is to plot the schedule against the airmass of the
targets.

.. code-block:: python

    >>> from astroplan.plots import plot_schedule_airmass
    >>> import matplotlib.pyplot as plt

    >>> plot_schedule_airmass(priority_schedule)
    >>> plt.legend(loc=0)
    >>> plt.show()

.. plot::

    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astropy.time import Time
    from astroplan import (Observer, FixedTarget, ObservingBlock, Transitioner, PriorityScheduler)
    from astroplan.constraints import AtNightConstraint, AirmassConstraint, TimeConstraint
    from astroplan.plots import plot_schedule_airmass
    import matplotlib.pyplot as plt

    targets = [FixedTarget.from_name('vega'),
               FixedTarget.from_name('Arcturus'),
               FixedTarget.from_name('Deneb'),
               FixedTarget.from_name('Polaris'),
               FixedTarget.from_name('eta umi'),
               FixedTarget.from_name('NLTT 42620'),
               FixedTarget.from_name('eps umi')]

    coords = SkyCoord('05h28m51.016s', '89d26m59.58s')
    interesting_target = FixedTarget(coords, name = 'UMi target')
    targets.append(interesting_target)

    start_time = Time('2015-06-16 03:00')
    end_time = Time('2015-06-16 14:00')
    apo = Observer.at_site('apo')

    global_constraints = [AirmassConstraint(max = 2.5, boolean_constraint = False),
                          AtNightConstraint()]
    rot = 1 * u.minute
    blocks = []
    time_constraint = TimeConstraint(start_time + 2*u.hour, end_time - 4*u.hour)
    for filter in ['g', 'r', 'i']:
        blocks.append(ObservingBlock.from_exposures(targets[-1], 0, 8*u.minute, 7, rot,
                                                    configuration = {'filter': filter},
                                                    constraints = [time_constraint]))
    for target in targets[4:7]:
        blocks.append(ObservingBlock.from_exposures(target, 1, 4*u.minute, 7, rot,
                                                    configuration = {'filter': 'g'}))
    for target in targets[:4]:
        blocks.append(ObservingBlock.from_exposures(target, 2, 2*u.minute, 7, rot,
                                                    configuration = {'filter': 'g'}))

    transitioner = Transitioner(.8*u.deg/u.second,
                                {'filter':{('g','i'): 10*u.second,
                                           'default': 30*u.second}})

    prior_scheduler = PriorityScheduler(start_time, end_time, constraints = global_constraints,
                                        observer = apo, transitioner = transitioner)
    priority_schedule = prior_scheduler(blocks)

    plot_schedule_airmass(priority_schedule)
    plt.legend(loc=0)
    plt.show()