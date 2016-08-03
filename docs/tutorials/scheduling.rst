.. _scheduling_tutorial:

***************************
Scheduling an Observing Run
***************************

Contents
========

* :ref:`scheduling-defining_targets`

* :ref:`scheduling-creating_constraints_and_observing_blocks`


.. _scheduling-defining_targets:

Defining Targets
================

Say we want to observe Vega, Arcturus, Deneb, and Polaris, as well as a few
other stars in Ursa Minor. For this we will be using the Apache Point Observatory.

First we define our `~astroplan.Observer` object (where we are observing from):

.. code-block:: python

    >>> from astroplan import Observer

    >>> apo = Observer.at_site('apo)

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

    >>> coords = SkyCoord('05h28m51.016s', '89h26m59.58s')
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
    >>> end_time = Time('20150=-06-16 14:00')

:ref:`Return to Top <scheduling_tutorial>`

.. _scheduling-creating_constraints_and_observing_blocks:

Creating Constraints and Observing Blocks
=========================================

We want to make sure that our targets will be high in the sky while observed
and that they will be observed during the night. We don't want any object to
be observed at an airmass greater than 2.5, but we prefer a better airmass.
Usually constraints are boolean, met or not met, but with `boolean_constraint = False`
the constraint will output floats that get better closer to the ideal.

.. code-block:: python

    >>> from astroplan.constraints import AtNightConstraint, AirmassConstraint

    >>> global_constraints = [AirmassConstraint(max = 2.5, boolean_constraint = False),
                       AtNightConstraint()]

Now that we have constraints that we will apply to every target, we need to
create `~astropy.ObservingBlock`s for each target. An observing block needs
a target, a duration, and a priority; configuration information can also be
given (i.e. filter, instrument, etc.). 
