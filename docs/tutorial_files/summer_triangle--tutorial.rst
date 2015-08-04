:orphan:

.. _summer_triangle_tutorial:

*****************************
Observing the Summer Triangle
*****************************

Defining Objects
================

Say we want to look at the Summer Triangle (Altair, Deneb, and Vega) using the
Subaru Telescope.

First, we define our `Observer` object::

    import astropy.units as u
    from astropy.coordinates import EarthLocation
    from pytz import timezone
    from astroplan import Observer

    longitude = '-155d28m48.900s'
    latitude = '+19d49m42.600s'
    elevation = 4163 * u.m
    location = EarthLocation.from_geodetic(longitude, latitude, elevation)

    subaru = Observer(name='Subaru Telescope',
                   location=location,
                   timezone=timezone('US/Hawaii'),
                   description="Subaru Telescope on Mauna Kea, Hawaii")

Then, we define our `Target` objects (`FixedTarget`'s in this case, since the
Summer Triangle is "fixed" with respect to the celestial sphere)::

    from astropy.coordinates import SkyCoord
    from astroplan import FixedTarget

    coordinates = SkyCoord('19h50m47.6s', '+08d52m12.0s', frame='icrs')
    altair = FixedTarget(name='Altair', coord=coordinates)

    coordinates = SkyCoord('18h36m56.5s', '+38d47m06.6s', frame='icrs')
    vega = FixedTarget(name='Vega', coord=coordinates)

    coordinates = SkyCoord('20h41m25.9s', '+45d16m49.3s', frame='icrs')
    deneb = FixedTarget(name='Deneb', coord=coordinates)

We also have to define a `Time` (in UTC) at which we wish to observe.  Here, I
pick 2AM local time, which is noon UTC during the summer::

    from astropy.time import Time

    time = Time('2015-06-16 12:00:00')

Observable?
===========

Next, it would be handy to know if our targets are visible from Subaru at the
time we settled on.  In other words--are they above the horizon while the Sun
is down?

.. code-block:: ipython

    In [162]: subaru.target_is_up(time, altair)
    Out[162]: array(True, dtype=bool)

    In [163]: subaru.target_is_up(time, vega)
    Out[163]: array(True, dtype=bool)

    In [164]: subaru.target_is_up(time, deneb)
    Out[164]: array(True, dtype=bool)

...They are!

Let's also pretend I'm not sure if the Sun is down at this time:

.. code-block:: ipython

    In [165]: subaru.is_night(time)
    Out[165]: array([ True], dtype=bool)

...It is!

However, we may want to find a window of time for tonight during which all
three of our targets are above the horizon *and* the Sun is below the horizon
(let's worry about light pollution from the Moon later).

Let's define the window of time during which all targets are above the horizon::

    altair_rise = subaru.target_rise_time(time, altair)
    altair_set = subaru.target_set_time(time, altair)

    vega_rise = subaru.target_rise_time(time, vega)
    vega_set = subaru.target_set_time(time, vega)

    deneb_rise = subaru.target_rise_time(time, deneb)
    deneb_set = subaru.target_set_time(time, deneb)

    import numpy as np

    all_up_start = np.max([altair_rise, vega_rise, deneb_rise])
    all_up_end = np.min([altair_set, vega_set, deneb_set])

Now, let's find sunset and sunrise for tonight (and confirm that they are
indeed those for tonight):

.. code-block:: ipython

    In [167]: sunset_tonight = subaru.sun_set_time(time, which='nearest')

    In [168]: sunset_tonight.iso
    Out[168]: '2015-06-16 04:59:12.610'

This is '2015-06-15 18:49:12.610' US/Hawaii.

.. code-block:: ipython

    In [168]: sunrise_tonight = subaru.sun_rise_time(time, which='nearest')

    In [169]: sunrise_tonight.iso
    Out[169]: '2015-06-16 15:47:36.466'

This is '2015-06-15 05:47:36.466' US/Hawaii.

Sunset and sunrise check out, so now we define the limits of our observation
window:

.. code-block:: ipython

    In [169]: start = np.max([sunset_tonight, all_up_start])

    In [170]: start.iso
    Out[170]: '2015-06-16 06:23:40.991'

    In [171]: end = np.min([sunrise_tonight, all_up_end])

    In [172]: end.iso
    Out[172]: '2015-06-16 15:47:36.466'

So, our targets will be visible (as we've defined it above) from
'2015-06-15 20:23:40.991' to '2015-06-16 05:47:36.466' US/Hawaii.  Depending on
our observation goals, this window of time may be good enough for preliminary
planning, or we may want to optimize our observational conditions.  If the
latter is the case, go on to Optimal Observation Time.

Optimal Observation Time
========================

There are a few things we can look at to find the best time to observe our
targets on a given night.

Airmass
-------

To get a general idea of our targets' airmass on the night of observation, we
can make a plot::

    from astroplan.plots import plot_airmass
    import matplotlib.pyplot as plt
    %matplotlib inline

    plot_airmass(altair, subaru, time)
    plot_airmass(vega, subaru, time)
    plot_airmass(deneb, subaru, time)

    plt.legend(bbox_to_anchor=(1, 1), loc=2)
    plt.show()

We want a minimum airmass when observing, and it looks like sometime between
9:00 and 15:00 UTC (or 23:00 on the 15th to 5:00 on the 16th, US/Hawaii) would
be the best time to observe all three targets.

However, if we want to define a more specific time window based on airmass, we
can calculate this quantity directly.

To get airmass measurements, we have to go through the `altaz` frame:

.. code-block:: ipython

    In [172]: subaru.altaz(time, altair).secz
    Out[172]: 1.030256

    In [173]: subaru.altaz(time, vega).secz
    Out[173]: 1.0690128

    In [174]: subaru.altaz(time, deneb).secz
    Out[174]: 1.1677464

Parallactic Angle
-----------------

To get a general idea of our targets' parallactic angle on the night of
observation, we can make a plot::

    from astroplan.plots import plot_parallactic

    plot_parallactic(altair, subaru, time)
    plot_parallactic(vega, subaru, time)
    plot_parallactic(deneb, subaru, time)

    plt.legend(bbox_to_anchor=(1, 1), loc=2)
    plt.show()

We can also calculate this quantity directly:

.. code-block:: ipython

    In [176]: subaru.parallactic_angle(time, altair)
    Out[176]: −0.640582rad

    In [177]: subaru.parallactic_angle(time, vega)
    Out[177]: −0.465298rad

    In [178]: subaru.parallactic_angle(time, deneb)
    Out[178]: 0.729871rad
