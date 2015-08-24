:orphan:

.. _summer_triangle_tutorial:

.. doctest-skip-all

.. todo::

    Add section on moon phases, illumination fraction, etc.

    Replace target construction with easier site name function.

    Update with constraints when available?

    Rise/set times currently return altitudes just below the horizon by a few
    arcseconds.  Need to clarify this and fix sky plot examples.

    Show users how to create/use Time objects in their own timezone.

*****************************
Observing the Summer Triangle
*****************************

Contents
========

* :ref:`summer_triangle-defining_objects`

* :ref:`summer_triangle-observable`

* :ref:`summer_triangle-optimal_observation`

* :ref:`summer_triangle-sky_charts`

.. _summer_triangle-defining_objects:

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

:ref:`Return to Top <summer_triangle_tutorial>`

.. _summer_triangle-observable:

Observable?
===========

Next, it would be handy to know if our targets are visible from Subaru at the
time we settled on.  In other words--are they above the horizon while the Sun
is down?

.. code-block:: python

    >>> subaru.target_is_up(time, altair)
    array(True, dtype=bool)

    >>> subaru.target_is_up(time, vega)
    array(True, dtype=bool)

    >>> subaru.target_is_up(time, deneb)
    array(True, dtype=bool)

...They are!

Let's also pretend I'm not sure if the Sun is down at this time:

.. code-block:: python

    >>> subaru.is_night(time)
    array([ True], dtype=bool)

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

.. code-block:: python

    >>> sunset_tonight = subaru.sun_set_time(time, which='nearest')

    >>> sunset_tonight.iso
    '2015-06-16 04:59:12.610'

This is '2015-06-15 18:49:12.610' US/Hawaii.

.. code-block:: python

    >>> sunrise_tonight = subaru.sun_rise_time(time, which='nearest')

    >>> sunrise_tonight.iso
    '2015-06-16 15:47:36.466'

This is '2015-06-16 05:47:36.466' US/Hawaii.

Sunset and sunrise check out, so now we define the limits of our observation
window:

.. code-block:: python

    >>> start = np.max([sunset_tonight, all_up_start])

    >>> start.iso
    '2015-06-16 06:23:40.991'

    >>> end = np.min([sunrise_tonight, all_up_end])

    >>> end.iso
    '2015-06-16 15:47:36.466'

So, our targets will be visible (as we've defined it above) from
'2015-06-15 20:23:40.991' to '2015-06-16 05:47:36.466' US/Hawaii.  Depending on
our observation goals, this window of time may be good enough for preliminary
planning, or we may want to optimize our observational conditions.  If the
latter is the case, go on to Optimal Observation Time.

:ref:`Return to Top <summer_triangle_tutorial>`

.. _summer_triangle-optimal_observation:

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

    plot_airmass(altair, subaru, time)
    plot_airmass(vega, subaru, time)
    plot_airmass(deneb, subaru, time)

    plt.legend(loc=1, bbox_to_anchor=(1, 1))
    plt.show()

.. plot::

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

    from astropy.coordinates import SkyCoord
    from astroplan import FixedTarget

    coordinates = SkyCoord('19h50m47.6s', '+08d52m12.0s', frame='icrs')
    altair = FixedTarget(name='Altair', coord=coordinates)

    coordinates = SkyCoord('18h36m56.5s', '+38d47m06.6s', frame='icrs')
    vega = FixedTarget(name='Vega', coord=coordinates)

    coordinates = SkyCoord('20h41m25.9s', '+45d16m49.3s', frame='icrs')
    deneb = FixedTarget(name='Deneb', coord=coordinates)

    from astropy.time import Time

    time = Time('2015-06-16 12:00:00')

    from astroplan.plots import plot_airmass
    import matplotlib.pyplot as plt

    plot_airmass(altair, subaru, time)
    plot_airmass(vega, subaru, time)
    plot_airmass(deneb, subaru, time)

    # Note that you don't need this code block to produce the plot.
    # It reduces the plot size for the documentation.
    ax = plt.gca()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height * 0.8])

    plt.legend(loc=1, bbox_to_anchor=(1.35, 1))
    plt.show()

We want a minimum airmass when observing, and it looks like sometime between
9:00 and 15:00 UTC (or 23:00 on the 15th to 5:00 on the 16th, US/Hawaii) would
be the best time to observe all three targets.

However, if we want to define a more specific time window based on airmass, we
can calculate this quantity directly.

To get airmass measurements, we have to go through the `altaz` frame:

.. code-block:: python

    >>> subaru.altaz(time, altair).secz
    1.030256

    >>> subaru.altaz(time, vega).secz
    1.0690128

    >>> subaru.altaz(time, deneb).secz
    1.1677464

Parallactic Angle
-----------------

To get a general idea of our targets' parallactic angle on the night of
observation, we can make a plot::

    from astroplan.plots import plot_parallactic

    plot_parallactic(altair, subaru, time)
    plot_parallactic(vega, subaru, time)
    plot_parallactic(deneb, subaru, time)

    plt.legend(loc=2)
    plt.show()

.. plot::

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

    from astropy.coordinates import SkyCoord
    from astroplan import FixedTarget

    coordinates = SkyCoord('19h50m47.6s', '+08d52m12.0s', frame='icrs')
    altair = FixedTarget(name='Altair', coord=coordinates)

    coordinates = SkyCoord('18h36m56.5s', '+38d47m06.6s', frame='icrs')
    vega = FixedTarget(name='Vega', coord=coordinates)

    coordinates = SkyCoord('20h41m25.9s', '+45d16m49.3s', frame='icrs')
    deneb = FixedTarget(name='Deneb', coord=coordinates)

    from astropy.time import Time

    time = Time('2015-06-16 12:00:00')

    from astroplan.plots import plot_parallactic
    import matplotlib.pyplot as plt

    plot_parallactic(altair, subaru, time)
    plot_parallactic(vega, subaru, time)
    plot_parallactic(deneb, subaru, time)

    plt.legend(loc=2)
    plt.show()

We can also calculate this quantity directly:

.. code-block:: python

    >>> subaru.parallactic_angle(time, altair)
    −0.640582rad

    >>> subaru.parallactic_angle(time, vega)
    −0.465298rad

    >>> subaru.parallactic_angle(time, deneb)
    0.729871rad

The Moon
--------

If you need to take the Moon into account when observing, you may want to know
when it rises, sets, what phase its in, etc.

Let's first find out if the Moon is out during the time we defined earlier:

.. warning::

    *moon_rise_time* and *moon_set_time* have not yet been implemented.

.. code-block:: python

    >>> #subaru.moon_rise_time(time)

    >>> #subaru.moon_set_time(time)

We could also look at the Moon's alt/az coordinates:

.. code-block:: python

    >>> subaru.moon_altaz(time).alt
    −45∘05′18.2435′′

    >>> subaru.moon_altaz(time).az
    34∘35′57.5413′′

It looks like the Moon is well below the horizon at the time we picked before,
but we should check to see if it will be out during the window of time our
targets will be visible (again--as defined at the beginning of this tutorial):

.. code-block:: python

    >>> visible_time = start + (end - start)*np.linspace(0, 1, 20)

    >>> subaru.moon_altaz(visible_time).alt
    [−24∘15′08.8308′′ −29∘49′04.6286′′ −35∘03′43.449′′ −39∘53′16.0653′′
    −44∘09′59.8904′′ −47∘44′08.5089′′ −50∘24′19.9784′′ −51∘59′18.4053′′
    −52∘20′53.9214′′ −51∘27′04.0998′′ −49∘22′46.0578′′ −46∘17′54.7431′′
    −42∘24′06.7653′′ −37∘52′10.4174′′ −32∘50′59.3228′′ −27∘27′24.8625′′
    −21∘46′34.5241′′ −15∘52′15.6116′′ −9∘47′16.3944′′ −2∘11′39.571′′]

Looks like the Moon will be below the horizon during the entire time.

:ref:`Return to Top <summer_triangle_tutorial>`

.. _summer_triangle-sky_charts:

Sky Charts
==========

Now that we've determined the best times to observe our targets on the night in
question, let's take a look at the positions of our objects in the sky.

We can use `plot_sky` as a sanity check on our target's positions or even just
to better visualize our observation run.

Let's take the `start` and `end` of the time window we determined
:ref:`earlier <summer_triangle-observable>` (using the most basic definition
of "visible" targets, above the horizon when the sun is down), and see where our
targets lay in the sky::

    from astroplan.plots import plot_sky
    import matplotlib.pyplot as plt

    altair_style = {'color': 'r'}
    deneb_style = {'color': 'g'}

    plot_sky(altair, subaru, start, style_kwargs=altair_style)
    plot_sky(vega, subaru, start)
    plot_sky(deneb, subaru, start, style_kwargs=deneb_style)

    plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
    plt.show()

    plot_sky(altair, subaru, end, style_kwargs=altair_style)
    plot_sky(vega, subaru, end)
    plot_sky(deneb, subaru, end, style_kwargs=deneb_style)

    plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
    plt.show()

.. plot::

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

    from astropy.coordinates import SkyCoord
    from astroplan import FixedTarget

    coordinates = SkyCoord('19h50m47.6s', '+08d52m12.0s', frame='icrs')
    altair = FixedTarget(name='Altair', coord=coordinates)

    coordinates = SkyCoord('18h36m56.5s', '+38d47m06.6s', frame='icrs')
    vega = FixedTarget(name='Vega', coord=coordinates)

    coordinates = SkyCoord('20h41m25.9s', '+45d16m49.3s', frame='icrs')
    deneb = FixedTarget(name='Deneb', coord=coordinates)

    from astropy.time import Time

    # Here we need to add a second to our start time so that all objects show up.
    start = Time('2015-06-16 06:23:40.991') + 1 * u.second
    end = Time('2015-06-16 15:47:36.466')

    from astroplan.plots import plot_sky
    import matplotlib.pyplot as plt

    altair_style = {'color': 'r'}
    deneb_style = {'color': 'g'}

    plot_sky(altair, subaru, start, style_kwargs=altair_style)
    plot_sky(vega, subaru, start)
    plot_sky(deneb, subaru, start, style_kwargs=deneb_style)

    # Note that you don't need this code block to produce the plot.
    # It reduces the plot size for the documentation.
    ax = plt.gca()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.75, box.height * 0.75])

    plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
    plt.show()

    plot_sky(altair, subaru, end, style_kwargs=altair_style)
    plot_sky(vega, subaru, end)
    plot_sky(deneb, subaru, end, style_kwargs=deneb_style)

    # Note that you don't need this code block to produce the plot.
    # It reduces the plot size for the documentation.
    ax = plt.gca()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.75, box.height * 0.75])

    plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
    plt.show()

We can also show how our targets move over time during the night in question::

    time_window = start + (end - start) * np.linspace(0, 1, 10) * u.hour

    plot_sky(altair, subaru, time_window, style_kwargs=altair_style)
    plot_sky(vega, subaru, time_window)
    plot_sky(deneb, subaru, time_window, style_kwargs=deneb_style)

    plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
    plt.show()

.. plot::

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

    from astropy.coordinates import SkyCoord
    from astroplan import FixedTarget

    coordinates = SkyCoord('19h50m47.6s', '+08d52m12.0s', frame='icrs')
    altair = FixedTarget(name='Altair', coord=coordinates)

    coordinates = SkyCoord('18h36m56.5s', '+38d47m06.6s', frame='icrs')
    vega = FixedTarget(name='Vega', coord=coordinates)

    coordinates = SkyCoord('20h41m25.9s', '+45d16m49.3s', frame='icrs')
    deneb = FixedTarget(name='Deneb', coord=coordinates)

    from astropy.time import Time
    from astroplan.plots import plot_sky
    import matplotlib.pyplot as plt

    # Here we need to add a second to our start time so that all objects show up.
    start = Time('2015-06-16 06:23:40.991') + 1 * u.second
    end = Time('2015-06-16 15:47:36.466')

    time_window = start + (end - start) * np.linspace(0, 1, 10)

    altair_style = {'color': 'r'}
    deneb_style = {'color': 'g'}

    plot_sky(altair, subaru, time_window, style_kwargs=altair_style)
    plot_sky(vega, subaru, time_window)
    plot_sky(deneb, subaru, time_window, style_kwargs=deneb_style)

    plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))

    plt.tight_layout()

    plt.show()
