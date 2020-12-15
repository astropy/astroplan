.. _summer_triangle_tutorial:

.. todo::

    Add section on moon phases, illumination fraction, etc.

    Replace target construction with easier site name function.

    Update with constraints when available?

    Show users how to create/use Time objects in their own timezone.

*****************************
Observing the Summer Triangle
*****************************

.. note::

    Your calculated rise/set and other times may differ slightly from those in
    this tutorial, on the order of ~1 second.  This is a normal variance in
    precision due to several factors, including varying :ref:`IERS tables <iers>`
    and machine architecture.

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

First, we define our `~astroplan.Observer` object:

.. code-block:: python

    >>> from astroplan import Observer

    >>> subaru = Observer.at_site('subaru')

Then, we define our `~astroplan.FixedTarget`'s, since the Summer Triangle is
fixed with respect to the celestial sphere (if we ignore the relatively small
proper motion). We will use the `~astroplan.FixedTarget.from_name` class method,
which queries the CDS name resolver for your target's coordinates (giving you
the power of SIMBAD!):

.. code-block:: python

    >>> from astropy.coordinates import SkyCoord
    >>> from astroplan import FixedTarget

    >>> altair = FixedTarget.from_name('Altair')
    >>> vega = FixedTarget.from_name('Vega')


For objects that can't be resolved with `~astroplan.FixedTarget.from_name`, you
can enter coordinates manually:

.. code-block:: python

    >>> coordinates = SkyCoord('20h41m25.9s', '+45d16m49.3s', frame='icrs')
    >>> deneb = FixedTarget(name='Deneb', coord=coordinates)

We also have to define a `~astropy.time.Time` (in UTC) at which we wish to
observe.  Here, we pick 2AM local time, which is noon UTC during the
summer::

    >>> from astropy.time import Time

    >>> time = Time('2015-06-16 12:00:00')

:ref:`Return to Top <summer_triangle_tutorial>`

.. _summer_triangle-observable:

Observable?
===========

Next, it would be handy to know if our targets are visible from Subaru at the
time we settled on.  In other words--are they above the horizon while the Sun
is down?

.. code-block:: python

    >>> subaru.target_is_up(time, altair)
    True

    >>> subaru.target_is_up(time, vega)
    True

    >>> subaru.target_is_up(time, deneb)
    True

...They are!

What if we weren't sure if the Sun is down at this time:

.. code-block:: python

    >>> subaru.is_night(time)
    True

...It is!

However, we may want to find a window of time for tonight during which all
three of our targets are above the horizon *and* the Sun is below the horizon
(let's worry about light pollution from the Moon later).

Let's define the window of time during which all targets are above the horizon.

Note that because of the precision limitations of rise/set calculations
(altitudes at these times won't equal precisely zero, but will be off by a few
arc seconds), we'll manually adjust rise/set times by a few minutes.

.. code-block:: python

    >>> import numpy as np
    >>> import astropy.units as u

    >>> altair_rise = subaru.target_rise_time(time, altair) + 5*u.minute
    >>> altair_set = subaru.target_set_time(time, altair) - 5*u.minute

    >>> vega_rise = subaru.target_rise_time(time, vega) + 5*u.minute
    >>> vega_set = subaru.target_set_time(time, vega) - 5*u.minute

    >>> deneb_rise = subaru.target_rise_time(time, deneb) + 5*u.minute
    >>> deneb_set = subaru.target_set_time(time, deneb) - 5*u.minute

    >>> all_up_start = np.max([altair_rise, vega_rise, deneb_rise])
    >>> all_up_end = np.min([altair_set, vega_set, deneb_set])

Now, let's find sunset and sunrise for tonight (and confirm that they are
indeed those for tonight):

.. code-block:: python

    >>> sunset_tonight = subaru.sun_set_time(time, which='nearest')

    >>> sunset_tonight.iso # doctest: +SKIP
    '2015-06-16 04:59:11.267'

This is '2015-06-15 18:59:11.267' in the Hawaii time zone (that's where Subaru
is).

.. code-block:: python

    >>> sunrise_tonight = subaru.sun_rise_time(time, which='nearest')

    >>> sunrise_tonight.iso # doctest: +SKIP
    '2015-06-16 15:47:35.822'

This is '2015-06-16 05:47:35.822' Hawaii time.

Sunset and sunrise check out, so now we define the limits of our observation
window:

.. code-block:: python

    >>> start = np.max([sunset_tonight, all_up_start])
    >>> start.iso # doctest: +SKIP
    '2015-06-16 06:28:40.126'

    >>> end = np.min([sunrise_tonight, all_up_end])
    >>> end.iso # doctest: +SKIP
    '2015-06-16 15:47:35.822'

So, our targets will be visible (as we've defined it above) from
'2015-06-15 20:28:40.126' to '2015-06-16 05:47:35.822' Hawaii time.  Depending
on our observation goals, this window of time may be good enough for preliminary
planning, or we may want to optimize our observational conditions.  If the
latter is the case, go on to the Optimal Observation Time section (immediately
below).

:ref:`Return to Top <summer_triangle_tutorial>`

.. _summer_triangle-optimal_observation:

Optimal Observation Time
========================

There are a few things we can look at to find the best time to observe our
targets on a given night.

Airmass
-------

To get a general idea of our targets' airmass on the night of observation, we
can plot it over the course of the night (for more on plotting see :doc:`plots`):

.. code-block:: python

    >>> from astroplan.plots import plot_airmass # doctest: +SKIP
    >>> import matplotlib.pyplot as plt # doctest: +SKIP

    >>> plot_airmass(altair, subaru, time) # doctest: +SKIP
    >>> plot_airmass(vega, subaru, time) # doctest: +SKIP
    >>> plot_airmass(deneb, subaru, time)  # doctest: +SKIP

    >>> plt.legend(loc=1, bbox_to_anchor=(1, 1)) # doctest: +SKIP
    >>> plt.show() # doctest: +SKIP

.. plot::

    import astropy.units as u
    from astropy.coordinates import EarthLocation, SkyCoord
    import matplotlib.pyplot as plt
    from pytz import timezone

    from astroplan import Observer, FixedTarget
    from astroplan.plots import plot_airmass

    longitude = '-155d28m48.900s'
    latitude = '+19d49m42.600s'
    elevation = 4163 * u.m
    location = EarthLocation.from_geodetic(longitude, latitude, elevation)

    subaru = Observer(name='Subaru Telescope',
                      location=location,
                      timezone=timezone('US/Hawaii'),
                      description="Subaru Telescope on Maunakea, Hawaii")

    coordinates = SkyCoord('19h50m47.6s', '+08d52m12.0s', frame='icrs')
    altair = FixedTarget(name='Altair', coord=coordinates)

    coordinates = SkyCoord('18h36m56.5s', '+38d47m06.6s', frame='icrs')
    vega = FixedTarget(name='Vega', coord=coordinates)

    coordinates = SkyCoord('20h41m25.9s', '+45d16m49.3s', frame='icrs')
    deneb = FixedTarget(name='Deneb', coord=coordinates)

    from astropy.time import Time

    time = Time('2015-06-16 12:00:00')

    plot_airmass(altair, subaru, time)
    plot_airmass(vega, subaru, time)
    plot_airmass(deneb, subaru, time)

    # Note that you don't need this code block to produce the plot.
    # It reduces the plot size for the documentation.
    ax = plt.gca()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height * 0.8])

    plt.legend(loc=1, bbox_to_anchor=(1.35, 1))
    plt.tight_layout()
    plt.show()

We want a minimum airmass when observing, and it looks like sometime between
9:00 and 15:00 UTC (or 23:00 on the 15th to 5:00 on the 16th, US/Hawaii) would
be the best time to observe all three targets.

However, if we want to define a more specific time window based on airmass, we
can calculate this quantity directly. To get airmass measurements, we need to
use the ``AltAz`` frame:

.. code-block:: python

    >>> subaru.altaz(time, altair).secz # doctest: +SKIP
    <Quantity 1.0302347952130682>

    >>> subaru.altaz(time, vega).secz # doctest: +SKIP
    <Quantity 1.0690421636016616>

    >>> subaru.altaz(time, deneb).secz # doctest: +SKIP
    <Quantity 1.167753811648361>

Behind the scenes here, ``subaru.altaz(time, altair)`` is actually creating
an `~astropy.coordinates.AltAz` object in the ``AltAz`` frame, so if you
know how to work with `~astropy.coordinates` objects, you can do lots more
than just computing airmass.

Parallactic Angle
-----------------

To get a general idea of our targets' parallactic angle on the night of
observation, we can make another plot (again, see :doc:`plots` for more on
customizing plots and the like):

.. code-block:: python

    >>> import matplotlib.pyplot as plt # doctest: +SKIP
    >>> from astroplan.plots import plot_parallactic # doctest: +SKIP

    >>> plot_parallactic(altair, subaru, time) # doctest: +SKIP
    >>> plot_parallactic(vega, subaru, time) # doctest: +SKIP
    >>> plot_parallactic(deneb, subaru, time) # doctest: +SKIP

    >>> plt.legend(loc=2) # doctest: +SKIP
    >>> plt.show() # doctest: +SKIP

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
                   description="Subaru Telescope on Maunakea, Hawaii")

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
    plt.tight_layout()
    plt.show()

We can also calculate the parallactic angle directly:

.. code-block:: python

    >>> subaru.parallactic_angle(time, altair) # doctest: +SKIP
    <Angle -0.6404957821112053 rad>

    >>> subaru.parallactic_angle(time, vega) # doctest: +SKIP
    <Angle -0.46542183982024 rad>

    >>> subaru.parallactic_angle(time, deneb) # doctest: +SKIP
    <Angle 0.7297067855978494 rad>

The `~astropy.coordinates.Angle` objects resulting from the calls to
``parallactic_angle()`` are subclasses of the `~astropy.units.Quantity`
class, so they can do everything a `~astropy.units.Quantity` can do -
basically they work like numbers with attached units, and keep track of
units so you don't have to.

For more on the many things you can do with these, take a look at the
`Astropy <https://astropy.org>`__ documentation or tutorials.  For now the most useful thing is to
know is that ``angle.degree``, ``angle.hourangle``, and ``angle.radian``
give you back Python floats (or `numpy` arrays) for the angle in degrees,
hours, or radians.

The Moon
--------

If you need to take the Moon into account when observing, you may want to know
when it rises, sets, what phase it's in, etc. Let's first find out if the Moon
is out during the time we defined earlier:

.. code-block:: python

    >>> subaru.moon_rise_time(time) # doctest: +SKIP
    <Time object: scale='utc' format='jd' value=2457190.1696768994>

    >>> subaru.moon_set_time(time) # doctest: +SKIP
    <Time object: scale='utc' format='jd' value=2457189.684134357>

We could also look at the Moon's alt/az coordinates:

.. code-block:: python

    >>> subaru.moon_altaz(time).alt # doctest: +SKIP
    <Latitude -45.08860929634166 deg>

    >>> subaru.moon_altaz(time).az # doctest: +SKIP
    <Longitude 34.605498354422686 deg>

It looks like the Moon is well below the horizon at the time we picked before,
but we should check to see if it will be out during the window of time our
targets will be visible (again--as defined at the beginning of this tutorial):

.. code-block:: python

    >>> visible_time = start + (end - start)*np.linspace(0, 1, 20)

    >>> subaru.moon_altaz(visible_time).alt # doctest: +SKIP
    <Latitude [-25.21127325,-30.68088873,-35.82145644,-40.53415037,
               -44.68898859,-48.12296182,-50.64971858,-52.08946099,
               -52.31849772,-51.31548444,-49.17038499,-46.04862654,
               -42.13887599,-37.61479774,-32.61875342,-27.26048709,
               -21.62215227,-15.76463668, -9.73313141, -2.19408792] deg>


Looks like the Moon will be below the horizon during the entire time.

:ref:`Return to Top <summer_triangle_tutorial>`

.. _summer_triangle-sky_charts:

Sky Charts
==========

Now that we've determined the best times to observe our targets on the night in
question, let's take a look at the positions of our objects in the sky.

We can use `~astroplan.plots.plot_sky` as a sanity check on our target's
positions or even just to better visualize our observation run.

Let's take the ``start`` and ``end`` of the time window we determined
:ref:`earlier <summer_triangle-observable>` (using the most basic definition
of "visible" targets, above the horizon when the sun is down), and see where our
targets lay in the sky:

.. code-block:: python

    >>> from astroplan.plots import plot_sky
    >>> import matplotlib.pyplot as plt  # doctest: +SKIP

    >>> altair_style = {'color': 'r'}
    >>> deneb_style = {'color': 'g'}

    >>> plot_sky(altair, subaru, start, style_kwargs=altair_style)  # doctest: +SKIP
    >>> plot_sky(vega, subaru, start)  # doctest: +SKIP
    >>> plot_sky(deneb, subaru, start, style_kwargs=deneb_style)  # doctest: +SKIP

    >>> plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))  # doctest: +SKIP
    >>> plt.show()  # doctest: +SKIP

    >>> plot_sky(altair, subaru, end, style_kwargs=altair_style)  # doctest: +SKIP
    >>> plot_sky(vega, subaru, end)  # doctest: +SKIP
    >>> plot_sky(deneb, subaru, end, style_kwargs=deneb_style)  # doctest: +SKIP

    >>> plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))  # doctest: +SKIP
    >>> plt.show()  # doctest: +SKIP

.. plot::


    from astropy.coordinates import EarthLocation, SkyCoord
    from astropy.time import Time
    import astropy.units as u
    import matplotlib.pyplot as plt
    from pytz import timezone

    from astroplan import Observer, FixedTarget
    from astroplan.plots import plot_sky

    longitude = '-155d28m48.900s'
    latitude = '+19d49m42.600s'
    elevation = 4163 * u.m
    location = EarthLocation.from_geodetic(longitude, latitude, elevation)

    subaru = Observer(name='Subaru Telescope',
                      location=location,
                      timezone=timezone('US/Hawaii'),
                      description="Subaru Telescope on Maunakea, Hawaii")

    coordinates = SkyCoord('19h50m47.6s', '+08d52m12.0s', frame='icrs')
    altair = FixedTarget(name='Altair', coord=coordinates)

    coordinates = SkyCoord('18h36m56.5s', '+38d47m06.6s', frame='icrs')
    vega = FixedTarget(name='Vega', coord=coordinates)

    coordinates = SkyCoord('20h41m25.9s', '+45d16m49.3s', frame='icrs')
    deneb = FixedTarget(name='Deneb', coord=coordinates)

    start = Time('2015-06-16 06:28:40.126')
    end = Time('2015-06-16 15:47:35.822')

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
    plt.tight_layout()
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
    plt.tight_layout()
    plt.show()

We can also show how our targets move over time during the night in question::

    >>> time_window = start + (end - start) * np.linspace(0, 1, 10)

    >>> plot_sky(altair, subaru, time_window, style_kwargs=altair_style)  # doctest: +SKIP
    >>> plot_sky(vega, subaru, time_window)  # doctest: +SKIP
    >>> plot_sky(deneb, subaru, time_window, style_kwargs=deneb_style)  # doctest: +SKIP

    >>> plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))  # doctest: +SKIP
    >>> plt.show()  # doctest: +SKIP

.. plot::

    import numpy as np
    from astropy.coordinates import EarthLocation, SkyCoord
    from astropy.time import Time
    import astropy.units as u
    import matplotlib.pyplot as plt
    from pytz import timezone

    from astroplan import Observer, FixedTarget
    from astroplan.plots import plot_sky

    longitude = '-155d28m48.900s'
    latitude = '+19d49m42.600s'
    elevation = 4163 * u.m
    location = EarthLocation.from_geodetic(longitude, latitude, elevation)

    subaru = Observer(name='Subaru Telescope',
                      location=location,
                      timezone=timezone('US/Hawaii'),
                      description="Subaru Telescope on Maunakea, Hawaii")

    coordinates = SkyCoord('19h50m47.6s', '+08d52m12.0s', frame='icrs')
    altair = FixedTarget(name='Altair', coord=coordinates)

    coordinates = SkyCoord('18h36m56.5s', '+38d47m06.6s', frame='icrs')
    vega = FixedTarget(name='Vega', coord=coordinates)

    coordinates = SkyCoord('20h41m25.9s', '+45d16m49.3s', frame='icrs')
    deneb = FixedTarget(name='Deneb', coord=coordinates)

    start = Time('2015-06-16 06:28:40.126')
    end = Time('2015-06-16 15:47:35.822')

    time_window = start + (end - start) * np.linspace(0, 1, 10)

    altair_style = {'color': 'r'}
    deneb_style = {'color': 'g'}

    plot_sky(altair, subaru, time_window, style_kwargs=altair_style)
    plot_sky(vega, subaru, time_window)
    plot_sky(deneb, subaru, time_window, style_kwargs=deneb_style)

    plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))

    plt.tight_layout()

    plt.show()
