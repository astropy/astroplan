:orphan:

.. _plots:

*****
Plots
*****

`astroplan` currently has convenience functions for making three different types
of plots: airmass vs time, parallactic angle vs time and sky charts.  While
`astroplan` requires `Matplotlib`_, the use of additional plotting packages
(such as `Seaborn`) is neither explicitly prohibited or supported.

All `astroplan` plots return a `matplotlib.axes.Axes` object (often referred
to as `ax`).  You can further manipulate the returned `ax` object (including
using it as input for an `astroplan` plotting function, see examples below)
or you can simply display/print the plot.

Contents
========

    * :ref:`plots_time_dependent`

    * :ref:`plots_sky_charts`

.. warning::

    All examples here assume you know how to and have already constructed
    `astroplan.Observer` and `astroplan.FixedTarget` objects.


.. _plots_time_dependent:

Time Dependent Plots
====================

Although all `astroplan` plots are time-dependent in some way, we label those
that have a time-based axis as "time-dependent".

`Astroplan` currently has two different types of "time-dependent" plots:

* Airmass vs. Time
* Parallactic Angle vs. Time

These take, at minimum, `astroplan.Observer`, `astroplan.FixedTarget` and
`astropy.time.Time` objects as input.

.. _plots_airmass:

Airmass plots are made with a command such as::

    plot_airmass(target, observer, time)

.. _plots_parallactic:

Parallactic angle plots::

    plot_parallactic(target, observer, time)

Below are general guidelines for working with time-dependent plots in
`astroplan`.  Examples use airmass, but the following works for parallactic
angle as well.

.. seealso:

    ???? `astroplan.Observer.altaz.secz` ?????
    `astroplan.Observer.parallactic_angle`

Making a quick plot
-------------------

Any plot function in `astroplan` with a time-based axis will allow you to make
a quick plot over a 24-hour period

After constructing `Observer` and `Target` objects, construct a
`astropy.time.Time` object with a single instance in time and issue the
plotting command::

    observe_time = Time('2015-06-15 23:30:00')

    plot_airmass(target, observer, observe_time)
    plt.show()

.. plot::

    import matplotlib.pyplot as plt
    import astropy.units as u
    from astropy.coordinates import EarthLocation, SkyCoord
    from pytz import timezone
    from astropy.time import Time

    from astroplan import Observer
    from astroplan import FixedTarget
    from astroplan.plots import plot_airmass

    longitude = '-155d28m48.900s'
    latitude = '+19d49m42.600s'
    elevation = 4163 * u.m
    location = EarthLocation.from_geodetic(longitude, latitude, elevation)

    observer = Observer(name='Subaru Telescope',
                   location=location,
                   pressure=0.615 * u.bar,
                   relative_humidity=0.11,
                   temperature=0 * u.deg_C,
                   timezone=timezone('US/Hawaii'),
                   description="Subaru Telescope on Mauna Kea, Hawaii")

    coordinates = SkyCoord('06h45m08.9173s', '-16d42m58.017s', frame='icrs')
    target = FixedTarget(name='Sirius', coord=coordinates)

    observe_time = Time('2015-06-15 23:30:00')

    plot_airmass(target, observer, observe_time)
    plt.show()

As you can see, the 24-hour plot is centered on the time input.  You can use
array `astropy.time.Time` objects for these quick plots as well--they just
can't contain more than one instance in time::

    Time(['2015-06-15 23:30:00'])

    [Time('2015-06-15 23:30:00')]

.. plots_time_window::

Specifying a time window
------------------------

If you want to see airmass plotted over a window that is not 24 hours long or
you want to control the precision of the plot, you must specify every time for
which you want to see an airmass plotted.

To quickly populate an `astropy.time.Time` object with many instances of time,
use `Numpy`_ and `astropy.units`.

Centering the window at some time
+++++++++++++++++++++++++++++++++

To center your time window at some instance in time::

    import numpy as np

    observe_time = Time('2015-06-15 23:30:00')
    observe_time = observe_time + np.linspace(-5, 5, 55)*u.hour

    plot_airmass(target, observer, observe_time)
    plt.show()

.. plot::

    import matplotlib.pyplot as plt
    import astropy.units as u
    from astropy.coordinates import EarthLocation, SkyCoord
    from pytz import timezone
    from astropy.time import Time

    from astroplan import Observer
    from astroplan import FixedTarget
    from astroplan.plots import plot_airmass

    # Set up Observer, Target and observation time objects.
    longitude = '-155d28m48.900s'
    latitude = '+19d49m42.600s'
    elevation = 4163 * u.m
    location = EarthLocation.from_geodetic(longitude, latitude, elevation)

    observer = Observer(name='Subaru Telescope',
                   location=location,
                   pressure=0.615 * u.bar,
                   relative_humidity=0.11,
                   temperature=0 * u.deg_C,
                   timezone=timezone('US/Hawaii'),
                   description="Subaru Telescope on Mauna Kea, Hawaii")

    coordinates = SkyCoord('06h45m08.9173s', '-16d42m58.017s', frame='icrs')
    target = FixedTarget(name='Sirius', coord=coordinates)

    import numpy as np

    observe_time = Time('2015-06-15 23:30:00')
    observe_time = observe_time + np.linspace(-5, 5, 55)*u.hour

    plot_airmass(target, observer, observe_time)
    plt.show()

Specify start and end times
+++++++++++++++++++++++++++

If you know the start and end times of your observation run, you can use a
`astropy.time.DeltaTime` object to create an array for time input::

    start_time = Time('2015-06-15 20:00:00')
    end_time = Time('2015-06-16 04:00:00')
    delta_t = end_time - start_time
    observe_time = start_time + delta_t*np.linspace(0, 1, 75)

    plot_airmass(target, observer, observe_time)
    plt.show()

.. plot::

    import matplotlib.pyplot as plt
    import astropy.units as u
    from astropy.coordinates import EarthLocation, SkyCoord
    from pytz import timezone
    from astropy.time import Time

    from astroplan import Observer
    from astroplan import FixedTarget
    from astroplan.plots import plot_airmass

    # Set up Observer, Target and observation time objects.
    longitude = '-155d28m48.900s'
    latitude = '+19d49m42.600s'
    elevation = 4163 * u.m
    location = EarthLocation.from_geodetic(longitude, latitude, elevation)

    observer = Observer(name='Subaru Telescope',
                   location=location,
                   pressure=0.615 * u.bar,
                   relative_humidity=0.11,
                   temperature=0 * u.deg_C,
                   timezone=timezone('US/Hawaii'),
                   description="Subaru Telescope on Mauna Kea, Hawaii")

    coordinates = SkyCoord('06h45m08.9173s', '-16d42m58.017s', frame='icrs')
    target = FixedTarget(name='Sirius', coord=coordinates)

    start_time = Time('2015-06-15 20:00:00')
    end_time = Time('2015-06-16 04:00:00')
    delta_t = end_time - start_time
    observe_time = start_time + delta_t*np.linspace(0, 1, 75)

    plot_airmass(target, observer, observe_time)
    plt.show()

Plotting a quantity for multiple targets
----------------------------------------

If you want to plot airmass information for multiple targets, simply reissue
the ``plot_airmass`` command, using a different `Target` object as input this
time. Repeat until you have as many targets on the plot as you wish.

When you're ready to make a different plot, use *ax.cla()* to clear the `ax`
object::

    coordinates = SkyCoord('02h31m49.09s', '+89d15m50.8s', frame='icrs')
    other_target = FixedTarget(name='Polaris', coord=coordinates)

    coordinates = SkyCoord('07h45m19.4s', '+28d01m35s', frame='icrs')
    third_target = FixedTarget(name='Pollux', coord=coordinates)

    observe_time = Time('2015-06-30 23:30:00') + np.linspace(-7.0, 5.5, 50)*u.hour

    plot_airmass(target, observer, observe_time)
    plot_airmass(other_target, observer, observe_time)
    plot_airmass(third_target, observer, observe_time)

    plt.legend(shadow=True, loc=2)
    plt.show()

.. plot::

    import matplotlib.pyplot as plt
    import astropy.units as u
    from astropy.coordinates import EarthLocation, SkyCoord
    from pytz import timezone
    from astropy.time import Time

    from astroplan import Observer
    from astroplan import FixedTarget
    from astroplan.plots import plot_airmass

    longitude = '-155d28m48.900s'
    latitude = '+19d49m42.600s'
    elevation = 4163 * u.m
    location = EarthLocation.from_geodetic(longitude, latitude, elevation)

    observer = Observer(name='Subaru Telescope',
                   location=location,
                   pressure=0.615 * u.bar,
                   relative_humidity=0.11,
                   temperature=0 * u.deg_C,
                   timezone=timezone('US/Hawaii'),
                   description="Subaru Telescope on Mauna Kea, Hawaii")

    coordinates = SkyCoord('06h45m08.9173s', '-16d42m58.017s', frame='icrs')
    target = FixedTarget(name='Sirius', coord=coordinates)

    coordinates = SkyCoord('02h31m49.09s', '+89d15m50.8s', frame='icrs')
    other_target = FixedTarget(name='Polaris', coord=coordinates)

    coordinates = SkyCoord('07h45m19.4s', '+28d01m35s', frame='icrs')
    third_target = FixedTarget(name='Pollux', coord=coordinates)

    observe_time = Time('2015-06-30 23:30:00') + np.linspace(-7.0, 5.5, 50)*u.hour

    plot_airmass(target, observer, observe_time)
    plot_airmass(other_target, observer, observe_time)
    plot_airmass(third_target, observer, observe_time)

    plt.legend(shadow=True, loc=2)
    plt.show()

Changing style options
----------------------

You can set the `Matplotlib` *linestyle* and *color* options by passing in a
dictionary with your preferences::

    sirius_styles = {'linestyle': '--', 'color': 'r'}
    polaris_styles = {'linestyle': '-', 'color': 'g'}

    plot_airmass(other_target, observer, observe_time, style_kwargs=sirius_styles)
    plot_airmass(third_target, observer, observe_time, style_kwargs=pollux_styles)

    plt.legend(shadow=True, loc=2)
    plt.show()

.. plot::

    import matplotlib.pyplot as plt
    import astropy.units as u
    from astropy.coordinates import EarthLocation, SkyCoord
    from pytz import timezone
    from astropy.time import Time

    from astroplan import Observer
    from astroplan import FixedTarget
    from astroplan.plots import plot_airmass

    longitude = '-155d28m48.900s'
    latitude = '+19d49m42.600s'
    elevation = 4163 * u.m
    location = EarthLocation.from_geodetic(longitude, latitude, elevation)

    observer = Observer(name='Subaru Telescope',
                   location=location,
                   pressure=0.615 * u.bar,
                   relative_humidity=0.11,
                   temperature=0 * u.deg_C,
                   timezone=timezone('US/Hawaii'),
                   description="Subaru Telescope on Mauna Kea, Hawaii")

    coordinates = SkyCoord('06h45m08.9173s', '-16d42m58.017s', frame='icrs')
    target = FixedTarget(name='Sirius', coord=coordinates)

    coordinates = SkyCoord('02h31m49.09s', '+89d15m50.8s', frame='icrs')
    other_target = FixedTarget(name='Polaris', coord=coordinates)

    coordinates = SkyCoord('07h45m19.4s', '+28d01m35s', frame='icrs')
    third_target = FixedTarget(name='Pollux', coord=coordinates)

    observe_time = Time('2015-06-30 23:30:00') + np.linspace(-7.0, 5.5, 50)*u.hour

    sirius_styles = {'linestyle': '--', 'color': 'r'}
    pollux_styles = {'linestyle': '-', 'color': 'g'}

    plot_airmass(other_target, observer, observe_time, style_kwargs=sirius_styles)
    plot_airmass(third_target, observer, observe_time, style_kwargs=pollux_styles)

    plt.legend(shadow=True, loc=2)
    plt.show()

:ref:`Return to Top <plots>`


.. _plots_sky_charts:

Sky Charts
==========

Many users planning an observation run will want to see the positions of targets
with respect to their local horizon, as well as the positions of familiar stars
or other objects to act as guides.

.. _plots_sky_charts::

`plot_sky` allows you to plot the positions of targets at a single instance or
over some window of time with a command such as::

    plot_sky(target, observer, time)

.. warning::

    Note that the time input for `plot_sky` has to either be an array of
    `astropy.time.Time` objects or has to be an `astropy.time.Time` object
    containing an array of times--in other words, it **cannot** be scalar.

.. warning::

    `plot_sky` currently produces polar plots in altitude/azimuth coordinates
    only.  Plots are centered on the observer's zenith.

.. seealso:

    `astroplan.Observer.altaz`

Making a plot for one instance in time
--------------------------------------

After constructing your `Observer` and `FixedTarget` objects, construct a time
input using an array of length 1.

That is, either an `astropy.time.Time` object with an array containing one time
value (e.g., ``Time(['2000-1-1'])``) or an array containing one scalar
`astropy.time.Time` object (e.g., ``[Time('2000-1-1')]``).

Let's say that you created `FixedTarget` objects for Polaris, Altair, Vega and
Deneb.  To plot a map of the sky::

    observe_time = Time(['2015-03-15 15:30:00'])

    polaris_style = {color': 'k'}
    altair_style = {'color': 'b'}
    vega_style = {'color': 'g'}
    deneb_style = {'color': 'r'}

    plot_sky(polaris, observer, observe_time, style_kwargs=polaris_style)
    plot_sky(altair, observer, observe_time, style_kwargs=altair_style)
    plot_sky(vega, observer, observe_time, style_kwargs=vega_style)
    plot_sky(deneb, observer, observe_time, style_kwargs=deneb_style)

    plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
    plt.show()

.. warning::

    Since `plot_sky` uses `matplotlib.pyplot.scatter` (which gives the same
    color to different plots made on the same set of axes), you have to specify
    the color for each target via a style dictionary if you don't want all
    targets to have the same color.

.. plot::

    import matplotlib.pyplot as plt
    import astropy.units as u
    from astropy.coordinates import EarthLocation, SkyCoord
    from pytz import timezone
    from astropy.time import Time

    from astroplan import Observer
    from astroplan import FixedTarget
    from astroplan.plots import plot_sky

    # Set up Observer, Target and observation time objects.
    longitude = '-155d28m48.900s'
    latitude = '+19d49m42.600s'
    elevation = 4163 * u.m
    location = EarthLocation.from_geodetic(longitude, latitude, elevation)

    observer = Observer(name='Subaru Telescope',
                   location=location,
                   pressure=0.615 * u.bar,
                   relative_humidity=0.11,
                   temperature=0 * u.deg_C,
                   timezone=timezone('US/Hawaii'),
                   description="Subaru Telescope on Mauna Kea, Hawaii")

    coordinates = SkyCoord('02h31m49.09s', '+89d15m50.8s', frame='icrs')
    polaris = FixedTarget(name='Polaris', coord=coordinates)
    polaris_style = {'marker': 'o', 'color': 'k'}

    coordinates = SkyCoord('19h50m47.6s', '+08d52m12.0s', frame='icrs')
    altair = FixedTarget(name='Altair', coord=coordinates)
    altair_style = {'marker': 'o', 'color': 'b'}

    coordinates = SkyCoord('18h36m56.5s', '+38d47m06.6s', frame='icrs')
    vega = FixedTarget(name='Vega', coord=coordinates)
    vega_style = {'marker': 'o', 'color': 'g'}

    coordinates = SkyCoord('20h41m25.9s', '+45d16m49.3s', frame='icrs')
    deneb = FixedTarget(name='Deneb', coord=coordinates)
    deneb_style = {'marker': 'o', 'color': 'r'}

    # Note that this is not a scalar.
    observe_time = Time(['2015-03-15 15:30:00'])

    plot_sky(polaris, observer, observe_time, style_kwargs=polaris_style)
    plot_sky(altair, observer, observe_time, style_kwargs=altair_style)
    plot_sky(vega, observer, observe_time, style_kwargs=vega_style)
    plot_sky(deneb, observer, observe_time, style_kwargs=deneb_style)

    plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
    plt.show()

Showing movement over time
--------------------------

If you want to see how your targets move over time, you need to explicitly
specify every instance in time.

Say I want to know how Altair moves in the sky over a 9-hour period::

    import numpy as np

    observe_time = Time('2015-03-15 17:00:00')
    observe_time = observe_time + np.linspace(-4, 5, 10)*u.hour

    plot_sky(altair, observer, observe_time)

    plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
    plt.show()

.. plot::

    import matplotlib.pyplot as plt
    import astropy.units as u
    from astropy.coordinates import EarthLocation, SkyCoord
    from pytz import timezone
    from astropy.time import Time

    from astroplan import Observer
    from astroplan import FixedTarget
    from astroplan.plots import plot_sky

    # Set up Observer, Target and observation time objects.
    longitude = '-155d28m48.900s'
    latitude = '+19d49m42.600s'
    elevation = 4163 * u.m
    location = EarthLocation.from_geodetic(longitude, latitude, elevation)

    observer = Observer(name='Subaru Telescope',
                   location=location,
                   pressure=0.615 * u.bar,
                   relative_humidity=0.11,
                   temperature=0 * u.deg_C,
                   timezone=timezone('US/Hawaii'),
                   description="Subaru Telescope on Mauna Kea, Hawaii")

    coordinates = SkyCoord('19h50m47.6s', '+08d52m12.0s', frame='icrs')
    altair = FixedTarget(name='Altair', coord=coordinates)

    import numpy as np

    observe_time = Time('2015-03-15 17:00:00')
    observe_time = observe_time + np.linspace(-4, 5, 10)*u.hour

    plot_sky(altair, observer, observe_time)

    plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
    plt.show()

For more examples on how to populate time objects, see `astropy.time.Time`
documentation, or :ref:`plots_time_window`.

.. warning::

    Note that in the case of an object being under the horizon (or having
    negative altitude) at any of the times in your **time** input, `plot_sky`
    will warn you.  Your object(s) will not show up on the plot for those
    particular times, but any positions above the horizon will still be plotted
    as normal.

Customizing your sky plot
-------------------------

Style options
+++++++++++++

Changing coordinate defaults
++++++++++++++++++++++++++++

Other tweaks
++++++++++++

:ref:`Return to Top <plots>`

Explicitly passing in a `matplotlib.axes.Axes` object
=====================================================
