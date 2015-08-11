:orphan:

.. _plots:

*****
Plots
*****

`Astroplan` currently has convenience functions for making three different types
of plots: airmass vs time, parallactic angle vs time and sky charts.  While
`Astroplan` requires `Matplotlib`, the use of other plotting packages
(such as `Seaborn`) is neither explicitly prohibited or supported.

.. warning::

    All examples here assume you know how to and have already constructed
    :ref:`observer` and :ref:`FixedTarget <targets_fixed_target_object>`
    objects.


.. _plots_airmass:

Airmass vs Time
===============

Airmass vs time plots are made with a call to the `plot_airmass` function.
It takes, at minimum, a `Target`, and `Observer` and a `Time` object as input.
Optional arguments include an `Axes` object and a plotting style dictionary.

plot_airmass will return an Axes object that contains airmass data for the
window of time specified by the Time object passed in. You can further
manipulate the returned Axes object or print the plot to your display/a file.

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
arrays for these quick plots as well--they just can't contain more than one
instance in time::

    Time(['2015-06-15 23:30:00'])

    [Time('2015-06-15 23:30:00')]

Specifying a time window
------------------------

If you want to see airmass plotted over a window that is not 24 hours long or
you want to control the precision of the plot, you must specify every time for
which you want to see an airmass plotted.

To quickly populate a ``Time`` object with many instances of time, use `Numpy`_
and ``astropy.units``.

Centering the window at some ``Time``
+++++++++++++++++++++++++++++++++++++

To center your time window at some ``Time``::

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
DeltaTime object::

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

Plotting airmass for multiple targets
-------------------------------------

If you want to plot airmass information for multiple targets, simply reissue
the ``plot_airmass command``, using a different Target object as input this
time. Repeat until you have as many targets on the plot as you wish.

When you're ready to make a different plot, use *ax.cla()* to clear the ``Axes``
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

You can set the *linestyle* and *color* options by passing in a style
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

.. _plots_parallactic_angle:

Parallactic Angle vs Time
=========================



.. _plots_sky_chart:

Sky Chart
=========
