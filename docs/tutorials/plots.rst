.. include:: ../references.txt

.. _plots:

.. doctest-skip-all

***********************
Plotting with Astroplan
***********************

`astroplan` currently has convenience functions for making three different types
of plots: airmass vs time, parallactic angle vs time and sky charts. This
plotting functionality in `astroplan` requires `Matplotlib`_ (although non-
plotting functionality will work even without `Matplotlib`_ ).  The use of
additional plotting packages (like `Seaborn
<http://stanford.edu/~mwaskom/software/seaborn/>`_) is not explicitly prevented,
but may or may not actually work.

All `astroplan` plots return a `~matplotlib.axes.Axes` object (which by
convention is assigned to the name ``ax`` in these tutorials).  You can further
manipulate the returned ``ax`` object, including using it as input for more
`astroplan` plotting functions, or you can simply display/print the plot.

Contents
========

    * :ref:`plots_time_dependent`

    * :ref:`plots_sky_charts`

    * :ref:`finder_image`


.. _plots_time_dependent:

Time Dependent Plots
====================

Although all `astroplan` plots are time-dependent in some way, we label those
that have a time-based axis as "time-dependent".

`astroplan` currently has a few different types of "time-dependent" plots,
for example `~astroplan.plots.plot_airmass`, `~astroplan.plots.plot_altitude`
and `~astroplan.plots.plot_parallactic`. These take, at minimum,
`~astroplan.Observer`, `~astroplan.FixedTarget` and `~astropy.time.Time` objects
as input.

.. _plots_airmass:

Airmass vs time plots are made the following way:

.. code-block:: python

    >>> from astroplan.plots import plot_airmass

    >>> plot_airmass(target, observer, time)

.. _plots_parallactic:

Parallactic angle vs time plots are made the following way:

.. code-block:: python

    >>> from astroplan.plots import plot_airmass

    >>> plot_parallactic(target, observer, time)

Below are general guidelines for working with time-dependent plots in
`astroplan`.  Examples use airmass but apply to parallactic angle as well.

.. seealso::

    `astropy.coordinates.AltAz.secz`

    `astroplan.Observer.parallactic_angle`

Making a quick plot
-------------------

Any plot function in `astroplan` with a time-based axis will allow you to make
a quick plot over a 24-hour period.

After constructing `~astroplan.Observer` and `~astroplan.FixedTarget`
objects, construct a `~astropy.time.Time` object with a single instance in
time and issue the plotting command.

.. code-block:: python

    >>> import matplotlib.pyplot as plt
    >>> from astropy.time import Time
    >>> from astroplan.plots import plot_airmass

    >>> observe_time = Time('2000-06-15 23:30:00')

    >>> plot_airmass(target, observer, observe_time)
    >>> plt.show()

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
                   description="Subaru Telescope on Maunakea, Hawaii")

    coordinates = SkyCoord('06h45m08.9173s', '-16d42m58.017s', frame='icrs')
    target = FixedTarget(name='Sirius', coord=coordinates)

    observe_time = Time('2000-06-15 23:30:00')

    plot_airmass(target, observer, observe_time)
    plt.tight_layout()
    plt.show()

As you can see, the 24-hour plot is centered on the *time* input.  You can also
use array `~astropy.time.Time` objects for these quick plots--they just
can't contain more than one instance in time.

For example, these are acceptable *time* inputs::

    Time(['2000-06-15 23:30:00'])

    [Time('2000-06-15 23:30:00')]

You can also add a second y-axis (on the right side) which shows the corresponding
altitude of the targets:

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
                   description="Subaru Telescope on Maunakea, Hawaii")

    coordinates = SkyCoord('06h45m08.9173s', '-16d42m58.017s', frame='icrs')
    target = FixedTarget(name='Sirius', coord=coordinates)

    observe_time = Time('2000-06-15 23:30:00')

    plot_airmass(target, observer, observe_time, altitude_yaxis=True)
    plt.tight_layout()
    plt.show()

You can make altitude the primary y-axis rather than airmass by using
`~astroplan.plots.plot_altitude`:

.. plot::

    import matplotlib.pyplot as plt
    import astropy.units as u
    from astropy.coordinates import EarthLocation, SkyCoord
    from pytz import timezone
    from astropy.time import Time

    from astroplan import Observer
    from astroplan import FixedTarget
    from astroplan.plots import plot_altitude

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
                   description="Subaru Telescope on Maunakea, Hawaii")

    coordinates = SkyCoord('06h45m08.9173s', '-16d42m58.017s', frame='icrs')
    target = FixedTarget(name='Sirius', coord=coordinates)

    observe_time = Time('2000-06-15 23:30:00')

    plot_altitude(target, observer, observe_time, airmass_yaxis=True)
    plt.tight_layout()
    plt.show()

.. _plots_time_window:

Specifying a time window
------------------------

If you want to see airmass plotted over a window that is not 24 hours long or
you want to control the precision of the plot, you must specify every time for
which you want to see an airmass plotted.  Therefore, an array
`~astropy.time.Time` object is necessary.

To quickly populate an `~astropy.time.Time` object with many instances of time,
use `Numpy`_ and `~astropy.units`.  See example below.

Centering the window at some time
+++++++++++++++++++++++++++++++++

To center your window at some instance in time:

 .. code-block:: python

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> import astropy.units as u
    >>> from astropy.time import Time
    >>> from astroplan.plots import plot_airmass

    >>> observe_time = Time('2000-06-15 23:30:00')
    >>> observe_time = observe_time + np.linspace(-5, 5, 55)*u.hour

    >>> plot_airmass(target, observer, observe_time)
    >>> plt.show()

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
                   description="Subaru Telescope on Maunakea, Hawaii")

    coordinates = SkyCoord('06h45m08.9173s', '-16d42m58.017s', frame='icrs')
    target = FixedTarget(name='Sirius', coord=coordinates)

    import numpy as np

    observe_time = Time('2000-06-15 23:30:00')
    observe_time = observe_time + np.linspace(-5, 5, 55)*u.hour

    plot_airmass(target, observer, observe_time)
    plt.tight_layout()
    plt.show()

Specify start and end times
+++++++++++++++++++++++++++

If you know the start and end times of your observation run, you can use a
`~astropy.time.TimeDelta` object to create an array for time input::

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from astropy.time import Time
    >>> from astroplan.plots import plot_airmass

    >>> start_time = Time('2000-06-15 20:00:00')
    >>> end_time = Time('2000-06-16 04:00:00')
    >>> delta_t = end_time - start_time
    >>> observe_time = start_time + delta_t*np.linspace(0, 1, 75)

    >>> plot_airmass(target, observer, observe_time)
    >>> plt.show()

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
                   description="Subaru Telescope on Maunakea, Hawaii")

    coordinates = SkyCoord('06h45m08.9173s', '-16d42m58.017s', frame='icrs')
    target = FixedTarget(name='Sirius', coord=coordinates)

    start_time = Time('2000-06-15 20:00:00')
    end_time = Time('2000-06-16 04:00:00')
    delta_t = end_time - start_time
    observe_time = start_time + delta_t*np.linspace(0, 1, 75)

    plot_airmass(target, observer, observe_time)
    plt.tight_layout()
    plt.show()

Plotting a quantity for multiple targets
----------------------------------------

If you want to plot airmass information for multiple targets, simply reissue
the `~astroplan.plots.plot_airmass` command, using a different
`~astroplan.FixedTarget` object as input this time. Repeat until you have as
many targets on the plot as you wish::

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from astropy.time import Time
    >>> from astroplan.plots import plot_airmass

    >>> observe_time = Time('2000-06-30 23:30:00') + np.linspace(-7.0, 5.5, 50)*u.hour

    >>> plot_airmass(target, observer, observe_time)
    >>> plot_airmass(other_target, observer, observe_time)
    >>> plot_airmass(third_target, observer, observe_time)

    >>> plt.legend(shadow=True, loc=2)
    >>> plt.show()

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
                   description="Subaru Telescope on Maunakea, Hawaii")

    coordinates = SkyCoord('06h45m08.9173s', '-16d42m58.017s', frame='icrs')
    target = FixedTarget(name='Sirius', coord=coordinates)

    coordinates = SkyCoord('02h31m49.09s', '+89d15m50.8s', frame='icrs')
    other_target = FixedTarget(name='Polaris', coord=coordinates)

    coordinates = SkyCoord('07h45m19.4s', '+28d01m35s', frame='icrs')
    third_target = FixedTarget(name='Pollux', coord=coordinates)

    observe_time = Time('2000-06-30 23:30:00') + np.linspace(-7.0, 5.5, 50)*u.hour

    plot_airmass(target, observer, observe_time)
    plot_airmass(other_target, observer, observe_time)
    plot_airmass(third_target, observer, observe_time)

    plt.legend(shadow=True, loc=2)
    plt.tight_layout()
    plt.show()

When you're ready to make a different plot, use ``ax.cla()`` to clear the
current `~matplotlib.axes.Axes` object.



.. _plots_style:

Changing style options
----------------------

The default line for time-dependent plots is solid and the default label
(should you choose to display a legend) is the name contained in the
`~astroplan.Target` object.  You can change the *linestyle*, *color*,
*label* and other plotting properties by setting the *style_kwargs* option.

.. code-block:: python

    >>> import matplotlib.pyplot as plt
    >>> from astroplan.plots import plot_airmass

    >>> sirius_styles = {'linestyle': '--', 'color': 'r'}
    >>> pollux_styles = {'color': 'g'}

    >>> plot_airmass(target, observer, observe_time, style_kwargs=sirius_styles)
    >>> plot_airmass(other_target, observer, observe_time, style_kwargs=pollux_styles)

    >>> plt.legend(shadow=True, loc=2)
    >>> plt.show()

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
                   description="Subaru Telescope on Maunakea, Hawaii")

    coordinates = SkyCoord('06h45m08.9173s', '-16d42m58.017s', frame='icrs')
    target = FixedTarget(name='Sirius', coord=coordinates)

    coordinates = SkyCoord('07h45m19.4s', '+28d01m35s', frame='icrs')
    other_target = FixedTarget(name='Pollux', coord=coordinates)

    observe_time = Time('2000-06-30 23:30:00') + np.linspace(-10, 10, 50)*u.hour

    sirius_styles = {'linestyle': '--', 'color': 'r'}
    pollux_styles = {'color': 'g'}

    plot_airmass(target, observer, observe_time, style_kwargs=sirius_styles)
    plot_airmass(other_target, observer, observe_time, style_kwargs=pollux_styles)

    plt.legend(shadow=True, loc=2)
    plt.tight_layout()
    plt.show()

See the `Matplotlib`_ documentation for information on plotting styles in line
plots.

Dark Theme Plots
++++++++++++++++

By default, `astroplan` uses the `Astropy <https://astropy.org>`__ style sheet for `Matplotlib`_ to
generate plots with more pleasing settings than provided for by the matplotlib
defaults. When using `astroplan` at night, you may prefer to make plots with
dark backgrounds, rather than the default white background, to preserve your
night vision. To do so, you may use the `astroplan` dark style sheet to
produce dark-themed plots by using the *style_sheet* keyword argument in
any plotting function:

.. code-block:: python

    >>> from astroplan.plots import dark_style_sheet, plot_airmass
    >>> plot_airmass(target, observatory, time, style_sheet=dark_style_sheet)


.. plot::

    from astropy.time import Time
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    from astroplan.plots import plot_airmass, dark_style_sheet
    from astroplan import Observer, FixedTarget

    import numpy as np
    import matplotlib.pyplot as plt

    vega = FixedTarget(coord=SkyCoord(ra=279.23473479*u.deg, dec=38.78368896*u.deg))
    apo = Observer.at_site('APO')

    plot_airmass(vega, apo, Time('2005-01-02 19:00') + np.linspace(-5, 5, 20)*u.hour,
                 style_sheet=dark_style_sheet)
    plt.tight_layout()
    plt.show()


Additional Options
++++++++++++++++++

You can also shade the background according to the darkness of the sky (light
shading for 0 degree twilight, darker shading for -18 degree twilight) with
the ``brightness_shading`` keyword, and display additional y-axis ticks on the
right side of the axis with the altitudes in degrees using the
``altitude_yaxis`` keyword:

.. code-block:: python

    >>> import matplotlib.pyplot as plt
    >>> from astropy.time import Time
    >>> from astroplan import FixedTarget, Observer
    >>> from astroplan.plots import plot_airmass

    >>> time = Time('2018-01-02 19:00')
    >>> target = FixedTarget.from_name('HD 189733')
    >>> apo = Observer.at_site('APO')
    >>> plot_airmass(target, apo, time, brightness_shading=True, altitude_yaxis=True)


.. plot::

    from astroplan.plots import plot_airmass
    import matplotlib.pyplot as plt
    from astropy.time import Time
    from astroplan import FixedTarget, Observer

    time = Time('2018-01-02 19:00')
    target = FixedTarget.from_name('HD 189733')
    apo = Observer.at_site('APO')
    plot_airmass(target, apo, time, brightness_shading=True, altitude_yaxis=True)
    plt.tight_layout()
    plt.show()



:ref:`Return to Top <plots>`


.. _plots_sky_charts:

Sky Charts
==========

Many users planning an observation run will want to see the positions of targets
with respect to their local horizon, as well as the positions of familiar stars
or other objects to act as guides.

`~astroplan.plots.plot_sky` allows you to plot the positions of targets at a
single instance or over some window of time. You make this plot the
following way:

.. code-block:: python

    >>> from astroplan.plots import plot_sky

    >>> plot_sky(target, observer, time)

.. note::

    `~astroplan.plots.plot_sky` currently produces polar plots in
    altitude/azimuth coordinates only.  Plots are centered on the observer's
    zenith.

.. seealso::

    `astroplan.Observer.altaz`

Making a plot for one instance in time
--------------------------------------

After constructing your `~astroplan.Observer` and `~astroplan.FixedTarget`
objects, construct a time input using an array of length 1.

That is, either an `~astropy.time.Time` object with an array containing one
time value (e.g., ``Time(['2000-1-1'])``) or an array containing one scalar
`~astropy.time.Time` object (e.g., ``[Time('2000-1-1')]``).

Let's say that you created `~astroplan.FixedTarget` objects for Polaris,
Altair, Vega and Deneb. To plot a map of the sky:

.. code-block:: python

    >>> import matplotlib.pyplot as plt
    >>> from astropy.time import Time
    >>> from astroplan.plots import plot_sky

    >>> observe_time = Time(['2000-03-15 15:30:00'])

    >>> polaris_style = {color': 'k'}
    >>> vega_style = {'color': 'g'}
    >>> deneb_style = {'color': 'r'}

    >>> plot_sky(polaris, observer, observe_time, style_kwargs=polaris_style)
    >>> plot_sky(altair, observer, observe_time)
    >>> plot_sky(vega, observer, observe_time, style_kwargs=vega_style)
    >>> plot_sky(deneb, observer, observe_time, style_kwargs=deneb_style)

    >>> plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
    >>> plt.show()

.. note::

    Since `~astroplan.plots.plot_sky` uses `~matplotlib.pyplot.scatter`
    (which gives the same color to different plots made on the same set of
    axes), you have to specify the color for each target via a style
    dictionary if you don't want all targets to have the same color.

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
                   description="Subaru Telescope on Maunakea, Hawaii")

    coordinates = SkyCoord('02h31m49.09s', '+89d15m50.8s', frame='icrs')
    polaris = FixedTarget(name='Polaris', coord=coordinates)
    polaris_style = {'color': 'k'}

    coordinates = SkyCoord('19h50m47.6s', '+08d52m12.0s', frame='icrs')
    altair = FixedTarget(name='Altair', coord=coordinates)

    coordinates = SkyCoord('18h36m56.5s', '+38d47m06.6s', frame='icrs')
    vega = FixedTarget(name='Vega', coord=coordinates)
    vega_style = {'color': 'g'}

    coordinates = SkyCoord('20h41m25.9s', '+45d16m49.3s', frame='icrs')
    deneb = FixedTarget(name='Deneb', coord=coordinates)
    deneb_style = {'color': 'r'}

    # Note that this is not a scalar.
    observe_time = Time(['2000-03-15 15:30:00'])

    plot_sky(polaris, observer, observe_time, style_kwargs=polaris_style)
    plot_sky(altair, observer, observe_time)
    plot_sky(vega, observer, observe_time, style_kwargs=vega_style)
    plot_sky(deneb, observer, observe_time, style_kwargs=deneb_style)

    # Note that you don't need this code block to produce the plot.
    # It reduces the plot size for the documentation.
    ax = plt.gca()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.75, box.height * 0.75])

    plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
    plt.tight_layout()
    plt.show()

Showing movement over time
--------------------------

If you want to see how your targets move over time, you need to explicitly
specify every instance in time.

Say I want to know how Altair moves in the sky over a 9-hour period:

.. code-block:: python

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> import astropy.units as u
    >>> from astroplan.plots import plot_sky

    >>> observe_time = Time('2000-03-15 17:00:00')
    >>> observe_time = observe_time + np.linspace(-4, 5, 10)*u.hour

    >>> plot_sky(altair, observer, observe_time)

    >>> plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
    >>> plt.show()

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
                   description="Subaru Telescope on Maunakea, Hawaii")

    coordinates = SkyCoord('19h50m47.6s', '+08d52m12.0s', frame='icrs')
    altair = FixedTarget(name='Altair', coord=coordinates)

    import numpy as np

    observe_time = Time('2000-03-15 17:00:00')
    observe_time = observe_time + np.linspace(-4, 5, 10)*u.hour

    plot_sky(altair, observer, observe_time)

    # Note that you don't need this code block to produce the plot.
    # It reduces the plot size for the documentation.
    ax = plt.gca()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.75, box.height * 0.75])

    plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
    plt.tight_layout()
    plt.show()

For more examples on how to populate time objects, see `~astropy.time.Time`
documentation, or :ref:`plots_time_window`.

.. note::

    Note that in the case of an object being under the horizon (or having
    negative altitude) at any of the times in your *time* input,
    `~astroplan.plots.plot_sky` will warn you.  Your object(s) will not show
    up on the plot for those particular times, but any positions above the
    horizon will still be plotted as normal.

Customizing your sky plot
-------------------------

`astroplan` plots use `Matplotlib`_ defaults, so you can customize your plots in
much the same way you tweak any `Matplotlib`_ plot.


Setting style options
+++++++++++++++++++++

The default marker for `~astroplan.plots.plot_sky` is a circle and the
default label (should you choose to display a legend) is the name contained
in the `~astroplan.Target` object.  You can change the *marker*, *color*,
*label* and other plotting properties by setting the *style_kwargs* option,
in the same way shown for the :ref:`time-dependent plots <plots_style>`.

One situation in which this is particularly useful is the plotting of guide
positions, such as a few familiar stars or any body used in calibrating your
telescope. You can also use this feature to set apart different types of targets
(e.g., high-priority, candidates for observing run, etc.).

See the `Matplotlib`_ documentation for information on plotting styles in
scatter plots.

Changing coordinate defaults
++++++++++++++++++++++++++++

As seen in the above examples, the default position of North is at the top of
the plot, and South at the bottom, with azimuth increasing counter-clockwise
(CCW), putting East to the left, and West to the right.

You can't change the position of North or South (either in the actual plotting
of the data, or the labels), but you can "flip" East/West by changing the
direction in which azimuth increases via the *north_to_east_ccw* option:

.. code-block:: python

    >>> import matplotlib.pyplot as plt
    >>> from astroplan.plots import plot_sky

    >>> guide_style = {'marker': '*'}

    >>> plot_sky(polaris, observer, observe_time, snorth_to_east_ccw=False, style_kwargs=guide_style)
    >>> plot_sky(altair, observer, observe_time, north_to_east_ccw=False)

    >>> plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
    >>> plt.show()

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
                   description="Subaru Telescope on Maunakea, Hawaii")

    coordinates = SkyCoord('02h31m49.09s', '+89d15m50.8s', frame='icrs')
    polaris = FixedTarget(name='Polaris', coord=coordinates)

    coordinates = SkyCoord('19h50m47.6s', '+08d52m12.0s', frame='icrs')
    altair = FixedTarget(name='Altair', coord=coordinates)

    import numpy as np

    observe_time = Time('2000-03-15 17:00:00') + np.linspace(-4, 5, 10)*u.hour

    guide_style = {'marker': '*'}

    plot_sky(polaris, observer, observe_time, north_to_east_ccw=False,
             style_kwargs=guide_style)
    plot_sky(altair, observer, observe_time, north_to_east_ccw=False)

    # Note that you don't need this code block to produce the plot.
    # It reduces the plot size for the documentation.
    ax = plt.gca()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.75, box.height * 0.75])

    plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
    plt.tight_layout()
    plt.show()

Some observatories may need to offset or rotate the azimuth labels due to their
particular telescope setup.

To do this, set *az_label_offset* equal to the number of degrees by which you
wish to rotate the labels.  By default, *az_label_offset* is set to 0 degrees.
A positive offset is in the same direction as azimuth increase (see the
*north_to_east_ccw* option):

.. code-block:: python

    >>> import matplotlib.pyplot as plt
    >>> from astroplan.plots import plot_sky

    >>> guide_style = {'marker': '*'}

    >>> plot_sky(polaris, observer, observe_time, style_kwargs=guide_style, az_label_offset=180.0*u.deg)
    >>> plot_sky(altair, observer, observe_time, az_label_offset=180.0*u.deg)

    >>> plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
    >>> plt.show()

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
                   description="Subaru Telescope on Maunakea, Hawaii")

    coordinates = SkyCoord('02h31m49.09s', '+89d15m50.8s', frame='icrs')
    polaris = FixedTarget(name='Polaris', coord=coordinates)

    coordinates = SkyCoord('19h50m47.6s', '+08d52m12.0s', frame='icrs')
    altair = FixedTarget(name='Altair', coord=coordinates)

    import numpy as np

    observe_time = Time('2000-03-15 17:00:00') + np.linspace(-4, 5, 10)*u.hour

    guide_style = {'marker': '*'}

    plot_sky(polaris, observer, observe_time, style_kwargs=guide_style,
             az_label_offset=180.0*u.deg)
    plot_sky(altair, observer, observe_time, az_label_offset=180.0*u.deg)

    # Note that you don't need this code block to produce the plot.
    # It reduces the plot size for the documentation.
    ax = plt.gca()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.75, box.height * 0.75])

    plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
    plt.tight_layout()
    plt.show()

.. note::

    The *az_label_offset* option does not rotate the actual positions on the
    plot, but simply the theta grid labels (which are drawn regardless of
    gridline presence).  Since labels are drawn with every call to
    `~astroplan.plots.plot_sky`, we recommend you use the same
    *az_label_offset* argument for every target on the same plot.

    It is not advised that most users change this option, as it may **appear**
    that your alt/az data does not coincide with the definition of altazimuth
    (local horizon) coordinate system.

You can turn off the grid lines by setting the *grid* option to *False*:

.. code-block:: python

    >>> import matplotlib.pyplot as plt
    >>> from astroplan.plots import plot_sky

    >>> guide_style = {'marker': '*'}

    >>> plot_sky(polaris, observer, observe_time, style_kwargs=guide_style, grid=False)
    >>> plot_sky(altair, observer, observe_time, grid=False)

    >>> plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
    >>> plt.show()

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
                   description="Subaru Telescope on Maunakea, Hawaii")

    coordinates = SkyCoord('02h31m49.09s', '+89d15m50.8s', frame='icrs')
    polaris = FixedTarget(name='Polaris', coord=coordinates)

    coordinates = SkyCoord('19h50m47.6s', '+08d52m12.0s', frame='icrs')
    altair = FixedTarget(name='Altair', coord=coordinates)

    import numpy as np

    observe_time = Time('2000-03-15 17:00:00') + np.linspace(-4, 5, 10)*u.hour

    guide_style = {'marker': '*'}

    plot_sky(polaris, observer, observe_time, style_kwargs=guide_style,
             grid=False)
    plot_sky(altair, observer, observe_time, grid=False)

    # Note that you don't need this code block to produce the plot.
    # It reduces the plot size for the documentation.
    ax = plt.gca()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.75, box.height * 0.75])

    plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
    plt.tight_layout()
    plt.show()

.. note::

    Since grids are redrawn with every call to `~astroplan.plots.plot_sky`,
    you must set ``grid=False`` for every target in the same plot.

Other tweaks
++++++++++++

You can easily change other plot attributes by acting on the returned
`matplotlib.axes.Axes` object or via `matplotlib.pyplot` calls (e.g.,
``plt.figure``, ``plt.rc``, etc.).

For instance, you can increase the size of your plot and its font:

.. code-block:: python

    >>> # Set the figure size/font before you issue the plotting command.
    >>> plt.figure(figsize=(8,6))
    >>> plt.rc('font', size=14)

    >>> plot_sky(polaris, observer, observe_time)
    >>> plot_sky(altair, observer, observe_time)

    >>> plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
    >>> plt.show()

    >>> # Change font size back to default once done plotting.
    >>> plt.rc('font', size=12)

:ref:`Return to Top <plots>`

Miscellaneous
+++++++++++++

The easiest way to reuse the `~matplotlib.axes.Axes` object that is the base of
your plots is to just let `astroplan`'s plotting functions take care of it in
the background.  You do, however, have the option of explicitly passing in a
named axis, assuming that you have created the appropriate type for the
particular plot you want.

We can explicitly give a name to the `~matplotlib.axes.Axes` object returned
by `~astroplan.plots.plot_sky` when plotting Polaris and reuse it to plot
Altair:

.. code-block:: python

    >>> my_ax = plot_sky(polaris, observer, observe_time, style_kwargs=guide_style)
    >>> plot_sky(altair, observer, observe_time, my_ax)

We can also create a `~matplotlib.axes.Axes` object entirely outside of
`~astroplan.plots.plot_sky`, then pass it in:

.. code-block:: python

    >>> my_ax = plt.gca(projection='polar')
    >>> plot_sky(polaris, observer, observe_time, my_ax)

Passing in named `matplotlib.axes.Axes` objects comes in handy when you want to
make multiple plots:

.. code-block:: python

    >>> from astroplan.plots import plot_sky
    >>> import matplotlib.pyplot as plt

    >>> my_ax = plot_sky(polaris, observer, observe_time, style_kwargs=polaris_style)
    >>> plot_sky(altair, observer, observe_time, my_ax, style_kwargs=altair_style)
    >>> plt.legend(loc='center left', bbox_to_anchor=(1.3, 0.5))

    >>> # Note that this plt.show (or another action, such as saving a figure) is critical in maintaining two separate plots.
    >>> plt.show()

    >>> other_ax = plot_sky(vega, observer, observe_time, style_kwargs=vega_style)
    >>> plot_sky(deneb, observer, observe_time, other_ax, style_kwargs=deneb_style)
    >>> plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
    >>> plt.show()

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
                   description="Subaru Telescope on Maunakea, Hawaii")

    coordinates = SkyCoord('02h31m49.09s', '+89d15m50.8s', frame='icrs')
    polaris = FixedTarget(name='Polaris', coord=coordinates)

    coordinates = SkyCoord('19h50m47.6s', '+08d52m12.0s', frame='icrs')
    altair = FixedTarget(name='Altair', coord=coordinates)
    altair_style = {'color': 'b'}

    coordinates = SkyCoord('18h36m56.5s', '+38d47m06.6s', frame='icrs')
    vega = FixedTarget(name='Vega', coord=coordinates)
    vega_style = {'color': 'g'}

    coordinates = SkyCoord('20h41m25.9s', '+45d16m49.3s', frame='icrs')
    deneb = FixedTarget(name='Deneb', coord=coordinates)
    deneb_style = {'color': 'r'}

    observe_time = Time(['2000-03-15 15:30:00'])

    my_ax = plot_sky(polaris, observer, observe_time)
    plot_sky(altair, observer, observe_time, my_ax, style_kwargs=altair_style)

    # Note that you don't need this code block to produce the plot.
    # It reduces the plot size for the documentation.
    ax = plt.gca()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.75, box.height * 0.75])

    plt.legend(loc='center left', bbox_to_anchor=(1.3, 0.5))

    # Note that this plt.show (or another action, such as saving a figure) is
    # critical in maintaining two separate plots.
    plt.tight_layout()
    plt.show()

    other_ax = plot_sky(vega, observer, observe_time, style_kwargs=vega_style)
    plot_sky(deneb, observer, observe_time, other_ax, style_kwargs=deneb_style)

    # Note that you don't need this code block to produce the plot.
    # It reduces the plot size for the documentation.
    ax = plt.gca()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.75, box.height * 0.75])

    plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
    plt.tight_layout()
    plt.show()

:ref:`Return to Top <plots>`

.. _finder_image:

Finder Chart/Image
==================

`astroplan` includes a function for generating quick finder images from
Python, `~astroplan.plots.plot_finder_image`, by querying for images from
sky surveys centered on a `~astroplan.FixedTarget`. This function depends on
`astroquery`_ (in addition to `Matplotlib`_). In this example, we'll quickly
make a finder image centered on The Crab Nebula (M1):

.. code-block:: python

    >>> from astroplan.plots import plot_finder_image
    >>> from astroplan import FixedTarget
    >>> import matplotlib.pyplot as plt

    >>> messier1 = FixedTarget.from_name("M1")
    >>> ax, hdu = plot_finder_image(messier1)
    >>> plt.show()

.. plot::

    from astroplan import FixedTarget
    from astroplan.plots import plot_finder_image
    import matplotlib.pyplot as plt

    messier1 = FixedTarget.from_name("M1")
    ax, hdu = plot_finder_image(messier1)
    plt.tight_layout()
    plt.show()


:ref:`Return to Top <plots>`
