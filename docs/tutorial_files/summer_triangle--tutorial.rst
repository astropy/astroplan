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

    time = Time(['2015-06-16 12:00:00'])

Observable?
===========

Next, it would be handy to know if our targets are visible from Subaru at the
time we settled on.  In other words--are they above the horizon?:

.. code-block:: ipython

    In [1]: subaru.can_see(time, altair)
    Out[1]: array([ True], dtype=bool)

    In [2]: subaru.can_see(time, vega)
    Out[2]: array([ True], dtype=bool)

    In [3]: subaru.can_see(time, deneb)
    Out[3]: array([ True], dtype=bool)

...They are!

However, we may want to know more details about their observability--when do
they rise/set, is the Sun or the Moon out at that time as well?

The best thing to do would probably be to find a window of time for tonight
during which all three of our targets are above the horizon **and** the Sun is
below the horizon (we can worry about the Moon later).

So, let's calculate our target rise and set times::

    altair_rise = subaru.calc_rise(time, altair)
    altair_set = subaru.calc_set(time, altair)

    vega_rise = subaru.calc_rise(time, vega)
    vega_set = subaru.calc_set(time, vega)

    deneb_rise = subaru.calc_rise(time, deneb)
    deneb_set = subaru.calc_set(time, deneb)
