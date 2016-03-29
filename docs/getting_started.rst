.. _getting_started:

***************
Getting Started
***************

General Guidelines
==================

`astroplan` is based on `Astropy`_ and was built around the creation of Python
objects that contain all the information needed to perform certain tasks.  You,
the user, will create and manipulate these objects to plan your observation. For
instance, an `~astroplan.Target` object contains information associated with
targets, such as right ascension, declination, etc.

Objects representing celestial bodies like stars (which, if we ignore proper
motion, are fixed on the celestial sphere) are created (or "instantiated") via
an `~astroplan.FixedTarget` object::

    from astropy.coordinates import SkyCoord
    from astroplan import FixedTarget

    coordinates = SkyCoord('19h50m47.6s', '+08d52m12.0s', frame='icrs')
    altair = FixedTarget(name='Altair', coord=coordinates)

Similarly, an `~astroplan.Observer` object contains information about the
observatory, telescope or place where you are observing, such as longitude,
latitude, elevation and other optional parameters.  You can initialize an
`~astroplan.Observer` object via the `~astroplan.get_site` function, which
accesses the site list built into `astroplan` (see
`~astroplan.get_site_names`)::

    from astroplan import get_site
    observer = get_site('subaru')

Or you can specify your own location parameters::

    import astropy.units as u
    from astropy.coordinates import EarthLocation
    from pytz import timezone
    from astroplan import Observer

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

`astroplan` makes heavy use of certain `Astropy`_ machinery, including the
`~astropy.coordinates` objects and transformations and
`~astropy.units`. Most importantly for basic use of `astroplan` is the
representation of dates/times as `~astropy.time.Time` objects (note that
these are in the UTC timezone by default)::

    from astropy.time import Time
    time = Time(['2015-06-16 06:00:00'])

Since `astroplan` objects are Python objects, manipulating them or accessing
attributes follows Python syntax and conventions.  See Python documentation on
`objects <https://docs.python.org/2/tutorial/classes.html#instance-objects>`_
for more information.

Doing More
==========

Now that you know the basics of working with `astroplan`, check out our
:ref:`tutorials` page for high-level examples of using `astroplan`, as well as
the :ref:`api` section for more exhaustive documentation and lower-level usage
examples.
