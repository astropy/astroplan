.. _getting_started:

***************
Getting Started
***************

General Guidelines
==================

`Astroplan` is based on `Astropy`_ and was built around the creation of "objects"
that contain all the information needed to perform certain tasks.

You, the user, will create and manipulate these objects to plan your
observation.  Don't worry--it's easier than it sounds!

For instance, a `Target` object contains all the data associated with observing
targets--RA, Dec, etc.::

    from astropy.coordinates import SkyCoord
    from astroplan import FixedTarget

    coordinates = SkyCoord('19h50m47.6s', '+08d52m12.0s', frame='icrs')
    altair = FixedTarget(name='Altair', coord=coordinates)

Similarly, an `Observer` object has longitude/latitude coordinates of the
observatory or telescope where you are planning on making your observations,
as well as elevation and other location information::

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

`Astroplan` makes heavy use of certain pieces of `Astropy` machinery, such as
the representation of dates/times in `Time` objects::

    from astropy.time import Time

    time = Time(['2015-06-15 22:00:00'])

To pull some piece of information from or to manipulate an `astroplan` object,
you will generally issue a command such as::

    right_ascension = altair.ra

or::

    p_angle = observer.parallactic_angle(time, altair)

Doing More
==========

Now that you know the basics of working with `astroplan`, check out our
:ref:`tutorials` and :ref:`how_to` pages for detailed examples.
