.. doctest-skip-all

.. _precision:

**************************************************************
Why is my target above/below the horizon at the rise/set time?
**************************************************************

Rise/set/meridian transit calculations in `astroplan` are designed to be fast
while achieving a precision comparable to what can be predicted given the
affects of the changing atmosphere. As a result, there may be some
counter-intuitive behavior in `astroplan` methods like
`astroplan.Observer.target_rise_time`, `astroplan.Observer.target_set_time` and
`astroplan.Observer.target_meridian_transit_time`, that can lead to small
changes in the numerical values of these computed timed (of order seconds).

For example, to calculate the rise time of Sirius, you might write::

    from astroplan import Observer, FixedTarget
    from astropy.time import Time

    # Set up observer, target, and time
    keck = Observer.at_site("Keck")
    sirius = FixedTarget.from_name("Sirius")
    time = Time('2010-05-11 06:00:00')

    # Find rise time of Sirius at Keck nearest to `time`
    rise_time = keck.target_rise_time(time, sirius)

You might expect the altitude of Sirius to be zero degrees at ``rise_time``,
i.e. Sirius will be on the horizon, but this is not the case::

    >>> altitude_at_rise = keck.altaz(rise_time, sirius).alt
    >>> print(altitude_at_rise.to('arcsec'))
    2.70185arcsec

The altitude that you compute on your machine may be different from the number
above by a small amount – for a detailed explanation on where the difference
arises from, see :doc:`iers`. The rise and set time methods use the following
approximation:

* A time series of altitudes for the target is computed at times near ``time``

* The two times when the target is nearest to the horizon are identified, and a
  linear interpolation is done between those times to find the horizon-crossing

This method has a precision of a few arcseconds, so your targets may be slightly
above or below the horizon at their rise or set times.

***********************************
How are sunrise and sunset defined?
***********************************

By default, `~astroplan.Observer.sun_set_time` will compute an approximation
to the moment when the center of the solar disk is on the horizon. This differs
slightly from the conventional definition used by USNO, for example, which
returns the time when the solar disk center is at -0.8333 degrees altitude to
account for the solar radius and atmospheric refraction. If you want to
accurately reproduce the sun rise and set times computed by USNO,
you can call `~astroplan.Observer.sun_set_time` with the keyword argument
``horizon=-0.8333*u.deg``, like so::

    >>> from astroplan import Observer
    >>> from astropy.time import Time
    >>> import astropy.units as u

    >>> mmt = Observer.at_site('mmt', pressure=0*u.bar)

    >>> # USNO time from the MMTO Almanac:
    >>> # http://www.mmto.org/sites/default/files/almanac_2019.pdf
    >>> usno_sunset = Time('2019-01-01 00:31')

    >>> # Compute equivalent time with astroplan
    >>> astroplan_sunset = mmt.sun_set_time(usno_sunset - 10*u.min,
    ...                                     horizon=-0.8333*u.deg, which='next')

    >>> abs(usno_sunset - astroplan_sunset) < 1 * u.min
    True
