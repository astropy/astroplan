.. _exoplanet_tutorial:

******************************************************
Observing Transiting Exoplanets and Eclipsing Binaries
******************************************************

.. note::
    The exoplanets module is new and under development. The API may change in
    upcoming versions of astroplan, and pull requests are welcome!


Contents
========

* :ref:`exoplanets-transit_times`
* :ref:`exoplanets-astroquery`

.. _exoplanets-transit_times:

Transit/secondary eclipse Times
===============================

We can define the properties of an eclipsing system, such as an eclipsing binary
or transiting exoplanet, using the `~astroplan.EclipsingSystem` object. Let's
make an instance for the transiting exoplanet HD 209458 b, which has a period
of 3.52474859 days, mid-transit time of JD=2452826.628514, and transit duration
of 0.1277:

.. code-block:: python

    >>> from astropy.time import Time
    >>> import astropy.units as u
    >>> from astroplan import EclipsingSystem

    >>> epoch = Time(2452826.628514, format='jd')
    >>> period = 3.52474859 * u.day
    >>> duration = 0.1277 * u.day

    >>> hd209458 = EclipsingSystem(epoch=epoch, period=period, duration=duration,
                                   name='HD 209458 b')

Let's say we're observing on 2016 January 1, 00:00 UTC. We can compute the next
transit and secondary eclipse using the
`~astroplan.EclipsingSystem.next_primary_eclipse_time` and
`~astroplan.EclipsingSystem.next_secondary_eclipse_time` methods, respectively:

.. code-block:: python

    >>> observing_time = Time('2016-01-01 00:00')
    >>> hd209458.next_primary_eclipse_time(observing_time)
    <Time object: scale='utc' format='iso' value=['2016-01-03 16:16:09.848']>

    >>> hd209458.next_secondary_eclipse_time(observing_time)
    <Time object: scale='utc' format='iso' value=['2016-01-01 21:58:20.708']>

You can compute the next ten mid-transit times with the ``n_eclipses`` keyword:

.. code-block:: python

    >>> hd209458.next_primary_eclipse_time(observing_time, n_eclipses=5)
    <Time object: scale='utc' format='iso' value=['2016-01-03 16:16:09.848' '2016-01-07 04:51:48.126'
                                                  '2016-01-10 17:27:26.404' '2016-01-14 06:03:04.682'
                                                  '2016-01-17 18:38:42.960']>

It's often useful to know the ingress and egress times of the next transits
when planning observations, which you can find with
`~astroplan.EclipsingSystem.next_primary_ingress_egress_time`

    >>> hd209458.next_primary_ingress_egress_time(observing_time, n_eclipses=3)
    <Time object: scale='utc' format='jd' value=[[ 2457391.11404175  2457391.24174175]
                                                 [ 2457394.63879034  2457394.76649034]
                                                 [ 2457398.16353893  2457398.29123893]]>

.. _exoplanets-astroquery:

Getting exoplanet properties via astroquery
===========================================

#As of `astroquery`_ version 0.X, users can query for transiting exoplanet
#properties from the Exoplanet Orbit Database (exoplanets.org) and the NASA
#Exoplanet Science Institute (NExScI) Exoplanet Archive.

