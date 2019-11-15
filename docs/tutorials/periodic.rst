.. include:: ../references.txt

.. doctest-skip-all

.. _periodic_tutorial:

******************************************************
Observing Transiting Exoplanets and Eclipsing Binaries
******************************************************

.. note::
    The ``periodic`` module is new and under development. The API may change in
    upcoming versions of astroplan, and pull requests are welcome!

.. warning::

    There are currently two major caveats in the implementation of
    `~astroplan.EclipsingSystem`. The secondary eclipse time approximation is
    only accurate when the orbital eccentricity is small, and the eclipse
    times are computed without any barycentric corrections. The current
    implementation should only be used for approximate mid-eclipse times for
    low eccentricity orbits, with event durations longer than the
    barycentric correction error (<=16 minutes).

Contents
========

* :ref:`periodic-transit_times`
* :ref:`periodic-transit_times_via_astroquery`
* :ref:`periodic-observable_transits`
* :ref:`periodic-phase_constraint`

.. _periodic-transit_times:

Transit/Primary and secondary eclipse times
===========================================

We can define the properties of an eclipsing system, such as an eclipsing binary
or transiting exoplanet, using the `~astroplan.EclipsingSystem` object. Let's
make an instance for the transiting exoplanet HD 209458 b, which has a period
of 3.52474859 days, mid-transit time of JD=2452826.628514, and transit duration
of 0.1277:

.. code-block:: python

    >>> from astropy.time import Time
    >>> import astropy.units as u
    >>> from astroplan import EclipsingSystem

    >>> primary_eclipse_time = Time(2452826.628514, format='jd')
    >>> orbital_period = 3.52474859 * u.day
    >>> eclipse_duration = 0.1277 * u.day

    >>> hd209458 = EclipsingSystem(primary_eclipse_time=primary_eclipse_time,
    ...                            orbital_period=orbital_period, duration=eclipse_duration,
    ...                            name='HD 209458 b')

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
`~astroplan.EclipsingSystem.next_primary_ingress_egress_time`:

.. code-block:: python

    >>> hd209458.next_primary_ingress_egress_time(observing_time, n_eclipses=3)
    <Time object: scale='utc' format='jd' value=[[ 2457391.11404175  2457391.24174175]
                                                 [ 2457394.63879034  2457394.76649034]
                                                 [ 2457398.16353893  2457398.29123893]]>

And remember - in the current implementation, all eclipse times are computed
without any barycentric corrections, and the secondary eclipse time
approximation is only accurate when the orbital eccentricity is small.

.. _periodic-transit_times_via_astroquery:

Transit times via astroquery
============================

The development version of `astroquery`_ allows users to query for properties of
known exoplanets with three different services:
`~astroquery.exoplanet_orbit_database`, `~astroquery.nasa_exoplanet_archive`,
and `~astroquery.open_exoplanet_catalogue`. In the example below, we will query
for the properties of the transiting exoplanet TRAPPIST-1 b with astroquery, and
calculate the times of the next three transits with
`~astroplan.EclipsingSystem`.

.. code-block:: python

    >>> # NASA Exoplanet Archive for planet properties
    >>> import astropy.units as u
    >>> from astropy.time import Time
    >>> from astroquery.nasa_exoplanet_archive import NasaExoplanetArchive
    >>> planet_properties = NasaExoplanetArchive.query_planet('TRAPPIST-1 b', all_columns=True)

    >>> # get relevant planet properties
    >>> epoch = Time(planet_properties['pl_tranmid'], format='jd')
    >>> period = planet_properties['pl_orbper']
    >>> transit_duration = planet_properties['pl_trandur'] * u.day

    >>> # Create an EclipsingSystem object for HD 209458
    >>> from astroplan import EclipsingSystem
    >>> trappist1b = EclipsingSystem(primary_eclipse_time=epoch, orbital_period=period,
    ...                              duration=transit_duration)

    >>> # Calculate next three mid-transit times which occur after ``obs_time``
    >>> obs_time = Time('2017-01-01 12:00')
    >>> trappist1b.next_primary_eclipse_time(obs_time, n_eclipses=3)
    <Time object: scale='utc' format='iso' value=['2017-01-02 15:17:40.205' '2017-01-04 03:33:19.443'
     '2017-01-05 15:48:58.681']>

.. _periodic-observable_transits:

When is the next observable transit?
====================================

Let's continue with the example from above, and now let's calculate all
mid-transit times of HD 209458 b which are observable from Apache Point
Observatory, when the target is above 30 degrees altitude, and in the "A" half
of the night (roughly between sunset and midnight). First we need to create a
`~astroplan.FixedTarget` object for the star, which contains the sky coordinate,
and the `~astroplan.EclipsingSystem` object, which defines the transit time,
period and duration:

.. code-block:: python

    >>> from astroplan import FixedTarget, Observer, EclipsingSystem
    >>> apo = Observer.at_site('APO', timezone='US/Mountain')
    >>> target = FixedTarget.from_name("HD 209458")

    >>> primary_eclipse_time = Time(2452826.628514, format='jd')
    >>> orbital_period = 3.52474859 * u.day
    >>> eclipse_duration = 0.1277 * u.day

    >>> hd209458 = EclipsingSystem(primary_eclipse_time=primary_eclipse_time,
    ...                            orbital_period=orbital_period, duration=eclipse_duration,
    ...                            name='HD 209458 b')

Then we compute a list of mid-transit times over the next year:

.. code-block:: python

    >>> n_transits = 100  # This is the roughly number of transits per year
    >>> obs_time = Time('2017-01-01 12:00')
    >>> midtransit_times = hd209458.next_primary_eclipse_time(obs_time, n_eclipses=n_transits)

Finally, we can check if the target is observable at each transit time, given
our constraints on the altitude of the target (`~astroplan.AltitudeConstraint`)
and the time of observations (`~astroplan.LocalTimeConstraint` and
`~astroplan.AtNightConstraint`) with the function`~astroplan.is_event_observable`:

.. code-block:: python

    >>> from astroplan import (PrimaryEclipseConstraint, is_event_observable
    ...                        AtNightConstraint, AltitudeConstraint, LocalTimeConstraint)
    >>> import datetime as dt
    >>> import astropy.units as u
    >>> min_local_time = dt.time(19, 0)  # 19:00 local time at APO (7pm)
    >>> max_local_time = dt.time(0, 0)  # 00:00 local time at APO (midnight)
    >>> constraints = [AtNightConstraint.twilight_civil(),
    ...                AltitudeConstraint(min=30*u.deg),
    ...                LocalTimeConstraint(min=min_local_time, max=max_local_time)]

    >>> is_event_observable(constraints, apo, target, times=midtransit_times)
    array([[ True, False,  True, ...,  True, False,  True, False]], dtype=bool)

In the above example, we only checked that the star is observable at the
mid-transit time. If you were planning to do transit photometry of HD 209458 b,
you might want to be sure that the entire transit is observable. Let's look
for only completely observable transits:

.. code-block:: python

    >>> ing_egr = hd209458.next_primary_ingress_egress_time(observing_time, n_eclipses=n_transits)
    >>> is_event_observable(constraints, apo, target, times_ingress_egress=ing_egr)
    array([[False, False, False, ...,  True, False, False, False]], dtype=bool)

Note that several of the transits that were observable at their mid-transit time
are not observable at both the ingress and egress times, and therefore are
not observable in the computation above.

.. _periodic-phase_constraint:

Orbital Phase Constraint
========================

It is often useful to plan observations as a function of orbital phase. You can
calculate the orbital phase of an eclipsing or non-eclipsing system with the
`~astroplan.PeriodicEvent` object, which you specify with an epoch and period.
Let's create a `~astroplan.PeriodicEvent` object for an imagined binary star:

.. code-block:: python

    >>> from astroplan import PeriodicEvent
    >>> import astropy.units as u
    >>> from astropy.time import Time

    >>> epoch = Time(2456001, format='jd')  # reference time of periodic event
    >>> period = 3.25 * u.day  # period of periodic event
    >>> duration = 2 * u.hour  # duration of event

    >>> binary_system = PeriodicEvent(epoch=epoch, period=period)

Now let's determine when we can observe the binary given some observing
constraints. We want to measure the binary's radial velocity at orbital phases
between 0.4 and 0.6, while observing between astronomical twilights, and while
the target is above 40 degrees altitude, for an observer in Greenwich, England
on the night of January 1, 2017. For this task we can use the
`~astroplan.PhaseConstraint` (learn more about the constraints module in
:doc:`constraints`):

.. code-block:: python

    >>> from astropy.coordinates import SkyCoord
    >>> from astroplan import FixedTarget, Observer, is_observable
    >>> target = FixedTarget(SkyCoord(ra=42*u.deg, dec=42*u.deg), name='Target')
    >>> greenwich = Observer.at_site("Greenwich")
    >>> start_time = Time('2017-01-01 01:00')
    >>> end_time = Time('2017-01-01 06:00')

    >>> from astroplan import PhaseConstraint, AtNightConstraint, AltitudeConstraint
    >>> constraints = [PhaseConstraint(binary_system, min=0.4, max=0.6),
    ...                AtNightConstraint.twilight_astronomical(),
    ...                AltitudeConstraint(min=40 * u.deg)]
    >>> is_observable(constraints, greenwich, target, time_range=[start_time, end_time])
    array([ True], dtype=bool)
