:orphan:

.. doctest-skip-all

**********************************
Using and updating IERS Bulletin A
**********************************

Contents
========

* :ref:`iers_a-Background`
* :ref:`iers_a-What_are_IERS_Bulletins_A_and_B`
* :ref:`iers_a-How_do_I_get_IERS_Bulletin_A_for_astroplan?`

.. _iers_a-Background:

Background
==========

Keeping accurate, consistent time systems is rough. We intuitively want
time to be described by a continuous cyclic clock like the clock on your
wall or the calendar, but the Earth's rate of rotation varies measurably
on human timescales. As a result, the time measured by an atomic clock and
the time that you would infer from measuring how long it has been since the
last solar noon are constantly becoming offset from one another with
stochastic changes in Earth's moment of inertia and slowly acting tidal forces.

For this reason, there is a time system that keeps track of seconds like an
atomic clock, UTC, one that matches up with the Earth's rotation, UT1. The
difference between UT1 and UTC is constantly changing as the Earth's rotation
changes and as leap seconds get added to UTC to compensate.

In order to accurately predict the apparent position of a celestial object from
the Earth, for example in altitude and azimuth coordinates, one needs to know
the orientation of the Earth at any point in time, and since the Earth's
rotation does not change in a predictable way, we must rely on observations
of the Earth's orientation as a function of time to accurately calculate
UT1-UTC.

.. _iers_a-What_are_IERS_Bulletins_A_and_B:

What are IERS Bulletins A and B?
================================

The International Earth Rotation and Reference Systems Service (IERS) is
responsible for measuring the Earth's orientation as a function of time, and
making predictions for the Earth's orientation in the near future (accounting
for scheduled leap seconds). The data products released by the IERS used by
astroplan are the IERS Bulletins A and B. IERS Bulletin B is a table with Earth
orientation observations from the recent few decades up through nearly the
present time. `Astropy` relies on IERS B to compute UT1-UTC, and by default
will throw an error if you try to compute UT1-UTC for a time that is outside
the bounds of the IERS Bulletin B table.

.. _iers_a-How_do_I_get_IERS_Bulletin_A_for_astroplan?:

How do I get IERS Bulletin A for astroplan?
===========================================

Without downloading IERS Bulletin A, astroplan simply approximates UT1-UTC=0,
which will lead to lower precision position and time calculations, on the order
of arcseconds or seconds.

To download the IERS Bulletin A table for the first time, or to refresh the
cached version that you already have, simply run::

    from astroplan import download_IERS_A
    download_IERS_A()

