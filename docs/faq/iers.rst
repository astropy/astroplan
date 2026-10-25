.. doctest-skip-all

.. _iers:

**************************************************
What are the IERS tables and how do I update them?
**************************************************

Contents
========

* :ref:`iers-Background`
* :ref:`iers-What_are_IERS_Bulletins_A_and_B`
* :ref:`iers-How_do_I_get_IERS_Bulletin_A_for_astroplan?`
* :ref:`iers-Older_astropy_versions`

.. _iers-Background:

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
atomic clock, TAI, a more commonly used one that works like an atomic clock
with a varying offset of some integer number of seconds, UTC, and one that
matches up with the Earth's rotation, UT1. The difference between UT1 and UTC is
constantly changing as the Earth's rotation changes and as leap seconds get
added to UTC to compensate.

In order to accurately predict the apparent position of a celestial object from
the Earth, for example in altitude and azimuth coordinates, one needs to know
the orientation of the Earth at any point in time, and since the Earth's
rotation does not change in a predictable way, we must rely on observations
of the Earth's orientation as a function of time to accurately calculate
UT1-UTC.

.. _iers-What_are_IERS_Bulletins_A_and_B:

What are IERS Bulletins A and B?
================================

The `International Earth Rotation and Reference Systems Service (IERS)
<https://www.iers.org/>`_ is responsible for measuring the Earth's orientation as
a function of time, and making predictions for the Earth's orientation in the
near future (accounting for scheduled leap seconds). The data products released
by the IERS used by astroplan are the IERS Bulletins A and B.

* IERS Bulletin B is a table with Earth orientation observations from the last
  few decades up through nearly the present time. `Astropy <https://astropy.org>`__ relies on IERS B to
  compute UT1-UTC, and by default will raise an error if you try to compute
  UT1-UTC for a time that is outside the bounds of the IERS Bulletin B table
  (see the `Astropy <https://astropy.org>`__ docs on the `UT1/UTC transformation offsets
  <http://astropy.readthedocs.io/en/latest/time/index.html?highlight=iers#transformation-offsets>`_
  for more details), like this::

    >>> from astropy.time import Time
    >>> Time('2040-01-01', scale='utc').ut1    # Convert from UTC to UT1
    IndexError: (some) times are outside of range covered by IERS table.

* IERS Bulletin A encompasses the observations of recent Earth orientation
  contained in Bulletin B, while also making extrapolations into the past,
  before the Bulletin B tables begin, and into the future, after the Bulletin
  B tables end. The future predictions include leap second additions scheduled
  to be added.

Let's plot the UT1-UTC from IERS Bulletins A and B to show the difference using
astropy's IERS machinery:

.. plot::

    import matplotlib.pyplot as plt
    import numpy as np
    import astropy.units as u
    from astropy.time import Time

    # Load IERS Bulletin A (bundled with astropy in astropy-iers-data, or
    # downloaded and cached by astropy if that is out of date) and the
    # Bulletin B table bundled with astropy
    from astropy.utils.iers import IERS_B, FROM_IERS_B, IERS_Auto
    iers_a = IERS_Auto.open()
    iers_b = IERS_B.open()

    # Get a range of times to plot from 1990-2022
    time_range = Time("1990-01-01") + np.arange(0, 30, 0.2)*u.year

    # Calculate the difference between UTC and UT1 at those times,
    # allowing times "outside of the table"
    DUT1_a, success_a = time_range.get_delta_ut1_utc(return_status=True,
                                                     iers_table=iers_a)
    DUT1_b, success_b = time_range.get_delta_ut1_utc(return_status=True,
                                                     iers_table=iers_b)

    # Compare input times to the times available in the table. For details, see
    # https://docs.astropy.org/en/stable/utils/iers.html
    measurements_from_b = (success_b == FROM_IERS_B)

    # Make a plot of the time difference
    fig, ax = plt.subplots(figsize=(10,8))
    ax.axhline(0, color='gray', ls='--', lw=2)

    ax.plot(time_range.plot_date,
            DUT1_a, '-', lw=2, label='IERS Bulletin A + extrapolation')
    ax.plot(time_range.plot_date[measurements_from_b],
            DUT1_b[measurements_from_b], 'r--', lw=2, label='IERS Bulletin B')
    ax.xaxis_date()
    ax.set(xlabel='Year', ylabel='UT1-UTC [seconds]')
    ax.legend(loc='upper right')
    plt.show()


.. _iers-How_do_I_get_IERS_Bulletin_A_for_astroplan?:

How do I get IERS Bulletin A for astroplan?
===========================================

You don't need to do anything specific to astroplan: astroplan uses whatever
IERS tables astropy uses. astropy ships both IERS bulletins in the
``astropy-iers-data`` package, which is installed with astropy and updated
often. To get the latest tables, upgrade it::

    pip install --upgrade astropy-iers-data

or, with conda::

    conda update astropy-iers-data

If a calculation needs Earth orientation predictions and the bundled
Bulletin A is more than ``astropy.utils.iers.conf.auto_max_age`` days (30 by
default) old, astropy downloads a newer one from the IERS and caches it. To
work offline, or to never download the tables, set::

    from astropy.utils import iers
    iers.conf.auto_download = False

See the astropy documentation on `IERS data access
<https://docs.astropy.org/en/stable/utils/iers.html>`_ for all the options.

.. note::

    ``astroplan.download_IERS_A`` is deprecated since astroplan 0.11 and
    does nothing. It was a workaround for limitations of astropy's IERS
    handling that astropy has since fixed.

.. _iers-Older_astropy_versions:

Older astropy versions
======================

astropy versions before 7.2.3, and 8.0.0, don't use the bundled Bulletin A
by default. Instead, they download it from the IERS the first time a
calculation needs it in each session. If the IERS servers can't be reached,
these versions:

* try again (and time out again) in every later calculation, which can make
  astroplan very slow;
* warn about each failed download;
* fall back to the bundled Bulletin B (astropy 6.x) or the bundled Bulletin A
  (astropy 7.x), so results for recent or future times are less precise and
  come with more warnings.

The best fix is to upgrade astropy. Otherwise, with astropy 7.x, setting
``iers.conf.auto_download = False`` as above makes astropy use the bundled
Bulletin A. With astropy 6.x, that falls back to Bulletin B instead, so use
the bundled Bulletin A directly::

    from astropy.utils import iers
    iers.earth_orientation_table.set(iers.IERS_A.open(iers.IERS_A_FILE))

Either way, keep ``astropy-iers-data`` up to date, since the bundled tables
are then the only ones used.
