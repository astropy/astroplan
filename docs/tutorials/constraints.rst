:orphan:

.. doctest-skip-all

******************************
Defining Observing Constraints
******************************

Frequently, we have a long list of targets that we want to observe, and we need
to know which ones are observable given a set of constraints imposed on our
observations by a wide range of limitations. For example, your telescope may
only point over a limited range of altitudes, your targets are only useful
in a range of airmasses, and they must be separated from the moon by some
large angle. The `astroplan.constraints` module is here to help!

Say we're planning to observe from Subaru Observatory in Hawaii on August 1,
2015 from 06:00-12:00 UTC. First, let's set up an `astroplan.Observer` object::

    from astroplan import Observer, FixedTarget
    from astropy.time import Time
    subaru = Observer.at_site("Subaru")
    time_range = Time(["2015-08-01 06:00", "2015-08-01 12:00"])

We're keeping a list of targets in a text file called `targets.txt`, which looks
like this::

    # name ra_degrees dec_degrees
    Polaris 37.95456067 89.26410897
    Vega 279.234734787 38.783688956
    Albireo 292.68033548 27.959680072
    Algol 47.042218553 40.955646675
    Rigel 78.634467067 -8.201638365
    Regulus 152.092962438 11.967208776

We'll read in this list of targets using `astropy.io.ascii`, and create a list
of `astroplan.FixedTarget` objects out of them::

    # Read in the table of targets
    from astropy.io import ascii
    target_table = ascii.read('targets.txt')

    # Create astroplan.FixedTarget objects for each one in the table
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    targets = [FixedTarget(coord=SkyCoord(ra=ra*u.deg, dec=dec*u.deg), name=name)
               for name, ra, dec in target_table]

We will build a bulleted list of our constraints first, then implement them in
code below.

* Our observations with Subaru can only occur between altitudes of ~10-80
  degrees, which we can define using the
  `astroplan.constraints.AltitudeConstraint` class.

* We place an upper limit on the airmass of each target during observations
  using the `astroplan.constraints.AirmassConstraint` class.

* Since we're optical observers, we only want to observe targets at night, so
  we'll also call the `astroplan.constraints.AtNightConstraint` class. We're
  not terribly worried about sky brightness for these bright stars, so we'll
  define "night" times as those between civil twilights by using the class
  method `AtNightConstraint.twilight_civil`::

    from astroplan import (AltitudeConstraint, AirmassConstraint,
                           AtNightConstraint)
    constraints = [AltitudeConstraint(10*u.deg, 80*u.deg),
                   AirmassConstraint(5), AtNightConstraint.twilight_civil()]

This list of constraints can now be applied to our target list to determine
whether:

* the targets are observable given the constraints at *any* times in the time
  range, using `astroplan.constraints.is_observable`,

* the targets are observable given the constraints at *all* times in the time
  range, using `astroplan.constraints.is_always_observable`::

    from astroplan import is_observable, is_always_observable
    # Are targets *ever* observable in the time range?
    ever_observable = is_observable(constraints, time_range, targets, subaru)

    # Are targets *always* observable in the time range?
    always_observable = is_always_observable(constraints, time_range, targets, subaru)

These two functions will return boolean arrays which tell you whether or not
each target is observable given your constraints. Let's print these results in
tabular form:

    >>> from astropy.table import Table
    >>> import numpy as np
    >>> observability_array = np.array([[t.name, ever, always]  for t, ever, always in
    ...                                 zip(targets, ever_observable, always_observable)])
    >>> observability_table = Table(observability_array,
    ...                             names=('Target', 'Ever Observable',
    ...                                    'Always Observable'))
    >>> print(observability_table)
    <Table length=6>
     Target  Ever Observable Always Observable
    unicode7     unicode7         unicode7
    -------- --------------- -----------------
     Polaris            True              True
        Vega            True              True
     Albireo            True             False
       Algol            True             False
       Rigel           False             False
     Regulus           False             False

Now we can see which targets are observable!
