from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import datetime as dt
import pytest

import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, get_sun

from ..moon import get_moon
from ..observer import Observer
from ..target import FixedTarget
from ..constraints import (AltitudeConstraint, AirmassConstraint, AtNightConstraint,
                           is_observable, is_always_observable, observability_table,
                           time_grid_from_range, SunSeparationConstraint,
                           MoonSeparationConstraint, MoonIlluminationConstraint,
                           LocalTimeConstraint)

try:
    import ephem
    HAS_PYEPHEM = True
except ImportError:
    HAS_PYEPHEM = False

vega = FixedTarget(coord=SkyCoord(ra=279.23473479*u.deg, dec=38.78368896*u.deg),
                   name="Vega")
rigel = FixedTarget(coord=SkyCoord(ra=78.63446707*u.deg, dec=8.20163837*u.deg),
                    name="Rigel")
polaris = FixedTarget(coord=SkyCoord(ra=37.95456067*u.deg,
                                     dec=89.26410897*u.deg), name="Polaris")

def test_at_night_basic():
    subaru = Observer.at_site("Subaru")
    time_ranges = [Time(['2001-02-03 04:05:06', '2001-02-04 04:05:06']), # 1 day
                   Time(['2007-08-09 10:11:12', '2007-08-09 11:11:12'])] # 1 hr
    targets = [vega, rigel, polaris]
    for time_range in time_ranges:
        # Calculate constraint using methods on astroplan.Observer:
        observer_is_night = [subaru.is_night(time)
                             for time in time_grid_from_range(time_range)]
        observer_is_night_any = any(observer_is_night)
        observer_is_night_all = all(observer_is_night)

        assert all(is_observable(AtNightConstraint(), subaru, targets,
                                 time_range=time_range) ==
                   len(targets)*[observer_is_night_any])

        assert all(is_always_observable(AtNightConstraint(), subaru, targets,
                                        time_range=time_range) ==
                   len(targets)*[observer_is_night_all])


def test_observability_table():
    subaru = Observer.at_site("Subaru")
    time_ranges = [Time(['2001-02-03 04:05:06', '2001-02-04 04:05:06']), # 1 day
                   Time(['2007-08-09 10:11:12', '2007-08-09 11:11:12'])] # 1 hr
    targets = [vega, rigel, polaris]

    time_range = Time(['2001-02-03 04:05:06', '2001-02-04 04:05:06'])
    constraints = [AtNightConstraint(), AirmassConstraint(3)]

    obstab = observability_table(constraints, subaru, targets,
                                 time_range=time_range)

    assert len(obstab) == 3
    assert np.all(obstab['target name'] == ['Vega', 'Rigel', 'Polaris'])

    assert 'times' in obstab.meta
    assert 'observer' in obstab.meta
    assert 'constraints' in obstab.meta
    assert len(obstab.meta['constraints']) == 2

    np.testing.assert_allclose(obstab['fraction of time observable'],
                               np.array([21, 22, 15])/48)

    #now compare to is_observable and is_always_observable
    is_obs = is_observable(constraints, subaru, targets, time_range=time_range)
    np.testing.assert_allclose(obstab['ever observable'], is_obs)
    all_obs = is_always_observable(constraints, subaru, targets,
                                   time_range=time_range)
    np.testing.assert_allclose(obstab['always observable'], all_obs)


def test_compare_altitude_constraint_and_observer():
    time = Time('2001-02-03 04:05:06')
    time_ranges = [Time([time, time+1*u.hour]) + offset
                   for offset in np.arange(0, 400, 100)*u.day]
    for time_range in time_ranges:
        subaru = Observer.at_site("Subaru")
        targets = [vega, rigel, polaris]

        min_alt = 40*u.deg
        max_alt = 80*u.deg
        # Check if each target meets altitude constraints using Observer
        always_from_observer = [all([min_alt < subaru.altaz(time, target).alt < max_alt
                                     for time in time_grid_from_range(time_range)])
                                for target in targets]
        # Check if each target meets altitude constraints using
        # is_always_observable and AltitudeConstraint
        always_from_constraint = is_always_observable(AltitudeConstraint(min_alt,
                                                                         max_alt),
                                                      subaru, targets,
                                                      time_range=time_range)
        assert all(always_from_observer == always_from_constraint)


def test_compare_airmass_constraint_and_observer():
    time = Time('2001-02-03 04:05:06')
    time_ranges = [Time([time, time+1*u.hour]) + offset
                   for offset in np.arange(0, 400, 100)*u.day]
    for time_range in time_ranges:
        subaru = Observer.at_site("Subaru")
        targets = [vega, rigel, polaris]

        max_airmass = 2
        # Check if each target meets airmass constraint in using Observer
        always_from_observer = [all([subaru.altaz(time, target).secz < max_airmass
                                     for time in time_grid_from_range(time_range)])
                                for target in targets]
        # Check if each target meets altitude constraints using
        # is_always_observable and AirmassConstraint
        always_from_constraint = is_always_observable(AirmassConstraint(max_airmass),
                                                      subaru, targets,
                                                      time_range=time_range)

        assert all(always_from_observer == always_from_constraint)


def test_sun_separation():
    time = Time('2003-04-05 06:07:08')
    apo = Observer.at_site("APO")
    sun = get_sun(time)
    one_deg_away = SkyCoord(ra=sun.ra, dec=sun.dec+1*u.deg)
    five_deg_away = SkyCoord(ra=sun.ra+5*u.deg, dec=sun.dec)
    twenty_deg_away = SkyCoord(ra=sun.ra+20*u.deg, dec=sun.dec)

    constraint = SunSeparationConstraint(min=2*u.deg, max=10*u.deg)
    is_constraint_met = constraint(apo, [one_deg_away, five_deg_away,
                                         twenty_deg_away], times=time)
    assert np.all(is_constraint_met == [[False], [True], [False]])

    constraint = SunSeparationConstraint(max=10*u.deg)
    is_constraint_met = constraint(apo, [one_deg_away, five_deg_away,
                                         twenty_deg_away], times=time)
    assert np.all(is_constraint_met == [[True], [True], [False]])

    constraint = SunSeparationConstraint(min=2*u.deg)
    is_constraint_met = constraint(apo, [one_deg_away, five_deg_away,
                                         twenty_deg_away], times=time)
    assert np.all(is_constraint_met == [[False], [True], [True]])


@pytest.mark.skipif('not HAS_PYEPHEM')
def test_moon_separation():
    time = Time('2003-04-05 06:07:08')
    apo = Observer.at_site("APO")
    altaz_frame = apo.altaz(time)
    moon = get_moon(time, apo.location, apo.pressure)
    one_deg_away = SkyCoord(az=moon.az, alt=moon.alt+1*u.deg, frame=altaz_frame)
    five_deg_away = SkyCoord(az=moon.az+5*u.deg, alt=moon.alt,
                             frame=altaz_frame)
    twenty_deg_away = SkyCoord(az=moon.az+20*u.deg, alt=moon.alt,
                               frame=altaz_frame)

    constraint = MoonSeparationConstraint(min=2*u.deg, max=10*u.deg)
    is_constraint_met = constraint(apo, [one_deg_away, five_deg_away,
                                         twenty_deg_away], times=time)
    assert np.all(is_constraint_met == [[False], [True], [False]])


    constraint = MoonSeparationConstraint(max=10*u.deg)
    is_constraint_met = constraint(apo, [one_deg_away, five_deg_away,
                                         twenty_deg_away], times=time)
    assert np.all(is_constraint_met == [[True], [True], [False]])

    constraint = MoonSeparationConstraint(min=2*u.deg)
    is_constraint_met = constraint(apo, [one_deg_away, five_deg_away,
                                         twenty_deg_away], times=time)
    assert np.all(is_constraint_met == [[False], [True], [True]])


@pytest.mark.skipif('not HAS_PYEPHEM')
def test_moon_illumination():
    times = Time(["2015-08-29 18:35", "2015-09-05 18:35", "2015-09-15 18:35"])
    lco = Observer.at_site("LCO")
    # At these times, moon illuminations are:
    # [ 0.99946328  0.46867661  0.05379006]

    constraint = MoonIlluminationConstraint(min=0.2, max=0.8)
    is_constraint_met = [constraint(lco, None, times=time) for time in times]
    assert np.all(is_constraint_met == [False, True, False])

    constraint = MoonIlluminationConstraint(min=0.2)
    is_constraint_met = [constraint(lco, None, times=time) for time in times]
    assert np.all(is_constraint_met == [True, True, False])

    constraint = MoonIlluminationConstraint(max=0.8)
    is_constraint_met = [constraint(lco, None, times=time) for time in times]
    assert np.all(is_constraint_met == [False, True, True])


def test_local_time_constraint_utc():
    time = Time('2001-02-03 04:05:06')
    subaru = Observer.at_site("Subaru")
    constraint = LocalTimeConstraint(min=dt.time(23,50), max=dt.time(4,8))
    is_constraint_met = constraint(subaru, None, times=time)
    assert is_constraint_met == [True]

    constraint = LocalTimeConstraint(min=dt.time(0,2), max=dt.time(4,3))
    is_constraint_met = constraint(subaru, None, times=time)
    assert is_constraint_met == [False]

    constraint = LocalTimeConstraint(min=dt.time(3,8), max=dt.time(5,35))
    is_constraint_met = constraint(subaru, None, times=time)
    assert is_constraint_met == [True]


def test_local_time_constraint_hawaii_tz():
    # Define timezone in Observer.timezone
    time = Time('2001-02-03 04:05:06')
    subaru = Observer.at_site("Subaru", timezone="US/Hawaii")
    constraint = LocalTimeConstraint(min=dt.time(23,50), max=dt.time(4,8))
    is_constraint_met = constraint(subaru, None, times=time)
    assert is_constraint_met == [True]

    constraint = LocalTimeConstraint(min=dt.time(0,2), max=dt.time(4,3))
    is_constraint_met = constraint(subaru, None, times=time)
    assert is_constraint_met == [False]

    constraint = LocalTimeConstraint(min=dt.time(3,8), max=dt.time(5,35))
    is_constraint_met = constraint(subaru, None, times=time)
    assert is_constraint_met == [True]


def test_docs_example():
    # Test the example in astroplan/docs/tutorials/constraints.rst
    target_table_string = """# name ra_degrees dec_degrees
    Polaris 37.95456067 89.26410897
    Vega 279.234734787 38.783688956
    Albireo 292.68033548 27.959680072
    Algol 47.042218553 40.955646675
    Rigel 78.634467067 -8.201638365
    Regulus 152.092962438 11.967208776"""

    from astroplan import Observer, FixedTarget
    from astropy.time import Time
    subaru = Observer.at_site("Subaru")
    time_range = Time(["2015-08-01 06:00", "2015-08-01 12:00"])

    # Read in the table of targets
    from astropy.io import ascii
    target_table = ascii.read(target_table_string)

    # Create astroplan.FixedTarget objects for each one in the table
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    targets = [FixedTarget(coord=SkyCoord(ra=ra*u.deg, dec=dec*u.deg), name=name)
               for name, ra, dec in target_table]

    from astroplan import Constraint, is_observable
    from astropy.coordinates import Angle

    class VegaSeparationConstraint(Constraint):
        """
        Constraint the separation from Vega
        """
        def __init__(self, min=None, max=None):
            """
            min : `~astropy.units.Quantity` or `None` (optional)
                Minimum acceptable separation between Vega and target. `None`
                indicates no limit.
            max : `~astropy.units.Quantity` or `None` (optional)
                Minimum acceptable separation between Vega and target. `None`
                indicates no limit.
            """
            self.min = min
            self.max = max

        def compute_constraint(self, times, observer, targets):

            # Vega's coordinate must be non-scalar for the dimensions
            # to work out properly when combined with other constraints which
            # test multiple times
            vega = SkyCoord(ra=[279.23473479]*u.deg, dec=[38.78368896]*u.deg)

            # Calculate separation between target and vega
            vega_separation = Angle([vega.separation(target.coord)
                                     for target in targets])

            # If a maximum is specified but no minimum
            if self.min is None and self.max is not None:
                mask = vega_separation < self.max

            # If a minimum is specified but no maximum
            elif self.max is None and self.min is not None:
                mask = self.min < vega_separation

            # If both a minimum and a maximum are specified
            elif self.min is not None and self.max is not None:
                mask = ((self.min < vega_separation) & (vega_separation < self.max))

            # Otherwise, raise an error
            else:
                raise ValueError("No max and/or min specified in "
                                 "VegaSeparationConstraint.")

            # Return an array that is True where the target is observable and
            # False where it is not
            return mask

    constraints = [VegaSeparationConstraint(min=5*u.deg, max=30*u.deg)]
    observability = is_observable(constraints, subaru, targets,
                                  time_range=time_range)

    assert all(observability == [False, False, True, False, False, False])
