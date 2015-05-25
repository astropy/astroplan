#!/usr/bin/python
# -*- coding: utf-8 -*-

# ============================
# Observation Planning Toolbox
# ============================
"""
The observation planning toolbox consists of a number of entities that
simplify the planning and scheduling of observations by abstracting much
of the boilerplate code necessary to accomplish the same thing using
core astropy or ephemerides objects.  These can be thought of as common
idioms or convenience classes and functions that would be recreated over
and over by astronomers and or software engineers building telescope
scheduling software.

Important assumptions:

- It should be convenient to specify dates and times in the local timezone
  of the observer.

- It should be convenient to do timezone conversions easily.

  This means methods that can take dates can either take astropy
  dates (as used internally) or python datetime objects (not naive, but
  tagged with timezone).  So dates can be passed in as timezone-aware
  datetime objects and results can be easily converted between such.

"""

# ================
# Observer objects
# ================
"""
An observer defines a single terrestrial location and applicable
environmental attributes thereof that affect observation calculations
at that location.
"""
# Define an observer by full specification.
# Latitude/Longitude take input as strings or quantities,
# elevation, pressure, temperature must be quantities.
# Accept `timezone` from pytz
from astroplan import Observer
import astropy.units as u
import pytz

obs = Observer(name='Subaru Telescope',
               longitude='-155:28:48.900',
               latitude='+19:49:42.600',
               elevation=4163 * u.meter,
               pressure=0.615 * u.bar,
               relative_humidity=0.11,
               temperature=0 * u.deg_C,
               timezone=pytz.timezone('US/Hawaii'),
               description="Subaru Telescope on Mauna Kea, Hawaii")

# Define the same observer with an instance of astropy.coordinates.EarthLocation
from astropy.coordinates import EarthLocation
obs = Observer(name='Subaru',
               location=EarthLocation(lon='-155:28:48.900',
                                      lat='+19:49:42.600'),
               elevation=4163 * u.meter,
               pressure=0.615 * u.bar,
               relative_humidity=0.11,
               temperature=0 * u.deg_C,
               timezone=pytz.timezone('US/Hawaii'),
               description="Subaru Telescope on Mauna Kea, Hawaii")

# longitude and latitude would be represented internally by astropy
# angles.Longitude/Latitude.  So any value that can naturally initialize
# one of those objects could be passed as the parameter to Observer

# It's also possible to initialize using astropy objects and other
# direct entities
from astropy.coordinates.angles import Longitude, Latitude

obs = Observer(name='Keck 1',
               longitude=Longitude(155.4750 * u.degree),
               latitude=Latitude(19.8264 * u.degree),
               elevation=4160 * u.m,
               pressure=0.615 * u.bar,
               relative_humidity=0.11,
               temperature=0.0 * u.deg_C,
               timezone=pytz.timezone('US/Hawaii'),
               description="W.M. Keck 1 Telescope on Mauna Kea, Hawaii")

# It would also be desirable to be able to have a small database of
# common telescopes.  Maybe this can at first simply take the form of
# a python module:
from astroplan import sites

obs = sites.Keck1

# Environmental conditions should be updatable.  Suggest using keyword
# parameters (this would make it easy to update an internal AltAz coordinate):

obs.set_environment(pressure=0.600 * u.bar, relative_humidity=0.2,
                    temperature=10 * u.deg_C)

# It's easy to construct date/times from the observer, since it knows about
# its timezone.  These dates can be passed directly to target check methods.

# dates default to midnight if not otherwise specified
# Date can be anything that astropy.time.Time will accept
day = obs.get_date("2015-05-01")

# Can get more and more specific
t1 = obs.get_date("2015-05-01 21:00")
t1 = obs.get_date("2015-05-01 21:33:41")

# default is current time
cur_time = obs.get_date()

# Dates are assumed to be in the observer's time zone unless a time zone
# specifier is given
t2_utc = obs.get_date("2015-05-01 21:33:41 UTC")

"""
We can ask about typical times of interest at this observing position
Returned dates are assumed to be in the observer's time zone for convenience.
The date parameter, if given, can be an astropy Time or a datetime instance.
If no parameter is given, the current date/time is assumed.
"""

# sunset
obs.next_sunset(date=day)
obs.previous_sunset(date=day)

# sunrise
obs.next_sunrise(date=day)
obs.previous_sunrise(date=day)

# moon rise
obs.next_moon_rise(date=day)
obs.previous_moon_rise(date=day)

# moon set
obs.next_moon_set(date=day)
obs.previous_moon_set(date=day)

# The above functions can be called with an `angle` keyword argument to specify
# a particular horizon angle for rising or setting, or can be called with
# convenience functions for particular morning/evening twilight:

obs.next_sunset(date=day, angle=18*u.degree) # (i.e., astronomical twilight)

# evening (astronomical) twilight
obs.evening_astronomical(date=day)

# evening (nautical) twilight
obs.evening_nautical(date=day)

# evening (civil) twilight
obs.evening_civil(date=day)

# morning (nautical) twilight
obs.morning_nautical(date=day)

# morning (civil) twilight
obs.morning_civil(date=day)

# morning (astronomical) twilight
obs.morning_astronomical(date=day)

# what is the moon illumination?
# returns a float, which is percentage of the moon illuminated
obs.moon_illumination(date=t1)

# what is the moon altitude and azimuth?
# returns an altaz coordinate
obs.moon_position(date=t1)

# Other sun-related convenience functions:
obs.noon(date=day)
obs.midnight(date=day)

# ==============
# Target objects
# ==============
'''
A target object defines a single observation target, with coordinates
and an optional name.
'''
from astroplan import FixedTarget, airmass_plot

# Define a target.
t1 = FixedTarget(name='Polaris', ra='02:31:49.09', dec='+89:15:50.8')

# initializing from astropy entities directly
from astropy.coordinates import SkyCoord

t1 = FixedTarget(name='Polaris',
               coord=SkyCoord('02h31m49.09s', '+89d15m50.8s', frame='icrs'))

# Leaves scope for NonSiderealTarget, etc.

# Convenience methods for looking up/constructing targets by name via
# astroquery?

# Plot the airmass for the specified target:
airmass_plot(t1, observer=obs, datetime=day)

# ================================
# Condition objects, observability
# ================================

"""
Q: Can I observe this target on the night of May 1, 2015 between 18:30
and 05:30 HST for 30 minutes somewhere between 15 degrees and 89 degrees
altitude?
"""

# first we define the conditions
from astroplan import Constraints, TimeWindow, AltitudeRange, AboveAirmass

# times can be passed in as strings (interpreted as for get_date()) or
# as astropy Time or datetime objects
cons = [TimeWindow("2015-05-01 18:30", "2015-05-02 05:30"),
        AltitudeRange(15.0*u.deg, 89.0*u.deg), AboveAirmass(1.2)]
# define a target
tgt = FixedTarget(name='S5', ra='14:20:00.00', dec='48:00:00.00')

# returns an object with information about the observability with these
# constraints
info = obs.Constraints(tgt, cons)

'''
attributes of return object:
  observable: a boolean -- whether target is observable at desired constraints
  rise_time: a datetime -- first time target is visible in at these conditions
  set_time: a datetime -- time target is no longer visible at these conditions

Use of a return object means that additional items can be added in the
future without breaking the API (e.g. subclasses are free to attach
additional items).  It can be also useful to attach an "explanation" for why
something cannot be observed.
'''

'''
Suggest that Constraints implement a time range, an altitude range and an
airmass limit at beginning.  Other things that could be added later:
 - azimuth range
 - airmass range

Idea is that Constraints can be subclassed to add additional constraints.
'''

# ======================================================
# Other useful calculations wrt an observer and a target
#=======================================================

# calculate the distance in alt and az degrees between two targets at
# the given time (e.g. to calculate slew time)
sf = FixedTarget(name='Sf', ra='09:40:00.00', dec='43:00:00.00')
sm = FixedTarget(name='Sm', ra='10:30:00.00', dec='36:00:00.00')
t = obs.get_date("2015-05-01 23:00")

# Coordinate arithmetic gives separations in RA, Dec, alt, az
dra, ddec = sf.ra - sm.ra, sf.dec - sm.dec
dalt = obs.altaz(sf, t).alt - obs.altaz(sm, t).alt
dazt = obs.altaz(sf, t).az - obs.altaz(sm, t).az

# tell me about object sm in relation to observer obs at time t
# returns an Observation object
cr = obs.calc(sm, t)

"""
Note use of Observation object--most computations are lazily delayed
until they are requested from the object.  This saves a lot of time if many
objects are evaluated wrt a given time and observation position and we only
need some of the calculation results (e.g. for scheduling).  For example, we
might evaluate 400 targets against a given time slot and we are only
interested in their altitude and airmass--calculation time dominates when
scheduling from a large number of potential targets.

Observation objects contain a number of properties, most of which
represent inter-related calculations that are not computed until necessary.
These are exposed via @property attributes:

cr.airmass -- airmass at observed position
cr.pang -- parallactic angle
cr.ha -- hour angle
cr.ut -- time of calculation result in ust
cr.lt -- time of calculation result in local time
cr.gmst -- time of calculation result in gmst
cr.lmst -- ditto local mean sidereal time
cr.moon_sep -- moon separation from the target
cr.altitude -- altitude of the object
cr.azimuth -- azimuth of the object
"""

# we can also ask for a calculation via the target (same result)
cr = sm.calc(obs, t)
