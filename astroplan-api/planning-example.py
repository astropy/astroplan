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
# Define an observer by full specification
#
# Because there may be many attributes for an observer, I suggest using
# keyword parameters for most of them.
from astroplan.entity import Observer

obs = Observer(name='Subaru Telescope',
               longitude='-155:28:48.900',
               latitude='+19:49:42.600',
               elevation=4163,
               pressure=615,
               relative_humidity=0.11,
               temperature=0,
               timezone='US/Hawaii',
               description="Subaru Telescope on Mauna Kea, Hawaii")

# longitude and latitude would be represented internally by astropy
# angles.Longitude/Latitude.  So any value that can naturally initialize
# one of those objects could be passed as the parameter to Observer
# Other items are represented internally using the appropriate astropy
# Quantities and Units and so can be initialized similarly.

# It's also possible to initialize using astropy objects and other
# direct entities
from astropy.coordinates.angles import Longitude, Latitude
import astropy.units as u
import pytz

obs = entity.Observer(name='Keck 1',
                      longitude=Longitude(155.4750 * u.degree),
                      latitude=Latitude(19.8264 * u.degree)
                      elevation=4160 * u.m,
                      pressure=615 * u.hPa,
                      relative_humidity=0.11,
                      temperature=0.0 * u.C,
                      timezone=pytz.timezone('US/Hawaii'),
                      description="W.M. Keck 1 Telescope on Mauna Kea, Hawaii")

# It would also be desirable to be able to have a small database of
# common telescopes.  Maybe this can at first simply take the form of
# a python module:
from astroplan import sites

obs = sites.get_observer('keck1')

# Environmental conditions should be updatable.  Suggest using keyword
# parameters (this would make it easy to update an internal AltAz coordinate):
#
obs.set_env_conditions(pressure=600 * u.hPa, relative_humidity=0.2,
                       temperature=10 * u.C)

# It's easy to construct date/times from the observer, since it knows about
# its timezone.  These dates can be passed directly to target check methods.

# dates default to midnight if not otherwise specified
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
obs.sunset(date=day)

# sunrise
obs.sunrise(date=day)

# moon rise
obs.moon_rise(date=day)

# moon set
obs.moon_set(date=day)

# evening (nautical) twilight
obs.evening_twilight_12(date=day)

# evening (civil) twilight
obs.evening_twilight_18(date=day)

# morning (nautical) twilight
obs.morning_twilight_12(date=day)

# morning (civil) twilight
obs.morning_twilight_18(date=day)

# what is the moon illumination?
# returns a float, which is percentage of the moon illuminated
obs.moon_illumination(date=t1)

# what is the moon altitude and azimuth?
# returns an altaz coordinate
obs.moon_position(date=t1)

# ==============
# Target objects
# ==============
'''
A target object defines a single observation target, with coordinates
and an optional name.
'''
from astroplan.entity import SiderealTarget

# Define a target.
t1 = SiderealTarget(name='Polaris', ra='02:31:49.09', dec='+89:15:50.8')

# initializing from astropy entities directly
from astropy.coordinates import SkyCoord

t1 = SiderealTarget(name='Polaris',
               coord=SkyCoord('02h31m49.09s', '+89d15m50.8s', frame='icrs'))

# Leaves scope for NonSiderealTarget, etc.

# Convenience methods for looking up/constructing targets by name via
# astroquery?


# ================================
# Condition objects, observability
# ================================

"""
Q: Can I observe this target on the night of May 1, 2015 between 18:30
and 05:30 HST for 30 minutes somewhere between 15 degrees and 89 degrees
altitude?
"""

# first we define the conditions
from astroplan.entity import Constraints

# times can be passed in as strings (interpreted as for get_date()) or
# as astropy Time or datetime objects
cts = Constraints(time_start="2015-05-01 18:30", time_end="2015-05-02 05:30",
                  alt_min=15.0*u.deg, alt_max=89.0*u.deg, duration=30*u.min)

# define a target
tgt = SiderealTarget(name='S5', ra='14:20:00.00', dec='48:00:00.00')

# returns an object with information about the observability with these
# constraints
info = obs.observable(tgt, cts)

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

# same question, but specify the lower elevation bound to be the
# telescope minimum and I will specify an airmass of 1.2 for the target
cond = Constraints(time_start="2015-05-01 18:30", time_end="2015-05-02 05:30",
                   alt_min=15.0*u.deg, alt_max=89.0*u.deg, duration=30*u.min,
                   airmass=1.2)
info = obs.observable(tgt, cond)

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
sf = SiderealTarget(name='Sf', ra='09:40:00.00', dec='43:00:00.00')
sm = SiderealTarget(name='Sm', ra='10:30:00.00', dec='36:00:00.00')
t = obs.get_date("2015-05-01 23:00")

# returns a vector or other astropy quantity specifying a delta in alt, az
obs.distance(sf, sm, t)

# tell me about object sm in relation to observer obs at time t
# returns a CalculationResult object
cr = obs.calc(sm, t)

"""
Note use of CalculationResult object--most computations are lazily delayed
until they are requested from the object.  This saves a lot of time if many
objects are evaluated wrt a given time and observation position and we only
need some of the calculation results (e.g. for scheduling).  For example, we
might evaluate 400 targets against a given time slot and we are only
interested in their altitude and airmass--calculation time dominates when
scheduling from a large number of potential targets.

CalculationResult objects contain a number of properties, most of which
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
