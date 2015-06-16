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

# Define the same observer with an instance of astropy.coordinates.EarthLocation,
# also use pytz.timezone() argument directly as `timezone` keyword's input
from astropy.coordinates import EarthLocation
obs = Observer(name='Subaru',
               location=EarthLocation(lon='-155:28:48.900',
                                      lat='+19:49:42.600'),
               elevation=4163 * u.meter,
               pressure=0.615 * u.bar,
               relative_humidity=0.11,
               temperature=0 * u.deg_C,
               timezone='US/Hawaii',
               description="Subaru Telescope on Mauna Kea, Hawaii")

# longitude and latitude would be represented internally by astropy
# coordinates.angles.Longitude/Latitude.  So any value that can naturally 
# initialize one of those objects could be passed as the parameter to Observer

# It's also possible to initialize using astropy objects and other
# direct entities
from astropy.coordinates import Longitude, Latitude

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

# Dates are assumed to be UTC by default

"""
We can ask about typical times of interest at this observing position
Returned dates are assumed to be in the observer's time zone for convenience.
The date parameter, if given, can be an astropy Time or a datetime instance.
If no parameter is given, the current date/time is assumed.
"""

# Define a time
from astropy.time import Time
time_obs = Time(2457189.500000, format='jd')

# sunset
obs.sunset(date=time_obs, which='nearest') # Default
obs.sunset(date=time_obs, which='next')
obs.sunset(date=time_obs, which='previous')

# sunrise
obs.sunrise(date=time_obs, which='nearest') # Default
obs.sunrise(date=time_obs, which='next')
obs.sunrise(date=time_obs, which='previous')

# moon rise
obs.moon_rise(date=time_obs, which='nearest')

# moon set
obs.moon_set(date=time_obs, which='nearest')

# The above functions can be called with an `angle` keyword argument to specify
# a particular horizon angle for rising or setting, or can be called with
# convenience functions for particular morning/evening twilight.
# For example, to compute astronomical twilight by specifying the `angle`:
obs.sunset(date=time_obs, which='next', angle=18*u.degree) 

# evening (astronomical) twilight
obs.evening_astronomical(date=time_obs)

# evening (nautical) twilight
obs.evening_nautical(date=time_obs)

# evening (civil) twilight
obs.evening_civil(date=time_obs)

# morning (nautical) twilight
obs.morning_nautical(date=time_obs)

# morning (civil) twilight
obs.morning_civil(date=time_obs)

# morning (astronomical) twilight
obs.morning_astronomical(date=time_obs)

# what is the moon illumination?
# returns a float, which is percentage of the moon illuminated
obs.moon_illumination(date=time_obs)

# what is the moon altitude and azimuth?
obs.moon_altaz(date=time_obs)

# Other sun-related convenience functions:
obs.noon(date=time_obs, which='nearest')
obs.midnight(date=time_obs, which='nearest')

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
# astroquery:
t1 = FixedTarget.from_name('Polaris')

# ================================
# Condition objects, observability
# ================================

"""
Q: Can I observe this target on the night of May 1, 2015 between 18:30
and 05:30 HST, above airmass 1.2, on a telescope that can observe between 
15 degrees and 89 degrees altitude?
"""

# first we define the conditions
from astroplan import TimeWindow, AltitudeWindow, AboveAirmass, is_observable

# times can be passed in as strings (interpreted as for get_date()) or
# as astropy Time or datetime objects
constraint_list = [TimeWindow("2015-05-01 18:30", "2015-05-02 05:30"),
                   AltitudeWindow(15.0*u.deg, 89.0*u.deg), AboveAirmass(1.2)]
# (AboveAirmass will be a subclass of AltitudeWindow)

# Define a target
tgt = FixedTarget(name='S5', ra='14:20:00.00', dec='48:00:00.00')

# Combine a list of constraints to run on Observer, FixedTarget, and time to 
# determine the observability of target
constraints = is_observable(constraint_list, obs, tgt, time_obs)

# Test only a single constraint:
constraints = is_observable(AboveAirmass(1.2), obs, tgt, time_obs)

# `constraints` will be a boolean where True=observable. For a list of 
# targets, observatories, or times, `constraints` may be a booleans array

# We will eventually need a more complicated method that minimizes a cost 
# function when optimizing an observing schedule given the results of 
# `is_observable`.

# ======================================================
# Other useful calculations wrt an observer and a target
#=======================================================

# calculate the distance in alt and az degrees between two targets at
# the given time (e.g. to calculate slew time)
sf = FixedTarget(name='Sf', ra='09:40:00.00', dec='43:00:00.00')
sm = FixedTarget(name='Sm', ra='10:30:00.00', dec='36:00:00.00')

# Coordinate arithmetic gives separations in RA, Dec, alt, az
dra, ddec = sf.ra - sm.ra, sf.dec - sm.dec
dalt = obs.altaz(sf, time_obs).alt - obs.altaz(sm, time_obs).alt
dazt = obs.altaz(sf, time_obs).az - obs.altaz(sm, time_obs).az
