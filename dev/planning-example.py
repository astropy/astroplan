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
from astropy.coordinates import EarthLocation
import pytz

# Define the observer with an instance of astropy.coordinates.EarthLocation,
# also use pytz.timezone() argument directly as `timezone` keyword's input
longitude = '-155d28m48.900s'
latitude = '+19d49m42.600s'
elevation = 4163 * u.m
location = EarthLocation.from_geodetic(longitude, latitude, elevation)

obs = Observer(name='Subaru Telescope',
               location=location,
               pressure=0.615 * u.bar,
               relative_humidity=0.11,
               temperature=0 * u.deg_C,
               timezone=pytz.timezone('US/Hawaii'),
               description="Subaru Telescope on Mauna Kea, Hawaii")

# It would also be desirable to be able to have a small database of
# common telescopes.  Maybe this can at first simply take the form of
# a python module:
from astroplan import sites

obs = sites.Keck1

# Environmental conditions should be updatable.
obs.pressure = 0.600 * u.bar
obs.relative_humidity = 0.2
obs.temperature = 10 * u.deg_C

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

# sun_set
obs.sun_set(time_obs, which='nearest') # Default
obs.sun_set(time_obs, which='next')
obs.sun_set(time_obs, which='previous')
obs.sun_set() # Will be use time=Time.now(), which='nearest'

# sun_rise
obs.sun_rise(time_obs, which='nearest') # Default
obs.sun_rise(time_obs, which='next')
obs.sun_rise(time_obs, which='previous')

# moon rise
obs.moon_rise(time_obs, which='nearest')

# moon set
obs.moon_set(time_obs, which='nearest')

# The above functions can be called with an `angle` keyword argument to specify
# a particular horizon angle for rising or setting, or can be called with
# convenience functions for particular morning/evening twilight.
# For example, to compute astronomical twilight by specifying the `angle`:
obs.sun_set(time_obs, which='next', angle=18*u.degree)

# evening (astronomical) twilight
obs.evening_astronomical(time_obs)

# evening (nautical) twilight
obs.evening_nautical(time_obs)

# evening (civil) twilight
obs.evening_civil(time_obs)

# morning (nautical) twilight
obs.morning_nautical(time_obs)

# morning (civil) twilight
obs.morning_civil(time_obs)

# morning (astronomical) twilight
obs.morning_astronomical(time_obs)

# what is the moon illumination?
# returns a float, which is percentage of the moon illuminated
obs.moon_illumination(time_obs)

# what is the moon altitude and azimuth?
obs.moon_altaz(time_obs)

# Other sun-related convenience functions:
obs.noon(time_obs, which='nearest')
obs.midnight(time_obs, which='nearest')

# ==============
# Target objects
# ==============
'''
A target object defines a single observation target, with coordinates
and an optional name.
'''
from astroplan import FixedTarget

# Define a target.

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
from astroplan import TimeRange, AltitudeRange, AirmassRange, is_observable
# `is_observable` is a temporary function which will eventually be a method of
# something to support caching

# Times in TimeRange can be passed in as strings, will be passed to the
# Time constructor.
constraint_list = [TimeRange("2015-05-01 18:30", "2015-05-02 05:30"),
                   AirmassRange(1.2), AltitudeRange(15.0*u.deg, 89.0*u.deg)]

# (AboveAirmass will be a subclass of AltitudeWindow)

# Combine a list of constraints to run on Observer, FixedTarget, and time to
# determine the observability of target
constraints = is_observable(constraint_list, obs, t1, time_obs)

# Test only a single constraint:
constraints = is_observable(AirmassRange(1.2), obs, t1, time_obs)

# AirmassRange can accept two bounding airmasses, assumes single argument is
# an upper limit, lower limit = 1.

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
sf = FixedTarget(SkyCoord('09d40m00.00s', '43d00m00.00s'), name='Sf')
sm = FixedTarget(SkyCoord('10d30m00.00s', '36d00m00.00s'), name='Sm')

# Coordinate arithmetic gives separations in RA, Dec, alt, az
dra, ddec = sf.ra - sm.ra, sf.dec - sm.dec
dalt = obs.altaz(time_obs, sf).alt - obs.altaz(time_obs, sm).alt
dazt = obs.altaz(time_obs, sf).az - obs.altaz(time_obs, sm).az
