#!/usr/bin/python

# ======================================
# Plots for Observation Planning Toolbox
# ======================================

# This is a strawman API for plotting with Astroplan.

# Must first construct an Observer object...

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

day = obs.get_date("2015-05-01")

# ...then define your target... 

from astroplan import FixedTarget

t1 = FixedTarget(name='Polaris', ra='02:31:49.09', dec='+89:15:50.8')

# ...once you've constructed your Observer and Target objects, you can plot any
#   any number of quantities.


# ===============
# Airmass vs Time
# ===============

from astroplan import plot_airmass

# The default airmass plot takes a Target, an Observer, and a datetime object
#   as arguments. 
# The default time range is six hours before to six hours after the datetime
#   specified, and a light background.
plot_airmasss(t1, observer=obs, datetime=day)

# You can also pass several targets to plot_airmass once you create them.

target_list = [t1, t2, t3]

plot_airmass(target_list, observer=obs, datetime=day)

# You can set several options, such as change the time start/end (relative
#   to the datetime object), set a dark background, etc., with optional keywords.

plot_airmass(target_list, observer=obs, datetime=day, time_window=(-8, +8), dark)

# You can leave the start or end time to the default simply by leaving that 
#   option blank.

plot_airmass(target_list, observer=obs, datetime=day, time_window=(-8,), dark)


# =========================
# Parallactic Angle vs Time
# =========================

from astroplan import plot_pang

# The default parallactic angle plot takes one or more Target objects, 
#   an Observer, and a datetime object.
# The defaults are: time range is six hours before to six hours after the datetime 
#   specified, a light background.
# You can set options the same way you set them for plot_airmass.

plot_pang(target_list, observer=obs, datetime=day)


# =========
# Sky Chart
# =========

from astroplan import plot_sky

# The default sky chart plot takes one or more Target objects, an Observer, 
#   and a datetime object. 
# The default uses alt/az coordinates, is centered on the zenith 
#   of the Observer at the datetime given, and has a light background.

plot_sky(target_list, observer=obs, datetime=day)

# You can set several options, including the use of RA/Dec coordinates, 
#   a dark background, turn markers on, etc.
# Markers are pulled from the targets provided. 

plot_sky(target_list, observer=obs, datetime=day, radec, dark, markers)

# You can show how your target objects move over time (with respect to your 
#   zenith) in one of two ways:

# 1) Provide multiple datetime objects after creating them:

datetime_list = [time1, time2, time3, time4]
plot_sky(target_list, observer=obs, datetime=datetime_list)

# 2) Provide a time window and number of "snapshots" you want displayed. 
#   Here, 8 hours before to 8 hours after the datetime object, 
#   and 5 snapshots evenly spaced throughout the time window (in addition to 
#   the datetime specified).

plot_sky(target_list, observer=obs, datetime=day, time_window=(-8, +8, 5)


# ===========
# Other plots
# ===========

# Will be added as come up in discussion.

