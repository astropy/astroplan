#!/usr/bin/python

# ======================================
# Plots for Observation Planning Toolbox
# ======================================

# This is a strawman API for plotting with Astroplan.

# Must first construct an Observer object using astropy.coordinates.

from astroplan import Observer
from astropy.coordinates import EarthLocation
import astropy.units as u
from pytz import timezone

my_location = EarthLocation.from_geodetic(
    '-155d28m48.900s',
    '+19d49m42.600s',
    4163 * u.m,
    )

obs = Observer(
    name='Subaru Telescope',
    location=my_location,
    pressure=0.615 * u.bar,
    relative_humidity=0.11,
    temperature=0 * u.deg_C,
    timezone=timezone('US/Hawaii'),
    description="Subaru Telescope on Mauna Kea, Hawaii"
    )

# Or call up one of of the sites in Astroplan's database.

from astroplan import sites

obs = sites.Keck1

# Then, define your target.

from astroplan import FixedTarget
from astroplan.coordinates import SkyCoord

my_target = FixedTarget(
    name='Polaris',
    coord=SkyCoord('02h31m49.09s', '+89d15m50.8s', frame='icrs')
    )

# You can also define an array with multiple targets.
targets = [
    FixedTarget(name='Polaris',
                coord=SkyCoord('02h31m49.09s', '+89d15m50.8s', frame='icrs')),
    FixedTarget(name='Sirius',
                coord=SkyCoord('06h45m08.9173s', '-16d42m58.017s', frame='icrs')),
    FixedTarget(name='Other Target',
                coord=SkyCoord('01h30m0s', '+10d15m20s', frame='icrs')),
    ]

# ...once you've constructed your Observer and Target objects, you can plot any
#   any number of quantities using the plotting functions and Time objects.


# ===============
# Airmass vs Time
# ===============

from astroplan import plot_airmass
from astropy.time import Time
import numpy as np

# plot_airmass takes an individual Target, an individual Observer,
#   and a Time object as arguments.

# Construct the Time object with as many individual times as you want data
#   points in the resulting airmass plot.

# One way to do this is with a TimeDelta object and numpy.linspace().

start_time = Time('2015-01-15T21:00:00')
end_time = Time('2015-01-16T04:00:00')
delta_t_airmass = end_time - start_time
times_airmass = start_time + delta_t_airmass*np.linspace(0.0, 1.0, 100)

plot_airmass(my_target, observer=obs, time=times_airmass)

# If you want to populate a Time object starting with just one time,
# you can use astropy.units.

times_airmass = start_time + np.linspace(-8, 8, 100)*u.hour

# If you want to plot airmass for multiple targets on the same plot,
#   use the new_plot=False option for each target after the first
#   (the default is True).

plot_airmass(targets[0], observer=obs, time=times_airmass)
plot_airmass(targets[1], observer=obs, time=times_airmass, new_plot=False)
plot_airmass(targets[2], observer=obs, time=times_airmass, new_plot=False)

# If you want to plot airmass for multiple targets on different plots,
#   simply issue the command as normal.

for i in range(0, len(targets)):
    plot_airmass(targets[i], observer=obs, time=times_airmass)

# You can also set a dark background with the dark_plot=True option.

plot_airmass(my_target, observer=obs, time=times_airmass, dark_plot=True)

# For custom line styles/colors, use the style_kwargs option in the following
#   manner.

plot_airmass(
    my_target,
    observer=obs,
    time=times_airmass,
    style_kwarg={'linestyle': '--', 'color': 'r'},
    )


# =========================
# Parallactic Angle vs Time
# =========================

from astroplan import plot_pang

# The plot_pang takes an individual Target, an individual Observer,
#   and a Time datetime object.

# Just as with plot_airmass, construct the Time object with as many
#   individual times as you want data points in the resulting plot.

start_time = Time('2015-01-15T21:00:00')
end_time = Time('2015-01-16T04:00:00')
delta_t_pang = end_time - start_time
times_pang = start_time + delta_t_pang*np.linspace(0.0, 1.0, 50)

plot_pang(my_target, observer=obs, time=times_pang)

# If you want to plot airmass for multiple targets on the same plot,
#   use the new_plot=False option for each target after the first
#   (the default is True).

plot_pang(targets[0], observer=obs, time=times_pang)
plot_pang(targets[1], observer=obs, time=times_pang, new_plot=False)
plot_pang(targets[2], observer=obs, time=times_pang, new_plot=False)

# If you want to plot airmass for multiple targets on different plots,
#   simply issue the command as normal.

for i in range(0, len(targets)):
    plot_pang(targets[i], observer=obs, time=times_pang)

# You can also set a dark background with the dark_plot=True option.

plot_pang(my_target, observer=obs, time=times_pang, dark_plot=True)

# For custom line styles/colors, use the style_kwargs option in the following
#   manner.

plot_pang(
    my_target,
    observer=obs,
    time=times_pang,
    style_kwarg={'linestyle': '--', 'color': 'r'},
    )


# =========
# Sky Chart
# =========

from astroplan import plot_sky

# plot_sky produces a sky chart relative to the observer, showing the
#   locations of Targets at user-specified times.

# plot_sky takes a Target, an Observer and a Time object.
# The default uses alt/az coordinates, is centered on the zenith
#   of the Observer at the datetime given, and has a light background.

# Most users will want to see the location of multiple Targets at a
#   given time.

some_time = Time('2015-01-15T22:00:00')

plot_sky(targets[0], observer=obs, time=some_time)
plot_sky(targets[1], observer=obs, time=some_time)
plot_sky(targets[2], observer=obs, time=some_time)

# But you can also plot guide objects as well, and set them apart by
#   specifying a a different plotting style.

guides = [
    FixedTarget(name='Guide 1',
                coord=SkyCoord('09h10m11.12s', '-50d40m30.20s', frame='icrs')),
    FixedTarget(name='Guide 2',
                coord=SkyCoord('06h45m00.0s', '-16d45m60.0s', frame='icrs')),
    FixedTarget(name='Guide 3',
                coord=SkyCoord('01h30m0s', '+10d15m20s', frame='icrs')),
    ]

for i in range(0, len(guides)):
    plot_sky(
        guides[i],
        observer=obs,
        time=some_time,
        style_kwargs={'marker': '*', 'color': 'k'},
        )

# You can also plot any number of objects on the same plot over a period
#   of time to see how they move over the course of a night.

# Just as with plot_airmass and plot_pang, construct the Time object
#   with as many individual times as you want data points in the resulting
#   plot.

start_time = Time('2015-01-15T21:00:00')
times_sky = start_time + np.linspace(-4.0, 4.0, 9)

for i in range(0, len(targets)):
    plot_sky(targets[i], observer=obs, time=times_sky)

# If you want to plot your objects on separate plots for each time,
#   (e.g., for making an animation) set new_plot=True
#   (the default is false).

for i in range(0, len(times_sky)):
    plot_sky(targets[0], observer=obs, time=times_sky[i], new_plot=True)
    for j in range(1, len(targets)):
        plot_sky(targets[j], observer=obs, time=times_sky[i])
    for k in range(0, len(guides)):
        plot_sky(guides[k], observer=obs, time=times_sky[i])

# You can also set a dark background with the dark_plot=True option.

plot_sky(targets[0], observer=obs, time=start_time, dark_plot=True)


# ===========
# Other plots
# ===========

# Will be added as come up in discussion.
