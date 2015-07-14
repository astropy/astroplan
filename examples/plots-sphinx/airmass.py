"""
============
plot_airmass
============

TODO:
    1) Link to modules (astroplan.Target, etc.)

Astroplan Plotting Examples - Airmass vs Time

To generate airmass vs. time plots with Astroplan, use the plot_airmass
function.

This function takes, at minimum, a Target, an Observer and a Time object as
input.
Optional arguments include an Axes object and a style dictionary.

plot_airmass will return an Axes object that contains airmass data for the
window of time specified by the Time object passed in. You can further
manipulate the returned Axes object, or simply print the plot to your display
or save it as a file.
"""

"""
Passing in a Time object with just one time.

Any plot function in Astroplan with a time-based axis will allow you to pass in
a Time object with just one time value in it. If you do this, the plotting
function will use a 24-hour window centered around the time you passed in.
"""
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import EarthLocation, SkyCoord
from pytz import timezone
from astropy.time import Time

from astroplan import Observer
from astroplan import FixedTarget
from astroplan.plots import plot_airmass

# Set up Observer, Target and observation time objects.
longitude = '-155d28m48.900s'
latitude = '+19d49m42.600s'
elevation = 4163 * u.m
location = EarthLocation.from_geodetic(longitude, latitude, elevation)

observer = Observer(name='Subaru Telescope',
                    location=location,
                    pressure=0.615 * u.bar,
                    relative_humidity=0.11,
                    temperature=0 * u.deg_C,
                    timezone=timezone('US/Hawaii'),
                    description="Subaru Telescope on Mauna Kea, Hawaii")

coordinates = SkyCoord('02h31m49.09s', '+89d15m50.8s', frame='icrs')
target = FixedTarget(name='Polaris', coord=coordinates)

observe_time = Time('2015-01-15 23:30:00')

plot_airmass(target, observer, observe_time)
plt.show()
