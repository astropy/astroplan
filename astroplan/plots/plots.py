# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
import random

# Originally authored by Jazmin Berlanga Medina
# (jazmin.berlanga@gmail.com), with input from
# Adrian Price-Whelan (email), Erik Tollerud (email),
# Christoph Deil (email), Eric Jeschke (email) and Brett Morris (email).


# Replace function call later with below.
def plot_airmass(target, observer, time, ax=None, style_kwargs=None, dark_plot=False):
    """
    TODO:
        1) Timezones?
        2) Limit airmass <= 3.
        3) Plot lines instead of points.
        4) Legend should not be transparent.
        5) Airmass equation? Should this be moved to a separate function?
        6) dark_plot option

    Returns an ax object with an airmass vs. time plot.

    If an ax object is passed in, plots an additional airmass plot on top.
    Otherwise, creates a new ax object with an airmass plot.

    Parameters
    ----------
    target: `~astroplan.FixedTarget`
        The celestial body of interest.

    observer: `~astroplan`
        The person, telescope, observatory, etc. doing the observing.

    time : `~astropy.time.Time`
        An collection of times to use in generating data to plot.
        Can be scalar (len(time)=1) or contain multiple times.

    ax : `~matplotlib.axes.Axes` or None
        Optional.

    style_kwargs : WHAT TYPE IS DICTIONARY? or None
        Optional.
        A dictionary with plotting style options.

    dark_plot: Boolean.
        Optional. False (default) or True.
    """

    # Set up plot axes and style if needed.
    if ax is None:
        ax = plt.gca()
    if style_kwargs is None:
        # Need to add more style options for ax's with existing plots?
        style_kwargs = {'linestyle': '-', 'color': 'r'}
    #if dark_plot is True:
        # Set dark plot style stuff here.

    # Populate time window if needed.
    if time.isscalar:
        time = time + np.linspace(-12, 12, 100)*u.hour

    # Calculate airmass
    altitude = observer.altaz(time, target).alt
    airmass = (1.0/np.cos(90*u.deg - altitude))

    # Some checks & info for labels.
    assert len(time) == len(airmass)

    if hasattr(observer, 'name') is False:
        observer_name = 'No Observer/Location Name'
    else:
        observer_name = observer.name

    if hasattr(target, 'name') is False:
        target_name = str('%.3f' % (random.random()))
    else:
        target_name = target.name

    observe_date = time[0].datetime.strftime('%Y-%m-%d')
    observe_timezone = 'UTC'

    # Plot data.
    ax.plot_date(
        time.plot_date,
        airmass,
        linestyle=style_kwargs['linestyle'],
        color=style_kwargs['color'],
        label=target_name,
        )
    ax.figure.autofmt_xdate()
    ax.invert_yaxis()

    # Set labels, title, legend, etc.
    ax.set_ylabel("Airmass")
    ax.set_xlabel("Hour - "+observe_timezone)
    ax.set_title("Airmass vs Time | "+observer_name+" | "+observe_date+" "+observe_timezone)
    ax.legend(shadow=True)

    return ax


def plot_parallactic(target, observer, time, ax=None, style_kwargs=None, dark_plot=False):
    """
    TODO:
        1) dark_plot style
        2) Parallactic angle equation? Move this to separate function?
        3) observe_timezone -- update with info from observer?
        4) Dealing with NaN or invalid values.

    Returns an ax object with a parallactic angle vs time plot.

    If an ax object is passed in, plots an additional parallactic angle
    plot on top.
    Otherwise, creates a new ax object with a parallactic angle plot.

    Parameters
    ----------
    target: `~astroplan.FixedTarget`
        The celestial body of interest.

    observer: `~astroplan.Observer`
        The person, telescope, observatory, etc. doing the observing.

    time : `~astropy.time.Time`
        An collection of times to use in generating data to plot.
        Can be scalar (or len(time)=1).

    ax : `~matplotlib.axes.Axes` or None
        Optional.

    style_kwargs : WHAT TYPE IS DICTIONARY? or None
        Optional.
        A dictionary with plotting style options.

    dark_plot: Boolean.
        Optional. False (default) or True.
    """

    # Set up axes & plot styles if needed.
    if ax is None:
        ax = plt.gca()
    if style_kwargs is None:
        # Need to add more style options for ax's with existing plots?
        style_kwargs = {'linestyle': '-', 'color': 'r'}
    #if dark_plot is True:
        #Add dark plot stuff here.

    # If Time object is scalar, pad time window.
    if time.isscalar:
        time = time + np.linspace(-12, 12, 100)*u.hour

    # Calculate parallactic angle.
    declination = target.dec
    latitude = observer.location.latitude
    azimuth = observer.altaz(time, target).az
    hour_angle = target.ra

    numerator = np.sin(azimuth) * np.cos(latitude)
    denominator = np.cos(declination) * (np.cos(azimuth)*np.cos(hour_angle) + np.sin(azimuth)*np.sin(hour_angle)*np.sin(latitude))
    parallactic_angle = np.arcsin(numerator/denominator)

    # Some checks & info for labels.
    assert len(time) == len(parallactic_angle)

    if hasattr(observer, 'name') is False:
        observer_name = 'No Observer/Location Name'
    else:
        observer_name = observer.name

    if hasattr(target, 'name') is False:
        target_name = str('%.3f' % (random.random()))
    else:
        target_name = target.name

    observe_date = time[0].datetime.strftime('%Y-%m-%d')
    observe_timezone = 'UTC'

    # Plot data.
    ax.plot_date(
        time.plot_date,
        parallactic_angle,
        linestyle=style_kwargs['linestyle'],
        color=style_kwargs['color'],
        label=target_name,
        )
    ax.figure.autofmt_xdate()

    # Set labels, title, legend, etc.
    ax.set_ylabel("Parallactic Angle")
    ax.set_xlabel("Hour - "+observe_timezone)
    ax.set_title("Parallactic Angle vs Time | "+observer_name+" | "+observe_date+" "+observe_timezone)
    ax.legend(shadow=True)

    return ax


def plot_sky(target, observer, time, ax=None, style_kwargs=None):
    """
    Returns an ax object with a sky plot.

    If an ax object is passed in, plots an additional sky plot on top.
    Otherwise, creates a new ax object with a sky plot.

    Parameters
    ----------
    target: `~astroplan.FixedTarget`
        The celestial body of interest.

    observer: `~astroplan.Observer`
        The person, telescope, observatory, etc. doing the observing.

    time : `~astropy.time.Time`
        An collection of times to use in generating data to plot.
        Can be scalar (len(time)=1).

    ax : `~matplotlib.axes.Axes` or None
        Optional.

    style_kwargs : WHAT TYPE IS DICTIONARY? or None
        Optional.
        A dictionary with plotting style options.
    """

    # Set up axes & plot styles if needed.
    if ax is None:
        # axisbg option sets background color within outer circle (plot bounds).
        # May need to set a color outside of plot bounds.
        #ax = plt.gca(projection='polar', axisbg='0.96')
        ax = plt.gca(polar=True, axisbg='0.96')
    if style_kwargs is None:
        # Need to add more style options for ax's with existing plots?
        style_kwargs = {'marker': 'o', 'color': 'r'}

    # Grab altitude and azimuth from Astroplan objects.
    # Note that values must be made dimensionless before plotting.
    altitude = observer.altaz(time, target).alt*(1/u.deg)
    azimuth = observer.altaz(time, target).az*(1/u.deg)

    # Some checks & info for labels.
    if hasattr(observer, 'name') is False:
        observer_name = 'No Observer/Location Name'
    else:
        observer_name = observer.name

    if hasattr(target, 'name') is False:
        target_name = str('%.3f' % (random.random()))
    else:
        target_name = target.name

    if time.isscalar is True:
        observe_date = time.datetime.strftime('%Y-%m-%d')
    else:
        observe_date = time[0].datetime.strftime('%Y-%m-%d')
    observe_timezone = 'UTC'

    # Set up figure.
    dimensions = (10, 10)
    plt.figure(figsize=dimensions)

    # More axes set-up.
    ax.set_theta_zero_location("N")
    ax.set_rmax(90)

    # Set cardinal directions--DON'T CHANGE.
    ax.annotate('N', (0.5, 1), (0.49, 1.08), textcoords='axes fraction',
                fontsize=16)
    ax.annotate('E', (0, 0.5), (-0.1, 0.49), textcoords='axes fraction',
                fontsize=16)
    ax.annotate('S', (0.5, 0), (0.49, -0.1), textcoords='axes fraction',
                fontsize=16)
    ax.annotate('W', (1, 0.5), (1.09, 0.49), textcoords='axes fraction',
                fontsize=16)

    # Plot star coordinates from dictionary or data object.

    # Coordinates MUST be given to plot() in radians.
    az_rad = (np.pi/180.0)*azimuth

    ax.plot(az_rad, altitude, marker=style_kwargs['marker'],
            color=style_kwargs['color'])

    # # Add time and date labels.
    # if time.isscalar is True:
    #     label = str(target_name)+' - '+str(time)
    #     ax.text(az_rad, altitude, label, withdash=True,
    #             dashdirection=0,
    #             dashlength=20,
    #             rotation=0,
    #             dashrotation=0,
    #             dashpush=10)
    # else:
    #     for t in time:
    #         if t == time[0]:
    #             ax.text(az_rad[t], altitude[t], target_name, withdash=True,
    #                     dashdirection=0,
    #                     dashlength=20,
    #                     rotation=0
    #                     dashrotation=0,
    #                     dashpush=10)
    #         ax.text(az_rad[t], alt[t], t, withdash=True,
    #                 dashdirection=0,
    #                 dashlength=20,
    #                 rotation=0,
    #                 dashrotation=0,
    #                 dashpush=10)

    # Grid, ticks & labels.
    # Set ticks and labels AFTER plotting points.
    ax.grid(True, which='major', color='orange', linewidth=1.5)
    degree_sign = u'\N{DEGREE SIGN}'

    r_labels = [
        '90'+degree_sign, '',
        '60'+degree_sign, '',
        '30'+degree_sign, '',
        '0'+degree_sign+' - Alt.',
        ]
    theta_labels = [
        '0'+degree_sign+' - Az.',
        '45'+degree_sign,
        '90'+degree_sign,
        '135'+degree_sign,
        '180'+degree_sign,
        '225'+degree_sign,
        '270'+degree_sign,
        '',
        ]

    # plt.rgrids(range(1, 105, 15), r_labels)
    # plt.thetagrids(range(0, 360, 45), theta_labels)

    # An alternate way to set ticks and labels.
    ax.set_rgrids(range(1, 105, 15), r_labels, angle=-45)
    ax.set_thetagrids(range(0, 360, 45), theta_labels)

    # Title, legend, etc.
    plt.title('Sky Chart | '+observer_name+' | '+observe_date, y=1.15, size=20)

    return ax
