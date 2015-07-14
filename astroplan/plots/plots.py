# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
import random


def plot_airmass(target, observer, time, ax=None, style_kwargs=None,
                 dark_plot=False):
    """
    TODO:
        1) Timezones?
        2) Limit airmass <= 3.
        3) Plot lines instead of points.
        4) Legend should not be transparent.
        5) dark_plot option

    Makes an airmass vs. time plot.

    If an ax object is passed in, plots an additional airmass plot on top.
    Otherwise, creates a new ax object with an airmass plot.

    When a `Time` object with a single instance in time is used (e.g., if
    len(time) = 1), a new `Time` object is created internally.  This new object
    contains a sampling of times within a window centered on the input `Time`.

    Parameters
    ----------
    target : `~astroplan.FixedTarget`
        The celestial body of interest.

    observer : `~astroplan.Observer`
        The person, telescope, observatory, etc. doing the observing.

    time : `~astropy.time.Time`
        Used in generating data to plot.
        Can be scalar (len(time)=1) or contain multiple times.

    ax : `~matplotlib.axes.Axes` or None, optional.
        The axes this plot will be added to. If none are passed in, new axes
        will be created.

    style_kwargs : Dictionary or Empty, optional.
        A dictionary of keywords passed into `matplotlib.pyplot.plot_date`
        to set plotting styles.

    dark_plot : Boolean, optional.
        Default is False.

    Returns
    -------
    ax :  `~matplotlib.axes.Axes`
        An axes object with added airmass vs. time plot.
    """

    # Set up plot axes and style if needed.
    if ax is None:
        ax = plt.gca()
    if dark_plot is True:
        raise NotImplementedError
    if style_kwargs is None:
        style_kwargs = {}
    style_kwargs = dict(style_kwargs)
    style_kwargs.setdefault('linestyle', '-')
    style_kwargs.setdefault('color', 'b')
    style_kwargs.setdefault('fmt', '-')

    # Populate time window if needed.
    if time.isscalar:
        time = time + np.linspace(-12, 12, 100)*u.hour

    # Calculate airmass
    airmass = observer.altaz(time, target).secz

    # Some checks & info for labels.
    assert len(time) == len(airmass)

    if hasattr(observer, 'name') is False:
        observer_name = ''
    else:
        observer_name = observer.name

    if hasattr(target, 'name') is False:
        target_name = str('%.3f' % (random.random()))
    else:
        target_name = target.name
    style_kwargs.setdefault('label', target_name)

    observe_date = time[0].datetime.strftime('%Y-%m-%d')
    observe_timezone = 'UTC'

    # Plot data.
    ax.plot_date(time.plot_date, airmass, **style_kwargs)
    ax.figure.autofmt_xdate()
    if plt.ylim()[1] > plt.ylim()[0]:
        ax.invert_yaxis()

    # Set labels, title, legend, etc.
    ax.set_ylabel("Airmass")
    ax.set_xlabel("Time - "+observe_timezone)
    ax.set_title("Airmass vs Time | " + observer_name + " | " + observe_date +
                 " " + observe_timezone)

    # Output.
    return ax


def plot_parallactic(target, observer, time, ax=None, style_kwargs=None,
                     dark_plot=False):
    """
    TODO:
        1) dark_plot style
        2) Parallactic angle equation? Move this to separate function?
        3) observe_timezone -- update with info from observer?
        4) Dealing with NaN or invalid values.

    Makes a parallactic angle vs time plot.

    If an ax object is passed in, plots an additional parallactic angle
    plot on top.
    Otherwise, creates a new ax object with a parallactic angle plot.

    When a `Time` object with a single instance in time is used (e.g., if
    len(time) = 1), a new `Time` object is created internally.  This new object
    contains a sampling of times within a window centered on the input `Time`.

    Parameters
    ----------
    target: `~astroplan.FixedTarget`
        The celestial body of interest.

    observer: `~astroplan.Observer`
        The person, telescope, observatory, etc. doing the observing.

    time : `~astropy.time.Time`
        Used in generating data to plot.
        Can be scalar (or len(time)=1) or contain multiple times.

    ax : `~matplotlib.axes.Axes` or None, optional.
        The axes this plot will be added to. If none are passed in, new axes
        will be created.

    style_kwargs : Dictionary or Empty, optional.
        A dictionary of keywords passed into `matplotlib.pyplot.plot_date`
        to set plotting styles.

    dark_plot: Boolean.
        Optional. False (default) or True.

    Returns
    -------
    ax :  `~matplotlib.axes.Axes`
        An axes object with added parallactic_angle vs. time plot.
    """

    # Set up plot axes and style if needed.
    if ax is None:
        ax = plt.gca()
    if dark_plot is True:
        raise NotImplementedError
    if style_kwargs is None:
        style_kwargs = {}
    style_kwargs = dict(style_kwargs)
    style_kwargs.setdefault('linestyle', '-')
    style_kwargs.setdefault('color', 'b')
    style_kwargs.setdefault('fmt', '-')

    # If Time object is scalar, pad time window.
    if time.isscalar:
        time = time + np.linspace(-12, 12, 100)*u.hour

    # Calculate parallactic angle.
    declination = target.dec.to(u.rad)
    latitude = observer.location.latitude.to(u.rad)
    azimuth = observer.altaz(time, target).az.to(u.rad)

    # Replace with Astroplan parallactic angle function call.
    longitude = observer.location.longitude.to(u.rad).value
    local_sidereal_time = time.sidereal_time('mean', longitude).to(u.rad).value
    hour_angle = local_sidereal_time - target.ra.to(u.rad).value

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
    style_kwargs.setdefault('label', target_name)

    observe_date = time[0].datetime.strftime('%Y-%m-%d')
    observe_timezone = 'UTC'

    # Plot data.
    ax.plot_date(time.plot_date, parallactic_angle, **style_kwargs)
    ax.figure.autofmt_xdate()

    # Set labels, title, legend, etc.
    ax.set_ylabel("Parallactic Angle")
    ax.set_xlabel("Time - "+observe_timezone)
    ax.set_title("Parallactic Angle vs Time | " + observer_name + " | " +
                 observe_date + " " + observe_timezone)

    return ax
