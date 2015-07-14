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
        3) Dark plot option.
        4) Dealing with no target name.
        5) ax.figure.autofmt_xdate() ?

    Plots airmass as a function of time for a given target.

    If an ax object already exists, an additional airmass plot will be
    "stacked" on it.
    Otherwise, creates a new ax object and plots airmass on top of that.

    When a scalar `Time` object is passed in (e.g., Time('2000-1-1')), the
    resulting plot will use a 24-hour window centered on the time indicated,
    with airmass sampled at regular intervals throughout.
    However, the user can control the exact number and frequency of airmass
    calculations used by passing in a non-scalar `Time` object. For instance,
    Time(['2000-1-1 23:00:00', '2000-1-1 23:30:00']) will result in a plot with
    only two airmass measurements.

    Parameters
    ----------
    target : `~astroplan.FixedTarget`
        The celestial body of interest.

    observer : `~astroplan.Observer`
        The person, telescope, observatory, etc. doing the observing.

    time : `~astropy.time.Time`
        Can be scalar (e.g., Time('2000-1-1')) or not (e.g., Time(['2000-1-1'])).

    ax : `~matplotlib.axes.Axes` or None, optional.
        The axes object to be drawn on.
        If None, use the current axes (matplotlib.pyplot.gca).

    style_kwargs : dict or Empty, optional.
        A dictionary of keywords passed into `matplotlib.pyplot.plot_date`
        to set plotting styles.

    dark_plot : bool, optional.
        Default is False.

    Returns
    -------
    ax :  `~matplotlib.axes.Axes`
        An axes object with added airmass vs. time plot.

    Notes
    -----
    y-axis defaults (if user wishes to change these, use `ax.<set attribute>`
    before drawing or saving plot):
        Inverted by default.
        Lower limit is 3.0.
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
    style_kwargs.setdefault('color', 'k')
    style_kwargs.setdefault('fmt', '-')

    # Populate time window if needed.
    if time.isscalar:
        time = time + np.linspace(-12, 12, 100)*u.hour

    # Calculate airmass
    airmass = observer.altaz(time, target).secz

    # Some checks & info for labels.
    assert len(time) == len(airmass)

    if not hasattr(target, 'name'):
        target_name = ''
    else:
        target_name = target.name
    style_kwargs.setdefault('label', target_name)

    observe_date = time[0].datetime.strftime('%Y-%m-%d')
    observe_timezone = 'UTC'

    # Plot data.
    ax.plot_date(time.plot_date, airmass, **style_kwargs)
    ax.figure.autofmt_xdate()
    if ax.get_ylim()[1] > ax.get_ylim()[0]:
        ax.invert_yaxis()
    if ax.get_ylim()[0] > 3.0:
        upper = ax.get_ylim()[1]
        ax.set_ylim([3, upper])

    # Set labels.
    ax.set_ylabel("Airmass")
    ax.set_xlabel("Time - "+observe_timezone)

    # Output.
    return ax


def plot_parallactic(target, observer, time, ax=None, style_kwargs=None,
                     dark_plot=False):
    """
    TODO:
        1) dark_plot style
        2) Use parallactic angle function.
        3) observe_timezone -- update with info from observer?
        4) Dealing with NaN or invalid values.
        5) Dealing with no target name.
        6) ax.figure.autofmt_xdate() ?

    Plots parallactic angle as a function of time for a given target.

    If an ax object already exists, an additional parallactic angle plot will
    be "stacked" on it.
    Otherwise, creates a new ax object and plots on top of that.

    When a scalar `Time` object is passed in (e.g., Time('2000-1-1')), the
    resulting plot will use a 24-hour window centered on the time indicated,
    with parallactic angle sampled at regular intervals throughout.
    However, the user can control the exact number and frequency of parallactic
    angle calculations used by passing in a non-scalar `Time` object. For
    instance, Time(['2000-1-1 23:00:00', '2000-1-1 23:30:00']) will result in a
    plot with only two parallactic angle measurements.

    Parameters
    ----------
    target: `~astroplan.FixedTarget`
        The celestial body of interest.

    observer: `~astroplan.Observer`
        The person, telescope, observatory, etc. doing the observing.

    time : `~astropy.time.Time`
        Can be scalar (e.g., Time('2000-1-1')) or not (e.g., Time(['2000-1-1'])).

    ax : `~matplotlib.axes.Axes` or None, optional.
        The axes object to be drawn on.
        If None, use the current axes (matplotlib.pyplot.gca).

    style_kwargs : dict or Empty, optional.
        A dictionary of keywords passed into `matplotlib.pyplot.plot_date`
        to set plotting styles.

    dark_plot: bool, optional.
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
    style_kwargs.setdefault('color', 'k')
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
    stuff_to_add = np.sin(azimuth)*np.sin(hour_angle)*np.sin(latitude)
    not_hair_product = np.cos(azimuth)*np.cos(hour_angle)
    denominator = np.cos(declination) * (not_hair_product + stuff_to_add)
    parallactic_angle = np.arcsin(numerator/denominator)

    # Some checks & info for labels.
    assert len(time) == len(parallactic_angle)

    if not hasattr(target, 'name'):
        target_name = ''
    else:
        target_name = target.name
    style_kwargs.setdefault('label', target_name)

    observe_date = time[0].datetime.strftime('%Y-%m-%d')
    observe_timezone = 'UTC'

    # Plot data.
    ax.plot_date(time.plot_date, parallactic_angle, **style_kwargs)
    ax.figure.autofmt_xdate()

    # Set labels.
    ax.set_ylabel("Parallactic Angle")
    ax.set_xlabel("Time - "+observe_timezone)

    return ax
