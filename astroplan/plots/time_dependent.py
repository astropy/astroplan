# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
import astropy.units as u
from astropy.time import Time
import warnings

from ..exceptions import PlotWarning

__all__ = ['plot_airmass', 'plot_parallactic']


def plot_airmass(target, observer, time, ax=None, style_kwargs=None):
    """
    Plots airmass as a function of time for a given target.

    If a `~matplotlib.axes.Axes` object already exists, an additional
    airmass plot will be "stacked" on it.  Otherwise, creates a new
    `~matplotlib.axes.Axes` object and plots airmass on top of that.

    When a scalar `~astropy.time.Time` object is passed in (e.g.,
    ``Time('2000-1-1')``), the resulting plot will use a 24-hour window
    centered on the time indicated, with airmass sampled at regular
    intervals throughout.
    However, the user can control the exact number and frequency of airmass
    calculations used by passing in a non-scalar `~astropy.time.Time`
    object. For instance, ``Time(['2000-1-1 23:00:00', '2000-1-1
    23:30:00'])`` will result in a plot with only two airmass measurements.

    Parameters
    ----------
    target : `~astroplan.FixedTarget`
        The celestial body of interest.

    observer : `~astroplan.Observer`
        The person, telescope, observatory, etc. doing the observing.

    time : `~astropy.time.Time`
        If scalar (e.g., ``Time('2000-1-1')``), will result in plotting target
        airmasses once an hour over a 24-hour window.
        If non-scalar (e.g., ``Time(['2000-1-1'])``, ``[Time('2000-1-1')]``,
        ``Time(['2000-1-1', '2000-1-2'])``),
        will result in plotting data at the exact times specified.

    ax : `~matplotlib.axes.Axes` or None, optional.
        The `~matplotlib.axes.Axes` object to be drawn on.
        If None, uses the current ``Axes``.

    style_kwargs : dict or None, optional.
        A dictionary of keywords passed into `~matplotlib.pyplot.plot_date`
        to set plotting styles.

    Returns
    -------
    ax : `~matplotlib.axes.Axes`
        An ``Axes`` object with added airmass vs. time plot.

    Notes
    -----
    y-axis is inverted and shows airmasses between 1.0 and 3.0 by default.
    If user wishes to change these, use ``ax.\<set attribute\>`` before drawing
    or saving plot:

    TODO:
        1) Timezones?
        2) Dark plot option.
    """
    import matplotlib.pyplot as plt

    # Set up plot axes and style if needed.
    if ax is None:
        ax = plt.gca()
    if style_kwargs is None:
        style_kwargs = {}
    style_kwargs = dict(style_kwargs)
    style_kwargs.setdefault('linestyle', '-')
    style_kwargs.setdefault('fmt', '-')

    # Populate time window if needed.
    time = Time(time)
    if time.isscalar:
        time = time + np.linspace(-12, 12, 100)*u.hour
    elif len(time) == 1:
        warnings.warn('You used a Time array of length 1.  You probably meant '
                      'to use a scalar. (Or maybe a list with length > 1?).',
                      PlotWarning)

    # Calculate airmass
    airmass = observer.altaz(time, target).secz

    # Some checks & info for labels.
    if not hasattr(target, 'name'):
        target_name = ''
    else:
        target_name = target.name
    style_kwargs.setdefault('label', target_name)

    observe_timezone = ''

    # Plot data.
    ax.plot_date(time.plot_date, airmass, **style_kwargs)
    ax.figure.autofmt_xdate()

    # Invert y-axis and set limits.
    if ax.get_ylim()[1] > ax.get_ylim()[0]:
        ax.invert_yaxis()
    ax.set_ylim([3, 1])

    # Set labels.
    ax.set_ylabel("Airmass")
    ax.set_xlabel("Time - "+observe_timezone)

    # Redraw figure for interactive sessions.
    ax.figure.canvas.draw()

    # Output.
    return ax


def plot_parallactic(target, observer, time, ax=None, style_kwargs=None):
    """
    Plots parallactic angle as a function of time for a given target.

    If a `~matplotlib.axes.Axes` object already exists, an additional
    parallactic angle plot will be "stacked" on it.  Otherwise, creates a
    new `~matplotlib.axes.Axes` object and plots on top of that.

    When a scalar `~astropy.time.Time` object is passed in (e.g.,
    ``Time('2000-1-1')``), the resulting plot will use a 24-hour window
    centered on the time indicated, with parallactic angle sampled at
    regular intervals throughout.
    However, the user can control the exact number and frequency of parallactic
    angle calculations used by passing in a non-scalar `~astropy.time.Time`
    object. For instance, ``Time(['2000-1-1 23:00:00', '2000-1-1 23:30:00'])``
    will result in a plot with only two parallactic angle measurements.

    Parameters
    ----------
    target : `~astroplan.FixedTarget`
        The celestial body of interest.

    observer : `~astroplan.Observer`
        The person, telescope, observatory, etc. doing the observing.

    time : `~astropy.time.Time`
        If scalar (e.g., ``Time('2000-1-1')``), will result in plotting target
        parallactic angle once an hour over a 24-hour window.
        If non-scalar (e.g., ``Time(['2000-1-1'])``, ``[Time('2000-1-1')]``,
        ``Time(['2000-1-1', '2000-1-2'])``),
        will result in plotting data at the exact times specified.

    ax : `~matplotlib.axes.Axes` or None, optional.
        The ``Axes`` object to be drawn on.
        If None, uses the current ``Axes``.

    style_kwargs : dict or None, optional.
        A dictionary of keywords passed into `~matplotlib.pyplot.plot_date`
        to set plotting styles.

    Returns
    -------
    ax :  `~matplotlib.axes.Axes`
        An ``Axes`` object with added parallactic angle vs. time plot.

    TODO:
        1) dark_plot style
        2) observe_timezone -- update with info from observer?
    """
    import matplotlib.pyplot as plt

    # Set up plot axes and style if needed.
    if ax is None:
        ax = plt.gca()
    if style_kwargs is None:
        style_kwargs = {}
    style_kwargs = dict(style_kwargs)
    style_kwargs.setdefault('linestyle', '-')
    style_kwargs.setdefault('fmt', '-')

    # Populate time window if needed.
    time = Time(time)
    if time.isscalar:
        time = time + np.linspace(-12, 12, 100)*u.hour
    elif len(time) == 1:
        warnings.warn('You used a Time array of length 1.  You probably meant '
                      'to use a scalar. (Or maybe a list with length > 1?).',
                      PlotWarning)

    # Calculate parallactic angle.
    p_angle = observer.parallactic_angle(time, target)

    # Some checks & info for labels.
    assert len(time) == len(p_angle)

    if not hasattr(target, 'name'):
        target_name = ''
    else:
        target_name = target.name
    style_kwargs.setdefault('label', target_name)

    observe_timezone = ''

    # Plot data.
    ax.plot_date(time.plot_date, p_angle, **style_kwargs)
    ax.figure.autofmt_xdate()

    # Set labels.
    ax.set_ylabel("Parallactic Angle - Radians")
    ax.set_xlabel("Time - "+observe_timezone)

    # Redraw figure for interactive sessions.
    ax.figure.canvas.draw()

    return ax
