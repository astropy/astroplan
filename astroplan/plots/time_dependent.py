# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import get_sun
import warnings

from ..exceptions import PlotWarning
from ..utils import _set_mpl_style_sheet

__all__ = ['plot_airmass', 'plot_parallactic']


def _secz_to_altitude(secant_z):
    """
    Convert airmass (approximated as the secant of the zenith angle) to
    an altitude (aka elevation) in degrees.

    Parameters
    ----------
    secant_z : float
        Secant of the zenith angle

    Returns
    -------
    altitude : float
        Altitude [degrees]
    """
    return np.degrees(np.pi/2 - np.arccos(1./secant_z))


def _has_twin(ax):
    """
    Solution for detecting twin axes built on `ax`. Courtesy of
    Jake Vanderplas http://stackoverflow.com/a/36209590/1340208
    """
    for other_ax in ax.figure.axes:
        if other_ax is ax:
            continue
        if other_ax.bbox.bounds == ax.bbox.bounds:
            return True
    return False


def plot_airmass(target, observer, time, ax=None, style_kwargs=None,
                 style_sheet=None, brightness_shading=False,
                 altitude_yaxis=False):
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

    For examples with plots, visit the astroplan Read the Docs
    documentation [1]_.

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

    style_sheet : dict or `None` (optional)
        matplotlib style sheet to use. To see available style sheets in
        astroplan, print *astroplan.plots.available_style_sheets*. Defaults
        to the light theme.

    brightness_shading : bool
        Shade background of plot to scale roughly with sky brightness. Dark
        shading signifies times when the sun is below the horizon. Default
        is `False`.

    altitude_yaxis : bool
        Add alternative y-axis on the right side of the figure with target
        altitude. Default is `False`.

    Returns
    -------
    ax : `~matplotlib.axes.Axes`
        An ``Axes`` object with added airmass vs. time plot.

    Notes
    -----
    y-axis is inverted and shows airmasses between 1.0 and 3.0 by default.
    If user wishes to change these, use ``ax.\<set attribute\>`` before drawing
    or saving plot:

    References
    ----------
    .. [1] astroplan plotting tutorial: https://astroplan.readthedocs.io/en/latest/tutorials/plots.html#time-dependent-plots

    """
    # Import matplotlib, set style sheet
    if style_sheet is not None:
        _set_mpl_style_sheet(style_sheet)

    import matplotlib.pyplot as plt
    from matplotlib import dates

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
    # Mask out nonsense airmasses
    masked_airmass = np.ma.array(airmass, mask=airmass < 1)

    # Some checks & info for labels.
    if not hasattr(target, 'name'):
        target_name = ''
    else:
        target_name = target.name
    style_kwargs.setdefault('label', target_name)

    # Plot data.
    ax.plot_date(time.plot_date, masked_airmass, **style_kwargs)

    # Format the time axis
    date_formatter = dates.DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(date_formatter)
    plt.setp(ax.get_xticklabels(), rotation=30, ha='right')

    # Shade background during night time
    if brightness_shading:
        for test_time in [time[0], time[-1]]:
            midnight = observer.midnight(test_time)
            previous_sunset = observer.sun_set_time(midnight, which='previous')
            next_sunrise = observer.sun_rise_time(midnight, which='next')

            previous_twilight = observer.twilight_evening_astronomical(midnight, which='previous')
            next_twilight = observer.twilight_morning_astronomical(midnight, which='next')

            ax.axvspan(previous_sunset.plot_date, next_sunrise.plot_date,
                       facecolor='k', alpha=0.05)
            ax.axvspan(previous_twilight.plot_date, next_twilight.plot_date,
                       facecolor='k', alpha=0.05)

    # Invert y-axis and set limits.
    if ax.get_ylim()[1] > ax.get_ylim()[0]:
        ax.invert_yaxis()
    ax.set_ylim([3, 1])
    ax.set_xlim([time[0].plot_date, time[-1].plot_date])

    # Set labels.
    ax.set_ylabel("Airmass")
    ax.set_xlabel("Time from {0} [UTC]".format(min(time).datetime.date()))

    if altitude_yaxis and not _has_twin(ax):
        altitude_ticks = np.array([90, 60, 50, 40, 30, 20])
        airmass_ticks = 1./np.cos(np.radians(90 - altitude_ticks))

        ax2 = ax.twinx()
        ax2.invert_yaxis()
        ax2.set_yticks(airmass_ticks)
        ax2.set_yticklabels(altitude_ticks)
        ax2.set_ylim(ax.get_ylim())
        ax2.set_ylabel('Altitude [degrees]')

    # Redraw figure for interactive sessions.
    ax.figure.canvas.draw()

    # Output.
    return ax


def plot_parallactic(target, observer, time, ax=None, style_kwargs=None,
                     style_sheet=None):
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

    For examples with plots, visit the astroplan Read the Docs
    documentation [1]_.

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

    style_sheet : dict or `None` (optional)
        matplotlib style sheet to use. To see available style sheets in
        astroplan, print *astroplan.plots.available_style_sheets*. Defaults
        to the light theme.

    Returns
    -------
    ax :  `~matplotlib.axes.Axes`
        An ``Axes`` object with added parallactic angle vs. time plot.

    References
    ----------
    .. [1] astroplan plotting tutorial: https://astroplan.readthedocs.io/en/latest/tutorials/plots.html#time-dependent-plots
    """
    # Import matplotlib, set style sheet
    if style_sheet is not None:
        _set_mpl_style_sheet(style_sheet)

    import matplotlib.pyplot as plt

    from matplotlib import dates

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

    # Plot data.
    ax.plot_date(time.plot_date, p_angle, **style_kwargs)

    # Format the time axis
    date_formatter = dates.DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(date_formatter)
    plt.setp(ax.get_xticklabels(), rotation=30, ha='right')

    # Set labels.
    ax.set_ylabel("Parallactic Angle - Radians")
    ax.set_xlabel("Time from {0} [UTC]".format(min(time).datetime.date()))

    # Redraw figure for interactive sessions.
    ax.figure.canvas.draw()

    return ax
