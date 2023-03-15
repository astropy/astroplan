# Licensed under a 3-clause BSD style license - see LICENSE.rst
import copy
from xmlrpc.client import boolean
import numpy as np
import astropy.units as u
from astropy.time import Time
from collections.abc import Sequence
import warnings

from ..exceptions import PlotWarning
from ..utils import _set_mpl_style_sheet
from ..constraints import Constraint

__all__ = ['plot_airmass', 'plot_schedule_airmass', 'plot_parallactic',
           'plot_altitude']


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


def _shade_nights(ax, observer, start, end, tzoffset=0*u.hour):
    """
    Shade the time between sunset and sunrise, including twilight.
    """
    start_jd = start.utc.jd
    end_jd = end.utc.jd
    if end_jd <= start_jd:
        return

    # Sample the Sun's altitude every five minutes.
    n_intervals = max(1, int(np.ceil((end_jd - start_jd) * 24 * 12)))
    sample_jds = np.linspace(start_jd, end_jd, n_intervals + 1)
    times = Time(sample_jds, format='jd', scale='utc')
    sun_altitudes = observer.sun_altaz(times).alt.to_value(u.deg)

    horizons = np.array([-18, -12, -6, 0])
    bands = np.searchsorted(horizons, sun_altitudes)
    changes = np.flatnonzero(np.diff(bands))

    # Interpolate the times when the Sun crosses a twilight boundary.
    crossing_altitudes = horizons[np.minimum(bands[changes], bands[changes + 1])]
    fractions = ((crossing_altitudes - sun_altitudes[changes]) /
                 (sun_altitudes[changes + 1] - sun_altitudes[changes]))
    crossings = (sample_jds[changes] +
                 fractions * (sample_jds[changes + 1] - sample_jds[changes]))

    boundaries = np.concatenate(([start_jd], crossings, [end_jd]))
    alphas = np.array([0.4, 0.3, 0.2, 0.1, 0.0])
    alphas = alphas[np.r_[bands[0], bands[changes + 1]]]
    plot_times = (Time(boundaries, format='jd', scale='utc') + tzoffset).plot_date

    for left, right, alpha in zip(plot_times[:-1], plot_times[1:], alphas):
        if alpha > 0:
            ax.axvspan(left, right, ymin=0, ymax=1, color='grey', alpha=alpha, linewidth=0)


def plot_airmass(targets, observer, time, ax=None, style_kwargs=None,
                 style_sheet=None, brightness_shading=False,
                 altitude_yaxis=False, min_airmass=1.0, min_region=None,
                 max_airmass=3.0, max_region=None, use_local_tz=False,
                 constraints=None):
    r"""
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

    For examples with plots, visit the documentation of
    :ref:`plots_time_dependent`.

    Parameters
    ----------
    targets : list of `~astroplan.FixedTarget` objects
        The celestial bodies of interest.
        If a single object is passed it will be converted to a list.

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
        A dictionary of keywords passed into `~matplotlib.pyplot.plot`
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

    min_airmass : float
        Lower limit of y-axis airmass range in the plot. Default is ``1.0``.

    max_airmass : float
        Upper limit of y-axis airmass range in the plot. Default is ``3.0``.

    min_region : float
        If set, defines an interval between ``min_airmass`` and ``min_region``
        that will be shaded. Default is `None`.

    max_region : float
        If set, defines an interval between ``max_airmass`` and ``max_region``
        that will be shaded. Default is `None`.

    use_local_tz : bool
        If the time is specified in a local timezone, the time will be plotted
        in that timezone.
    
    constraints : Constraint or list of Constraint objects (optional)
        If a list of boolean constraints is provided, the objects will have
        blackened curves outside the constraints.

    Returns
    -------
    ax : `~matplotlib.axes.Axes`
        An ``Axes`` object with added airmass vs. time plot.

    Notes
    -----
    y-axis is inverted and shows airmasses between 1.0 and 3.0 by default.
    If user wishes to change these, use ``ax.<set attribute>`` before drawing
    or saving plot:

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

    if 'lw' not in style_kwargs:
        style_kwargs.setdefault('linewidth', 1.5)
    if 'ls' not in style_kwargs and 'linestyle' not in style_kwargs:
        style_kwargs.setdefault('linestyle', '-')

    if hasattr(time, 'utcoffset') and use_local_tz:
        tzoffset = time.utcoffset()
        tzname = time.tzname()
    else:
        tzoffset = 0 * u.hour
        tzname = 'UTC'
    # Populate time window if needed.
    # (plot against local time if that's requested)
    time_ut = Time(time)
    if time_ut.isscalar:
        time_ut = time_ut + np.linspace(-12, 12, 100)*u.hour
    elif len(time_ut) == 1:
        warnings.warn('You used a Time array of length 1.  You probably meant '
                      'to use a scalar. (Or maybe a list with length > 1?).',
                      PlotWarning)
    timetoplot = time_ut + tzoffset

    if not isinstance(targets, Sequence):
        targets = [targets]
    
    # Firstly, check that the constraints are valid if any.
    if constraints is not None:
        if isinstance(constraints, Constraint):
            constraints = [constraints]
    
        assert isinstance(constraints, (list, tuple, np.ndarray)), "constraints must be an array-like object or a single Constraint instance."
        for cons in constraints:
            assert isinstance(cons, Constraint), "All constraint objects must be Constraints"

            # Make sure they're boolean
            assert hasattr(cons, 'boolean_constraint'), "Each constraints element must have the 'boolean_constraint attribute."
            assert cons.boolean_constraint, "Cannot work with non-boolean constraints."
            # Then make sure they're not some weird user-generated constraint. 
            assert hasattr(cons, 'compute_constraint'), "Constraints must have a 'compute_constraint' method."
            


    for target in targets:
        # Calculate airmass
        airmass = observer.altaz(time_ut, target).secz
        # Mask out nonsense airmasses
        masked_airmass = np.ma.array(airmass, mask=airmass < 1)

        # Some checks & info for labels.
        try:
            target_name = target.name
        except AttributeError:
            target_name = ''

        # Plot data (against timezone-offset time)
        ax.plot(timetoplot.plot_date, masked_airmass, label=target_name, **style_kwargs)

        # Plot black curves for where the target is outside constraints
        if constraints:
            cons_mask = np.ones(len(airmass), dtype=bool)
            for cons in constraints:
                cons_mask &= cons(observer, target, time_ut)

            # Handle discontinuous curves:
            time_vals = np.where((~cons_mask)&(airmass>1),timetoplot.plot_date, np.nan)
            am_vals = np.where((~cons_mask)&(airmass>1), airmass, np.nan)
            ax.plot_date(time_vals, am_vals, "k", lw=3)

        # Plot black curves for where the target is outside constraints
        if constraints:
            cons_mask = np.ones(len(airmass), dtype=bool)
            for cons in constraints:
                cons_mask &= cons(observer, target, time_ut)

            # Handle discontinuous curves:
            time_vals = np.where((~cons_mask)&(airmass>1),timetoplot.plot_date, np.nan)
            am_vals = np.where((~cons_mask)&(airmass>1), airmass, np.nan)
            ax.plot_date(time_vals, am_vals, "k", lw=3)

    # Format the time axis
    xlo, xhi = (timetoplot[0]), (timetoplot[-1])
    ax.xaxis_date()
    ax.set_xlim([xlo.plot_date, xhi.plot_date])
    date_formatter = dates.DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(date_formatter)
    plt.setp(ax.get_xticklabels(), rotation=30, ha='right')

    # Shade background during night time
    if brightness_shading:
        _shade_nights(ax, observer, time_ut.min(), time_ut.max(), tzoffset=tzoffset)

    # Invert y-axis and set limits.
    y_lim = ax.get_ylim()
    if y_lim[1] > y_lim[0]:
        ax.invert_yaxis()
    ax.set_ylim([max_airmass, min_airmass])

    # Draw lo/hi limit regions, if present
    ymax, ymin = ax.get_ylim()       # should be (hi_limit, lo_limit)

    if max_region is not None:
        ax.axhspan(ymax, max_region, facecolor='#F9EB4E', alpha=0.10)
    if min_region is not None:
        ax.axhspan(min_region, ymin, facecolor='#F9EB4E', alpha=0.10)

    # Set labels.
    ax.set_ylabel("Airmass")
    ax.set_xlabel("Time from {0} [{1}]".format(min(timetoplot).datetime.date(), tzname))

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


def plot_altitude(targets, observer, time, ax=None, style_kwargs=None,
                  style_sheet=None, brightness_shading=False,
                  airmass_yaxis=False, min_altitude=0, min_region=None,
                  max_altitude=90, max_region=None):
    r"""
    Plots altitude as a function of time for a given target.

    If a `~matplotlib.axes.Axes` object already exists, an additional
    altitude plot will be "stacked" on it.  Otherwise, creates a new
    `~matplotlib.axes.Axes` object and plots altitude on top of that.

    When a scalar `~astropy.time.Time` object is passed in (e.g.,
    ``Time('2000-1-1')``), the resulting plot will use a 24-hour window
    centered on the time indicated, with altitude sampled at regular
    intervals throughout.
    However, the user can control the exact number and frequency of altitude
    calculations used by passing in a non-scalar `~astropy.time.Time`
    object. For instance, ``Time(['2000-1-1 23:00:00', '2000-1-1
    23:30:00'])`` will result in a plot with only two altitude measurements.

    For examples with plots, visit the documentation of
    :ref:`plots_time_dependent`.

    Parameters
    ----------
    targets : list of `~astroplan.FixedTarget` objects
        The celestial bodies of interest.
        If a single object is passed it will be converted to a list.

    observer : `~astroplan.Observer`
        The person, telescope, observatory, etc. doing the observing.

    time : `~astropy.time.Time`
        If scalar (e.g., ``Time('2000-1-1')``), will result in plotting target
        altitudes once an hour over a 24-hour window.
        If non-scalar (e.g., ``Time(['2000-1-1'])``, ``[Time('2000-1-1')]``,
        ``Time(['2000-1-1', '2000-1-2'])``),
        will result in plotting data at the exact times specified.

    ax : `~matplotlib.axes.Axes` or None, optional.
        The `~matplotlib.axes.Axes` object to be drawn on.
        If None, uses the current ``Axes``.

    style_kwargs : dict or None, optional.
        A dictionary of keywords passed into `~matplotlib.pyplot.plot`
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

    min_altitude : float
        Lower limit of y-axis altitude range in the plot. Default is ``1.0``.

    max_altitude : float
        Upper limit of y-axis altitude range in the plot. Default is ``3.0``.

    min_region : float
        If set, defines an interval between ``min_altitude`` and ``min_region``
        that will be shaded. Default is `None`.

    max_region : float
        If set, defines an interval between ``max_altitude`` and ``max_region``
        that will be shaded. Default is `None`.

    Returns
    -------
    ax : `~matplotlib.axes.Axes`
        An ``Axes`` object with added altitude vs. time plot.
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
    if 'ls' not in style_kwargs and 'linestyle' not in style_kwargs:
        style_kwargs.setdefault('linestyle', '-')
    if 'lw' not in style_kwargs:
        style_kwargs.setdefault('linewidth', 1.5)

    # Populate time window if needed.
    time = Time(time)
    if time.isscalar:
        time = time + np.linspace(-12, 12, 100)*u.hour
    elif len(time) == 1:
        warnings.warn('You used a Time array of length 1.  You probably meant '
                      'to use a scalar. (Or maybe a list with length > 1?).',
                      PlotWarning)

    if not isinstance(targets, Sequence):
        targets = [targets]

    for target in targets:
        # Calculate airmass
        altitude = observer.altaz(time, target).alt
        # Mask out nonsense airmasses
        masked_altitude = np.ma.array(altitude, mask=altitude < 0)

        # Some checks & info for labels.
        try:
            target_name = target.name
        except AttributeError:
            target_name = ''

        # Plot data
        ax.plot(time.plot_date, masked_altitude, label=target_name, **style_kwargs)

    # Format the time axis
    ax.xaxis_date()
    ax.set_xlim([time[0].plot_date, time[-1].plot_date])
    date_formatter = dates.DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(date_formatter)
    plt.setp(ax.get_xticklabels(), rotation=30, ha='right')

    # Shade background during night time
    if brightness_shading:
        _shade_nights(ax, observer, time.min(), time.max())

    # Invert y-axis and set limits.
    # y_lim = ax.get_ylim()
    # if y_lim[1] > y_lim[0]:
    #     ax.invert_yaxis()
    ax.set_ylim([min_altitude, max_altitude])

    # Draw lo/hi limit regions, if present
    ymax, ymin = ax.get_ylim()       # should be (hi_limit, lo_limit)

    if max_region is not None:
        ax.axhspan(ymax, max_region, facecolor='#F9EB4E', alpha=0.10)
    if min_region is not None:
        ax.axhspan(min_region, ymin, facecolor='#F9EB4E', alpha=0.10)

    # Set labels.
    ax.set_ylabel("Altitude")
    ax.set_xlabel("Time from {0} [UTC]".format(min(time).datetime.date()))

    if airmass_yaxis and not _has_twin(ax):
        # altitude_ticks = np.array([90, 60, 50, 40, 30, 20])
        # airmass_ticks = 1./np.cos(np.radians(90 - altitude_ticks))

        airmass_ticks = np.array([1, 2, 3])
        altitude_ticks = 90 - np.degrees(np.arccos(1/airmass_ticks))

        ax2 = ax.twinx()
        # ax2.invert_yaxis()
        ax2.set_yticks(altitude_ticks)
        ax2.set_yticklabels(airmass_ticks)
        ax2.set_ylim(ax.get_ylim())
        ax2.set_ylabel('Airmass')

    # Redraw figure for interactive sessions.
    ax.figure.canvas.draw()

    # Output.
    return ax


def plot_schedule_airmass(schedule, show_night=False):
    """
    Plots when observations of targets are scheduled to occur superimposed
    upon plots of the airmasses of the targets.

    Parameters
    ----------
    schedule : `~astroplan.Schedule`
        a schedule object output by a scheduler
    show_night : bool
        Shades the night-time on the plot

    Returns
    -------
    ax :  `~matplotlib.axes.Axes`
        An ``Axes`` object with added airmass and schedule vs. time plot.
    """
    import matplotlib.pyplot as plt
    from matplotlib import dates
    blocks = copy.copy(schedule.scheduled_blocks)
    sorted_blocks = sorted(schedule.observing_blocks, key=lambda x: x.priority)
    targets = [block.target for block in sorted_blocks]
    ts = (schedule.start_time +
          np.linspace(0, (schedule.end_time - schedule.start_time).value, 100) * u.day)
    targ_to_color = {}
    color_idx = np.linspace(0, 1, len(targets))

    if show_night:
        ax = plt.gca()
        ax.xaxis_date()
        _shade_nights(ax, schedule.observer, schedule.start_time, schedule.end_time)
        if not targets:
            ax.set_xlim(schedule.start_time.plot_date, schedule.end_time.plot_date)
            ax.xaxis.set_major_formatter(dates.DateFormatter('%H:%M'))
            plt.setp(ax.get_xticklabels(), rotation=30, ha='right')
            ax.set_ylim([3.0, 1.0])
            ax.set_ylabel('Airmass')
            ax.set_xlabel('Time from {0} [UTC]'.format(schedule.start_time.datetime.date()))

    # lighter, bluer colors indicate higher priority
    for target, ci in zip(set(targets), color_idx):
        plot_airmass(target, schedule.observer, ts,
                     style_kwargs=dict(color=plt.cm.cool(ci)))
        targ_to_color[target.name] = plt.cm.cool(ci)

    for block in blocks:
        if hasattr(block, 'target'):
            plt.axvspan(block.start_time.plot_date, block.end_time.plot_date,
                        fc=targ_to_color[block.target.name], lw=0, alpha=.6)
        else:
            plt.axvspan(block.start_time.plot_date, block.end_time.plot_date,
                        color='k')
    plt.axhline(3, color='k', label='Transitions')
    # TODO: make this output a `axes` object


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

    For examples with plots, visit the documentation of
    :ref:`plots_time_dependent`.

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
        A dictionary of keywords passed into `~matplotlib.pyplot.plot`
        to set plotting styles.

    style_sheet : dict or `None` (optional)
        matplotlib style sheet to use. To see available style sheets in
        astroplan, print *astroplan.plots.available_style_sheets*. Defaults
        to the light theme.

    Returns
    -------
    ax :  `~matplotlib.axes.Axes`
        An ``Axes`` object with added parallactic angle vs. time plot.

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
    if 'ls' not in style_kwargs and 'linestyle' not in style_kwargs:
        style_kwargs.setdefault('linestyle', '-')

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
    ax.plot(time.plot_date, p_angle, **style_kwargs)

    # Format the time axis
    ax.xaxis_date()
    date_formatter = dates.DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(date_formatter)
    plt.setp(ax.get_xticklabels(), rotation=30, ha='right')

    # Set labels.
    ax.set_ylabel("Parallactic Angle - Radians")
    ax.set_xlabel("Time from {0} [UTC]".format(min(time).datetime.date()))

    # Redraw figure for interactive sessions.
    ax.figure.canvas.draw()

    return ax
