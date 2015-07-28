# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u


@u.quantity_input(az_label_offset=u.deg)
def plot_sky(target, observer, time, ax=None, style_kwargs=None,
             north_to_east_ccw=True, grid=True, az_label_offset=0.0*u.deg):
    """
    Plots target positions in the sky with respect to the observer's location.

    If an ax object already exists, plots an additional target position mapped
    on top.
    Otherwise, creates a new ax object with a sky plot.

    Can pass in a scalar `Time` object (e.g., Time('2000-1-1')) or a non-scalar
    one (e.g., Time(['2000-1-1'])).
    If pass in `Time` objects with multiple instances of time
    (e.g., Time(['2000-1-1 20:00:00', '2000-1-1 20:30:00'])), target's position
    will be shown at each of these times.

    Parameters
    ----------
    target: `~astroplan.FixedTarget`
        The celestial body of interest.

    observer: `~astroplan.Observer`
        The person, telescope, observatory, etc. doing the observing.

    time : `~astropy.time.Time`
        Can be scalar (e.g., Time('2000-1-1')) or not
        (e.g., Time(['2000-1-1'])).

    ax : `~matplotlib.axes.Axes` or None, optional.
        The axes object to be drawn on.
        If None, use the current axes (`matplotlib.pyplot.gca`).

    style_kwargs : dict or Empty, optional.
        A dictionary of keywords passed into `matplotlib.pyplot.plot_date`
        to set plotting styles.

    north_to_east_ccw: bool, optional.
        True by default, meaning that azimuth is shown increasing
        counter-clockwise (CCW), or with North at top, East at left, etc.
        To show azimuth increasing clockwise (CW), set to False.

    grid: bool, optional.
        True by default, meaning that grid is drawn.

    az_label_offset: `~astropy.units.degree`, optional.
        DANGER: It is not recommended that you change the default behavior,
        as to do so makes it seem as if N/E/S/W are being decoupled from the
        definition of azimuth (North from az = 0 deg., East from az = 90 deg.,
        etc.).
        An offset for azimuth labels from the North label.  A positive
        offset will increase in the same direction as azimuth
        (see `north_to_east_ccw` option).

    Returns
    -------
    An `Axes` object (ax) with a map of the sky.

    Notes
    -----
    Coordinate defaults:

        Altazimuth (local horizon) coordinate system.  North is always at top
        of plot, South is always at the bottom, E/W can be right or left
        depending on the `north_to_east_cw` option.

        Altitude: 90 degrees (zenith) is at plot origin (center) and 0 degrees
        (horizon) is at plot edge.  This cannot be changed by user.

        Azimuth: 0 degrees is at North (top of plot), 90 degrees at East, etc.
        DANGER: Azimuth labels can be changed by user via the `az_label_offset`
        option, but it is not recommended, as to do so makes it seem as if
        N/E/S/W are being decoupled from the definition of azimuth
        (North from az = 0 deg., East from az = 90 deg., etc.).


    TODO:
        1) Add time/date and target name labels?
    """

    # Set up axes & plot styles if needed.
    if ax is None:
        # axisbg option sets background color within plot bounds (circle).
        # May need to set a color outside of plot bounds.
        ax = plt.gca(projection='polar')
        #ax = plt.gca(polar=True)
    if style_kwargs is None:
        style_kwargs = {}
    style_kwargs = dict(style_kwargs)
    style_kwargs.setdefault('marker', 'o')

    # Grab altitude and azimuth from Astroplan objects.
    # Note that values must be made dimensionless before plotting.
    # Modifying altitude is easier than inverting r-axis.
    altitude = (91 * u.deg - observer.altaz(time, target).alt) *(1/u.deg)
    # Azimuth MUST be given to plot() in radians.
    azimuth = observer.altaz(time, target).az * (1/u.deg) * (np.pi/180.0)

    # Some checks & info for labels.
    if not hasattr(target, 'name'):
        target_name = '?'
    else:
        target_name = target.name
    style_kwargs.setdefault('label', target_name)

    # We only want to plot positions above the horizon.
    az_plot = None
    for alt in range(0, len(altitude)):
        if altitude[alt] > 91.0:
            print("Warning: Target " + str(target_name) +
                  " is below the horizon at time: " + str(time[alt]))
        else:
            if az_plot is None:
                az_plot = np.array([azimuth[alt]])
            else:
                az_plot = np.append(az_plot, azimuth[alt])
    alt_plot = altitude[altitude <= 91.0]
    if az_plot is None:
        az_plot = []

    # More axes set-up.
    # Position of azimuth = 0 (data, not label).
    ax.set_theta_zero_location('N')

    # Direction of azimuth increase. Clockwise is -1
    if north_to_east_ccw is False:
        ax.set_theta_direction(-1)

    # Plot target coordinates.
    # ax.plot connects dots.
    #ax.plot(az_plot, alt_plot, **style_kwargs)
    # ax.scatter does not connect dots.
    ax.scatter(az_plot, alt_plot, **style_kwargs)
    #ax.polar(az_plot, alt_plot, **style_kwargs)

    # Set radial limits.
    ax.set_rlim(1, 91)

    # Grid, ticks & labels.
    # May need to set ticks and labels AFTER plotting points.
    if grid is True:
        ax.grid(True, which='major')
    if grid is False:
        ax.grid(False)
    degree_sign = u'\N{DEGREE SIGN}'

    # For positively-increasing range (e.g., range(1, 90, 15)),
    # labels go from middle to outside.
    r_labels = [
        '90' + degree_sign,
        '',
        '60' + degree_sign,
        '',
        '30' + degree_sign,
        '',
        '0' + degree_sign + ' Alt.',
        ]

    theta_labels = []
    for chunk in range(0, 7):
        label_angle = (az_label_offset*(1/u.deg)) + (chunk*45.0)
        while label_angle >= 360.0:
            label_angle = label_angle - 360.0
        if chunk == 0:
            theta_labels.append('N ' + '\n' + str(label_angle) + degree_sign
                                + ' Az')
        elif chunk == 2:
            theta_labels.append('E' + '\n' + str(label_angle) + degree_sign)
        elif chunk == 4:
            theta_labels.append('S' + '\n' + str(label_angle) + degree_sign)
        elif chunk == 6:
            theta_labels.append('W' + '\n' + str(label_angle) + degree_sign)
        else:
            theta_labels.append(str(label_angle) + degree_sign)

    # Set ticks and labels.
    ax.set_rgrids(range(1, 106, 15), r_labels, angle=-45)
    ax.set_thetagrids(range(0, 360, 45), theta_labels, frac=1.2)

    # Below commands don't seem to work.
    #ax.rgrids(range(1, 91, 15), r_labels, angle=-45)
    #ax.thetagrids(range(0, 360, 45), theta_labels)

    # Redraw the figure for interactive sessions.
    ax.figure.canvas.draw()

    return ax
