# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from astropy.visualization import astropy_mpl_style

__all__ = ['astroplan_light_style', 'astroplan_dark_style',
           'available_style_sheets']

available_style_sheets = ["astroplan_light_style", "astroplan_dark_style"]

# One update specifically for astroplan.plots: raise bottom up from default 0.1
astropy_mpl_style.update({'figure.subplot.bottom': 0.15})
astroplan_light_style = astropy_mpl_style

astroplan_dark_style = {

    # Lines
    'lines.linewidth': 1.7,
    'lines.antialiased': True,

    # Patches
    'patch.linewidth': 1.0,
    'patch.facecolor': 'k',
    'patch.edgecolor': 'k',
    'patch.antialiased': True,

    # images
    'image.cmap': 'gist_heat',
    'image.origin': 'upper',

    # Font
    'font.size': 12.0,

    # Axes
    'axes.facecolor': '#000000',
    'axes.edgecolor': '#F2F2F2',
    'axes.linewidth': 1.0,
    'axes.grid': True,
    'axes.titlesize': 'x-large',
    'axes.labelsize': 'large',
    'axes.labelcolor': 'w',
    'axes.axisbelow': True,
    'axes.color_cycle': ['#00FFFF',   # blue
                         '#E066E0',   # purple
                         '#FF8080',   # red
                         '#467821',   # green
                         '#CF4457',   # pink
                         '#188487',   # turquoise
                         '#E24A33'],  # orange

    # Ticks
    'xtick.major.size': 0,
    'xtick.minor.size': 0,
    'xtick.major.pad': 6,
    'xtick.minor.pad': 6,
    'xtick.color': '#F2F2F2',
    'xtick.direction': 'in',
    'ytick.major.size': 0,
    'ytick.minor.size': 0,
    'ytick.major.pad': 6,
    'ytick.minor.pad': 6,
    'ytick.color': '#F2F2F2',
    'ytick.direction': 'in',

    # Legend
    'legend.fancybox': True,
    'legend.loc': 'best',

    # Figure
    'figure.figsize': [8, 6],
    'figure.facecolor': 'k',
    'figure.edgecolor': 'k',
    'figure.subplot.hspace': 0.5,
    'figure.subplot.bottom': 0.15,

    # Saved figures
    'savefig.dpi': 72,
    'savefig.facecolor': 'k',
    'savefig.edgecolor': 'k',
    'savefig.transparent': False,

    # Other
    'grid.color': 'w',
    'text.color': 'w'
}