# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import copy

from astropy.visualization import astropy_mpl_style
from astropy.utils import minversion

# This returns False if matplotlib cannot be imported
MATPLOTLIB_GE_1_5 = minversion('matplotlib', '1.5')


__all__ = ['light_style_sheet', 'dark_style_sheet',
           'available_style_sheets']

available_style_sheets = ["light_style_sheet", "dark_style_sheet"]
"""Matplotlib style sheets available within astroplan."""

light_style_sheet = copy.deepcopy(astropy_mpl_style)

# One update specifically for astroplan.plots: raise bottom up from default 0.1
light_style_sheet.update({'figure.subplot.bottom': 0.15})

# Define dark style sheet starting from the astropy_mpl_style sheet
dark_style_sheet = {

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

color_cycle = ['#00CCFF',   # blue
               '#FF0000',   # red
               '#40C000',   # green
               '#E066E0',   # purple
               '#CF4457',   # pink
               '#188487',   # turquoise
               '#E24A33']   # orange

if MATPLOTLIB_GE_1_5:
    # This is a dependency of matplotlib, so should be present.
    from cycler import cycler
    dark_style_sheet['axes.prop_cycle'] = cycler('color', color_cycle)

else:
    dark_style_sheet['axes.color_cycle'] = color_cycle
