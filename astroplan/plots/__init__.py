"""
`astroplan.plots` is a module for making plots of important quantities within
the `astroplan` toolbox, such as airmass vs. time, etc.

Note that, at minimum, you must construct `Observer`, `Target` and
`astropy.time.Time` objects as input for all plots.
"""

from .time_dependent import *
from .sky import *
