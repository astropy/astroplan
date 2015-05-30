# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ..extern import six

import sys
from math import sqrt, pi, exp, log, floor
from abc import ABCMeta, abstractmethod

import numpy as np

# Originally authored by Jazmin Berlanga Medina (jazmin.berlanga@gmail.com),
# and modified by Adrian Price-Whelan (email) and Erik Tollerud (email).

__all__ = ["Observer", "Date", "Target", "FixedTarget", "NonFixedTarget",
           "Constraints", "Observation"]

#__doctest_requires__ = {'*': ['scipy.integrate']}


class Observer(object):
    """
    Some comments.
    """
    #def __init__(self, name, longitude, latitude, elevation, pressure,
    #             rel_humidity, temp, timezone, description=None):
    """
    Comments.
    """

    #def set_environment():
    """
    Comments.
    """

    raise NotImplementedError()


class Date(Observer):
    """
    Inherits from Observer.
    """
    #def get_date(self):
    """
    Comments.
    """

    #def morning_astronomical(Date):
    """
    Comments.
    """

    #def next_sunset(Date):
    """
    Comments.
    """

    #def previous_sunset(Date):
    """
    Comments.
    """

    #def next_sunrise(Date):
    """
    Comments.
    """

    #def previous_sunrise(Date):
    """
    Comments.
    """

    #def next_moonrise(Date):
    """
    Comments.
    """

    #def previous_moonrise(Date):
    """
    Comments.
    """

    #def next_moonset(Date):
    """
    Comments.
    """

    #def previous_moonset(Date):
    """
    Comments.
    """

    #def evening_nautical(Date):
    """
    Comments.
    """

    #def evening_civil(Date):
    """
    Comments.
    """

    #def morning_nautical(Date):
    """
    Comments.
    """

    #def morning_civil(Date):
    """
    Comments.
    """

    #def moon_illumination(Date):
    """
    Comments.
    """

    #def moon_position(Date):
    """
    Comments.
    """

    #def noon(Date):
    """
    Comments.
    """

    #def midnight(Date):
    """
    Comments.
    """

    raise NotImplementedError()


@six.add_metaclass(ABCMeta)
class Target(object):
    """ 
    This is an abstract base class -- you can't instantiate
    examples of this class, but must work with one of its
    subclasses such as `FixedTarget` or `NonFixedTarget` 
    """
    #def __init__(self, name, ra, dec, marker=none):

    raise NotImplementedError()


class FixedTarget(Target):
    """
    Comments.
    """
    raise NotImplementedError()


class NonFixedTarget(Target):
    """
    Placeholder for future function.
    """


class Constraints(object):
    """
    object/methods need clarification
    """
    raise NotImplementedError()

class Observation(object):
    """
    Comments.    
    """

    def __init__(self, target, date):

    @property
    def airmass(self):
        """
        Airmass.
        """

    @property
    def pang(self):
        """
        Parallactic angle. 
        """

    @property
    def ha(self):
        """ 
        Hour angle.
        """

    @property
    def ut(self):
        """
        Time of observation in UST.
        """

    @property
    def lt(self):
        """
        Time of observation in local time.
        """

    @property
    def gmst(self):
        """
        Time of observation in GMST.
        """

    @property
    def lmst(self):
        """
        Time of observation in local mean sidereal time.
        """

    @property
    def moon_sep(self):
        """
        Separation between moon and object at time of observation.
        """

    @property
    def alt(self):
        """
        Altitude at time of observation.
        """

    @property
    def az(self):
        """
        Azimuth at time of observation.
        """

