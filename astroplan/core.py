# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

#from ..extern import six

import sys
#from math import sqrt, pi, exp, log, floor
#from abc import ABCMeta, abstractmethod

import numpy as np

# Originally authored by Jazmin Berlanga Medina (jazmin.berlanga@gmail.com),
# and modified by Adrian Price-Whelan (email) and Erik Tollerud (email).

__all__ = ["Observer", "Target", "FixedTarget", "NonFixedTarget",
           "Constraints", "Observation"]

#__doctest_requires__ = {'*': ['scipy.integrate']}


class Observer(object):
    """
    Some comments.
    """

    def __init__(self, name, longitude, latitude, elevation, pressure,
                 rel_humidity, temp, timezone, description=None):
        """
        Comments.
        """
        raise NotImplementedError()

    def set_environment():
        """
        Comments.
        """
        raise NotImplementedError()

    # Methods to 

    def get_date(self):
        """
        Compatible with astropy.time.Time.

        Must first create date object with get_date() to then feed to all other methods below.
        """
        raise NotImplementedError()

    # Sun-related methods.

    def noon(Date):
        """
        Comments.
        """
        raise NotImplementedError()

    def midnight(Date):
        """
        Comments.
        """
        raise NotImplementedError()

    def sunset(Date, **kwargs):
        """
        Comments.

        Keywords: 
            previous
            next

        Default: 
            ?
        """
        raise NotImplementedError()

    def sunrise(Date, **kwargs):
        """
        Comments.
        
        Keywords: 
            previous
            next

        Default: 
            ?
        """
        raise NotImplementedError()

    # Moon-related methods.

    def moonrise(Date, **kwargs):
        """
        Comments.
        
        Keywords: 
            previous
            next

        Default: 
            ?
        """
        raise NotImplementedError()

    def moonset(Date, **kwargs):
        """
        Comments.

        Keywords: 
            previous
            next

        Default: 
            ?
        """
        raise NotImplementedError()

    def moon_illumination(Date):
        """
        Comments.
        """
        raise NotImplementedError()

    def moon_position(Date):
        """
        Comments.
        """
        raise NotImplementedError()

    # Other time-related methods.

    def evening_nautical(Date):
        """
        Comments.
        """
        raise NotImplementedError()

    def evening_civil(Date):
        """
        Comments.
        """
        raise NotImplementedError()

    def morning_astronomical(Date):
        """
        Comments.
        """
        raise NotImplementedError()

    def morning_nautical(Date):
        """
        Comments.
        """
        raise NotImplementedError()

    def morning_civil(Date):
        """
        Comments.
        """
        raise NotImplementedError()


#@six.add_metaclass(ABCMeta)
class Target(object):
    """ 
    This is an abstract base class -- you can't instantiate
    examples of this class, but must work with one of its
    subclasses such as `FixedTarget` or `NonFixedTarget` 

    Would need to import six, abc to make this a metaclass?
    """
    def __init__(self, name, ra, dec, marker=None):
        raise NotImplementedError()


class FixedTarget(Target):
    """
    Comments.
    """


class NonFixedTarget(Target):
    """
    Placeholder for future function.
    """


class Constraints(object):
    """
    object/methods need clarification
    """
    

class Observation(object):
    """
    Comments.    
    """

    def __init__(self, target, date):
        raise NotImplementedError()

    # Observability properties.

    @property
    def alt(self):
        """
        Altitude at time of observation.
        """
        raise NotImplementedError()

    @property
    def az(self):
        """
        Azimuth at time of observation.
        """
        raise NotImplementedError()

    @property
    def airmass(self):
        """
        Airmass.
        """
        raise NotImplementedError()

    @property
    def pang(self):
        """
        Parallactic angle. 
        """
        raise NotImplementedError()

    @property
    def ha(self):
        """ 
        Hour angle.
        """
        raise NotImplementedError()

    @property
    def moon_sep(self):
        """
        Separation between moon and object at time of observation.
        """
        raise NotImplementedError()

    # Time properties.

    @property
    def ut(self):
        """
        Time of observation in UST.
        """
        raise NotImplementedError()

    @property
    def lt(self):
        """
        Time of observation in local time.
        """
        raise NotImplementedError()

    @property
    def gmst(self):
        """
        Time of observation in GMST.
        """
        raise NotImplementedError()

    @property
    def lmst(self):
        """
        Time of observation in local mean sidereal time.
        """
        raise NotImplementedError()

