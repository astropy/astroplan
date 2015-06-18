# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

#from ..extern import six

#import sys
#from math import sqrt, pi, exp, log, floor
from abc import ABCMeta, abstractmethod

#import numpy as np

# Originally authored by Jazmin Berlanga Medina (jazmin.berlanga@gmail.com),
# and modified by Adrian Price-Whelan (email) and Erik Tollerud (email).

__all__ = ["Observer", "Target", "FixedTarget", "NonFixedTarget",
           "Constraint", "TimeWindow", "AltitudeRange",
           "AboveAirmass", "Observation"]

#__doctest_requires__ = {'*': ['scipy.integrate']}


class Observer(object):
    """
    Some comments.
    """

    def __init__(self, name, longitude, latitude, elevation, pressure,
                 relative_humidity, temperature, timezone, description=None):
        """
        Initializes an Observer object.

        <longer description>

        Parameters
        ----------
        name : str
            A short name for the telescope, observatory or location.

        longitude : str or `~astropy.units.quantity`?
            The longitude of the observing location.

        latitude : str or `~astropy.units.quantity`?
            The latitude of the observing location.

        elevation : `~astropy.units.Quantity`
            The elevation of the observing location, with respect to sea
            level.

        pressure : `~astropy.units.Quantity`
            The ambient pressure.

        relative_humidity : float
            The ambient relative humidity.

        temperature : `~astropy.units.Quantity`
            The ambient temperature.

        timezone : WHAT TYPE IS pytz.timezone ?
            The local timezone.

        description : str
            A short description of the telescope, observatory or observing
            location.
        """
        raise NotImplementedError()

    def set_environment(self, pressure, relative_humidity, temperature):
        """
        Updates the Observer object with environmental conditions.

        <longer description>

        Parameters
        ----------
        pressure : `~astropy.units.Quantity`
            The ambient pressure.

        relative_humidity : `~astropy.units.Quantity`
            The ambient relative humidity.

        temperature :
            The ambient temperature.
        """
        raise NotImplementedError()

    def get_date(self, input_date_time, timezone=None):
        """
        Builds an object containing date and time information.

        Default time zone is UTC.
        If time zone (local or otherwise) requested, then uses `datetime.datetime`
        and pytz.timezone to convert from UTC.

        Must first create date_time object with get_date() before using any other 
        Observer methods, all of which require one as an argument.

        Parameters
        ----------
        input_date_time : WHAT TYPE IS astropy.time OBJECT?

        """
        raise NotImplementedError()

    # Sun-related methods.

    def noon(date_time):
        """
        Returns the local, solar noon time.

        Parameters
        ----------
        date_time : WHAT TYPE IS date_time OBJECT?
        """
        raise NotImplementedError()

    def midnight(date_time):
        """
        Returns the local, solar midnight time.

        Parameters
        ----------
        date_time : WHAT TYPE IS date_time OBJECT?
        """
        raise NotImplementedError()

    def sunset(date_time, **kwargs):
        """
        Returns the local sunset time.

        The default sunset returned is the next one to occur.

        Parameters
        ----------
        date_time : WHAT TYPE IS date_time OBJECT?

        Keywords: str, optional
            previous
            next
        """
        raise NotImplementedError()

    def sunrise(date_time, **kwargs):
        """
        Returns the local sunrise time.

        The default sunrise returned is the next one to occur.

        Parameters
        ----------
        date_time : WHAT TYPE IS date_time OBJECT?

        Keywords: str, optional
            previous
            next
        """
        raise NotImplementedError()

    # Moon-related methods.

    def moonrise(date_time, **kwargs):
        """
        Returns the local moonrise time.

        The default moonrise returned is the next one to occur.

        Parameters
        ----------
        date_time : WHAT TYPE IS date_time OBJECT?

        Keywords: str, optional
            previous
            next
        """
        raise NotImplementedError()

    def moonset(date_time, **kwargs):
        """
        Returns the local moonset time.

        The default moonset returned is the next one to occur.

        Parameters
        ----------
        date_time : WHAT TYPE IS date_time OBJECT?

        Keywords: str, optional
            previous
            next
        """
        raise NotImplementedError()

    def moon_illumination(date_time):
        """
        Returns a float giving the percent illumation.

        Parameters
        ----------
        date_time : WHAT TYPE IS date_time OBJECT?
        """
        raise NotImplementedError()

    def moon_position(date_time):
        """
        Returns the position of the moon in alt/az.

        Parameters
        ----------
        date_time : WHAT TYPE IS date_time OBJECT?
        """
        raise NotImplementedError()

    # Other time-related methods.

    def evening_nautical(date_time):
        """
        Returns the local evening (nautical) twilight time.

        Parameters
        ----------
        date_time : WHAT TYPE IS date_time OBJECT?
        """
        raise NotImplementedError()

    def evening_civil(date_time):
        """
        Returns the local evening (civil) twilight time.

        Parameters
        ----------
        date_time : WHAT TYPE IS date_time OBJECT?
        """
        raise NotImplementedError()

    def morning_astronomical(date_time):
        """
        Returns the local morning (astronomical) twilight time.

        Parameters
        ----------
        date_time : WHAT TYPE IS date_time OBJECT?
        """
        raise NotImplementedError()

    def morning_nautical(date_time):
        """
        Returns the local morning (nautical) twilight time.

        Parameters
        ----------
        date_time : WHAT TYPE IS date_time OBJECT?
        """
        raise NotImplementedError()

    def morning_civil(date_time):
        """
        Returns the local morning (civil) twilight time.

        Parameters
        ----------
        date_time : WHAT TYPE IS date_time OBJECT?
        """
        raise NotImplementedError()


class Target(object):
    """
    This is an abstract base class -- you can't instantiate
    examples of this class, but must work with one of its
    subclasses such as `FixedTarget` or `NonFixedTarget`.

    Would need to import six, abc to make this a metaclass?
    """
    __metaclass__ = ABCMeta

    def __init__(self, name, ra, dec, marker=None):
        """
        Defines a single observation target.

        Parameters
        ----------
        name : str, optional

        ra : WHAT TYPE IS ra ?

        dec : WHAT TYPE IS dec ?

        marker : str, optional
            User-defined markers to differentiate between different types
            of targets (e.g., guides, high-priority, etc.).
        """
        raise NotImplementedError()

    @property
    def ra(self):
        """
        Right ascension.
        """
        raise NotImplementedError

    @property
    def dec(self):
        """
        Declination.
        """
        raise NotImplementedError


class FixedTarget(Target):
    """
    An object that is "fixed" with respect to the celestial sphere.
    """


class NonFixedTarget(Target):
    """
    Placeholder for future function.
    """


class Constraint(object):
    """
    An object containing observational constraints.

    A Constraints object is used in conjunction with a Target
    and an Observer object (via the apply_constraints method) to find out
    if a particular target is visible to the observer.
    """
    __metaclass__ = ABCMeta

    @abstractmethod
    def apply_constraints(self, target, observer, constraint_list):
        """
        Returns information on a target's visibility.

        Finds out if a Target is observable by an Observer given a list
        of Constraint objects.  The list must contain at least one
        Constraint object.

        Parameters
        ----------
        target : WHAT TYPE IS Target OBJECT ?

        constraint_list : WHAT TYPE IS constraint_list ? `numpy.array` ??
        """
        raise NotImplementedError


class TimeWindow(Constraint):
    """
    An object containing start and end times for an observation.
    """

    def __init__(self, start, end):
        """
        Initializes a TimeWindow object.

        Parameters
        ----------
        start : STRING OR astropy.time OBJECT ?

        end : STRING OR astropy.time OBJECT ?
        """
        raise NotImplementedError


class AltitudeRange(Constraint):
    """
    An object containing upper and lower altitude limits.
    """

    def __init__(self, low, high):
        """
        Initializes an AltitudeRange object.

        Parameters
        ----------
        low : `~astropy.units.Quantity`

        high : `~astropy.units.Quantity`
        """
        raise NotImplementedError


class AboveAirmass(Constraint):
    """
    An object containing an airmass lower limit.
    """

    def __init__(self, low):
        """
        Initializes an AboveAirmass object.

        Parameters
        ----------
        low : float
        """
        raise NotImplementedError


class Observation(object):
    """
    Comments.
    """

    def __init__(self, target, date_time):
        """
        Initializes an Observation object.

        Parameters
        ----------
        target : WHAT TYPE IS Target OBJECT ?

        date : WHAT TYPE IS date OBJECT ?
        """
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
