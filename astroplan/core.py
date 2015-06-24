# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from astropy.coordinates import (EarthLocation, Latitude, Longitude, SkyCoord,
                                 AltAz)
import astropy.units as u
from astropy.units import Quantity

import pytz

################################################################################
# TODO: Temporary solution to IERS tables problems
from astropy.utils.data import download_file
from astropy.utils import iers
import datetime
iers.IERS.iers_table = iers.IERS_A.open(download_file(iers.IERS_A_URL, cache=True))
################################################################################

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

    def __init__(self, name=None, location=None, latitude=None, longitude=None,
                 elevation=None, timezone='UTC', pressure=None,
                 relative_humidity=None, temperature=None, description=None):
        """
        Initializes an Observer object.

        TODO: <longer description>

        Parameters
        ----------
        name : str
            A short name for the telescope, observatory or location.

        pressure : `~astropy.units.Quantity`
            The ambient pressure.

        location : `~astropy.coordinates.EarthLocation`
            The location (latitude, longitude, elevation) of the observatory.

        longitude : str or `~astropy.units.Quantity`
            The longitude of the observing location. If str, should be a string
            that initializes a `~astropy.coordinates.Longitude` object with
            units in degrees.

        latitude : str or `~astropy.units.Quantity`
            The latitude of the observing location. If str, should be a string
            that initializes a `~astropy.coordinates.Latitude` object with
            units in degrees.

        elevation : `~astropy.units.Quantity`
            The elevation of the observing location, with respect to sea
            level.

        relative_humidity : float
            The ambient relative humidity.

        temperature : `~astropy.units.Quantity`
            The ambient temperature.

        timezone : str or `datetime.tzinfo`
            The local timezone, as either an instance of `pytz.timezone()` or
            the string accepted by `pytz.timezone()`.

        description : str
            A short description of the telescope, observatory or observing
            location.
        """

        if name is not None:
            self.name = name

        if pressure is None:
            pressure = 0

        if temperature is None:
            temperature = 0 * u.deg_C

        if relative_humidity is None:
            relative_humidity = 0

        self.pressure = pressure
        self.temperature = temperature
        self.relative_humidity = relative_humidity

        # If lat/long given instead of EarthLocation, convert them
        # to EarthLocation
        if location is None and (latitude is not None and longitude is not None):
            accepted_latlon_types = [str, Quantity]
            if (type(latitude) in accepted_latlon_types and
                type(longitude) in accepted_latlon_types):
                latitude = Latitude(latitude, unit=u.degree)
                longitude = Longitude(longitude, unit=u.degree)

            self.location = EarthLocation.from_geodetic(longitude, latitude,
                                                        elevation)

        elif isinstance(location, EarthLocation):
            self.location = location

        else:
            raise TypeError('Observatory location must be specified with '
                            'either (1) latitude and longitude in degrees as '
                            'accepted by astropy.coordinates.Latitude and '
                            'astropy.coordinates.Latitude, or (2) as an '
                            'instance of astropy.coordinates.EarthLocation.')

        # Accept various timezone inputs, default to UTC
        if isinstance(timezone, datetime.tzinfo):
            self.timezone = timezone
        elif isinstance(timezone, str) or isinstance(timezone, unicode):
            self.timezone = pytz.timezone(timezone)
        else:
            raise TypeError('timezone keyword should be a string, or an '
                            'instance of datetime.tzinfo')

    def altaz(self, time, target=None, obswl=None):
        """
        Returns an instance of `~astropy.coordinates.SkyCoord` with altitude and
        azimuth of the `FixedTarget` called `target` at time `time`.

        Parameters
        ----------
        time : `~astropy.time.Time`
            Astropy time object.

        target : None (default) or `~astroplan.FixedTarget` or `~astropy.coordinates.SkyCoord`
            Celestial object of interest. If `target`=None, return just the
            `~astropy.coordinates.AltAz` frame without coordinates.

        """
        if obswl is None:
            obswl = 1*u.micron

        altaz_frame = AltAz(location=self.location, obstime=time,
                            pressure=self.pressure, obswl=obswl,
                            temperature=self.temperature,
                            relative_humidity=self.relative_humidity)

        if target is None:
            return altaz_frame
        else:
            if not (isinstance(target, FixedTarget) or
                        isinstance(target, SkyCoord)):
                raise TypeError('The target must be an instance of FixedTarget '
                                'or SkyCoord.')

            if hasattr(target, 'coord'):
                coordinate = target.coord
            else:
                coordinate = target
            return coordinate.transform_to(altaz_frame)

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

    def __init__(self, name=None, ra=None, dec=None, marker=None):
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
        if isinstance(self, FixedTarget):
            return self.coord.ra
        raise NotImplementedError()

    @property
    def dec(self):
        """
        Declination.
        """
        if isinstance(self, FixedTarget):
            return self.coord.dec
        raise NotImplementedError()


class FixedTarget(Target):
    """
    An object that is "fixed" with respect to the celestial sphere.
    """
    def __init__(self, coord, **kwargs):
        '''
        TODO: Docstring.
        '''
        if not isinstance(coord, SkyCoord):
            raise TypeError('Coordinate must be a SkyCoord.')

        self.name = kwargs.get('name', None)
        self.coord = coord

    @classmethod
    def from_name(cls, query_name, **kwargs):
        # Allow manual override for name keyword so that the target name can
        # be different from the query name, otherwise assume name=queryname.
        name = kwargs.pop('name', query_name)
        return cls(SkyCoord.from_name(query_name), name=name, **kwargs)

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
