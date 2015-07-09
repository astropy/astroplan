# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from astropy.coordinates import (EarthLocation, Latitude, Longitude, SkyCoord,
                                 AltAz, Angle)
import astropy.units as u
from astropy.units import Quantity

import pytz
import numpy as np
from astropy.time import Time

################################################################################
# TODO: Temporary solution to IERS tables problems
from astropy.utils.data import download_file
from astropy.utils import iers
import datetime
iers.IERS.iers_table = iers.IERS_A.open(download_file(iers.IERS_A_URL,
                                                      cache=True))
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
    A container class for information about an observer's location and
    environment.
    """

    def __init__(self, name=None, location=None, latitude=None, longitude=None,
                 elevation=0*u.m, timezone='UTC', pressure=None,
                 relative_humidity=None, temperature=None, description=None):
        """
        Parameters
        ----------
        name : str
            A short name for the telescope, observatory or location.

        location : `~astropy.coordinates.EarthLocation`
            The location (latitude, longitude, elevation) of the observatory.

        longitude : float, str, `~astropy.units.Quantity` (optional)
            The longitude of the observing location. Should be valid input for
            initializing a `~astropy.coordinates.Longitude` object.

        latitude : float, str, `~astropy.units.Quantity` (optional)
            The latitude of the observing location. Should be valid input for
            initializing a `~astropy.coordinates.Latitude` object.

        elevation : `~astropy.units.Quantity` (optional), default = 0 meters
            The elevation of the observing location, with respect to sea
            level. Defaults to sea level.

        pressure : `~astropy.units.Quantity` (optional)
            The ambient pressure. Defaults to zero (i.e. no atmosphere).

        relative_humidity : float (optional)
            The ambient relative humidity.

        temperature : `~astropy.units.Quantity` (optional)
            The ambient temperature.

        timezone : str or `datetime.tzinfo` (optional)
            The local timezone to assume. If a string, it will be passed through
            `pytz.timezone()` to produce the timezone object.

        description : str (optional)
            A short description of the telescope, observatory or observing
            location.
        """

        self.name = name
        self.pressure = pressure
        self.temperature = temperature
        self.relative_humidity = relative_humidity

        # If lat/long given instead of EarthLocation, convert them
        # to EarthLocation
        if location is None and (latitude is not None and
                                 longitude is not None):
            self.location = EarthLocation.from_geodetic(longitude, latitude,
                                                        elevation)

        elif isinstance(location, EarthLocation):
            self.location = location

        else:
            raise TypeError('Observatory location must be specified with '
                            'either (1) an instance of '
                            'astropy.coordinates.EarthLocation or (2) '
                            'latitude and longitude in degrees as '
                            'accepted by astropy.coordinates.Latitude and '
                            'astropy.coordinates.Latitude.')

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
        If `target` is None, generates an altitude/azimuth frame. Otherwise,
        calculates the transformation to that frame for the requested `target`.

        Parameters
        ----------
        time : `~astropy.time.Time`
            The time at which the observation is taking place. Will be used as
            the `obstime` attribute in the resulting frame or coordinate.

        target : `~astroplan.FixedTarget`, `~astropy.coordinates.SkyCoord`, defaults to `None` (optional)
            Celestial object of interest. If `target` is `None`, return the
            `~astropy.coordinates.AltAz` frame without coordinates.

        obswl : `~astropy.units.Quantity` (optional)
            Wavelength of the observation used in the calculation.

        Returns
        -------
        If `target` is `None`, returns `~astropy.coordinates.AltAz` frame. If
        `target` is not `None`, returns the `target` transformed to the
        `~astropy.coordinates.AltAz` frame.

        """

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

    def parallactic_angle(self, time, target):
        '''
        Calculate the parallactic angle.

        Parameters
        ----------
        time : `~astropy.time.Time`
            Observation time.

        target : `~astroplan.FixedTarget` or `~astropy.coordinates.SkyCoord`
            Celestial object of interest.

        Returns
        -------
        `~astropy.coordinates.Angle`
            Parallactic angle
        '''

        if not (hasattr(target, 'ra') or hasattr(target, 'dec') or
                hasattr(target, 'coord')):
            raise TypeError('The target must be a coordinate (i.e. a '
                            'FixedTarget or SkyCoord).')

        if hasattr(target, 'coord'):
            coordinate = target.coord
        else:
            coordinate = target

        # Eqn (14.1) of Meeus' Astronomical Algorithms

        H = (time.sidereal_time('mean', longitude=self.location.longitude)
             - coordinate.ra).radian
        q = np.arctan(np.sin(H) /
                      (np.tan(self.location.latitude.radian)*
                       np.cos(coordinate.dec.radian) -
                       np.sin(coordinate.dec.radian)*np.cos(H)))*u.rad

        if q < 0*u.rad:
            q += np.pi*u.rad
        return Angle(q)
    # Sun-related methods.

    def noon(self, time):
        """
        Returns the local, solar noon time.

        Parameters
        ----------
        time : `~astropy.time.Time`
        """
        raise NotImplementedError()

    def midnight(self, time):
        """
        Returns the local, solar midnight time.

        Parameters
        ----------
        time : `~astropy.time.Time`
        """
        raise NotImplementedError()

    def sunset(self, time, **kwargs):
        """
        Returns the local sunset time.

        The default sunset returned is the next one to occur.

        Parameters
        ----------
        time : `~astropy.time.Time`

        Keywords: str, optional
            previous
            next
        """
        raise NotImplementedError()

    def sunrise(self, time, **kwargs):
        """
        Returns the local sunrise time.

        The default sunrise returned is the next one to occur.

        Parameters
        ----------
        time : `~astropy.time.Time`

        Keywords: str, optional
            previous
            next
        """
        raise NotImplementedError()

    # Moon-related methods.

    def moonrise(self, time, **kwargs):
        """
        Returns the local moonrise time.

        The default moonrise returned is the next one to occur.

        Parameters
        ----------
        time : `~astropy.time.Time`

        Keywords: str, optional
            previous
            next
        """
        raise NotImplementedError()

    def moonset(self, time, **kwargs):
        """
        Returns the local moonset time.

        The default moonset returned is the next one to occur.

        Parameters
        ----------
        time : `~astropy.time.Time`

        Keywords: str, optional
            previous
            next
        """
        raise NotImplementedError()

    def moon_illumination(self, time):
        """
        Returns a float giving the percent illumation.

        Parameters
        ----------
        time : `~astropy.time.Time`
        """
        raise NotImplementedError()

    def moon_position(self, time):
        """
        Returns the position of the moon in alt/az.

        Parameters
        ----------
        time : `~astropy.time.Time`
        """
        raise NotImplementedError()

    # Other time-related methods.

    def evening_nautical(self, time):
        """
        Returns the local evening (nautical) twilight time.

        Parameters
        ----------
        time : `~astropy.time.Time`
        """
        raise NotImplementedError()

    def evening_civil(self, time):
        """
        Returns the local evening (civil) twilight time.

        Parameters
        ----------
        time : `~astropy.time.Time`
        """
        raise NotImplementedError()

    def morning_astronomical(self, time):
        """
        Returns the local morning (astronomical) twilight time.

        Parameters
        ----------
        time : `~astropy.time.Time`
        """
        raise NotImplementedError()

    def morning_nautical(self, time):
        """
        Returns the local morning (nautical) twilight time.

        Parameters
        ----------
        time : `~astropy.time.Time`
        """
        raise NotImplementedError()

    def morning_civil(self, time):
        """
        Returns the local morning (civil) twilight time.

        Parameters
        ----------
        time : `~astropy.time.Time`
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
    def __init__(self, coord, name=None, **kwargs):
        '''
        TODO: Docstring.
        '''
        if not (hasattr(coord, 'transform_to') and
                hasattr(coord, 'represent_as')):
            raise TypeError('`coord` must be a coordinate object.')

        self.name = name
        self.coord = coord

    @classmethod
    def from_name(cls, query_name, name=None, **kwargs):
        '''
        Initialize a `FixedTarget` by querying for a name, using the machinery
        in `~astropy.coordinates.SkyCoord.from_name`.
        '''
        # Allow manual override for name keyword so that the target name can
        # be different from the query name, otherwise assume name=queryname.
        if name is None:
            name = query_name
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

    def __init__(self, target, time):
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
