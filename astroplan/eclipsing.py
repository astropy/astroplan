# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np
import astropy.units as u
from astropy.time import Time

__all__ = ['PeriodicEvent', 'EclipsingSystem']


class PeriodicEvent(object):
    """
    A periodic event defined by an epoch and period.
    """
    @u.quantity_input(period=u.day, duration=u.day)
    def __init__(self, epoch, period, duration=None, name=None):
        """

        Parameters
        ----------
        epoch : `~astropy.time.Time`
            Time of event
        period : `~astropy.units.Quantity`
            Period of event
        duration : `~astropy.units.Quantity` (optional)
            Duration of event
        name : str (optional)
            Name of target/event
        Returns
        -------
        """
        self.epoch = epoch
        self.period = period
        self.name = name
        self.duration = duration

    def phase(self, time):
        """
        Phase of periodic event, on range [0, 1). For example, the phase
        could be an orbital phase for an eclipsing binary system.

        Parameters
        ----------
        time : `~astropy.time.Time`
            Evaluate the phase at this time or times

        Returns
        -------
        phase_array : `~numpy.ndarray`
            Phase at each ``time``, on range [0, 1)
        """
        return ((time - self.epoch).to(u.day).value %
                self.period.to(u.day).value) / self.period.to(u.day).value


class EclipsingSystem(PeriodicEvent):
    """
    Define some parameters for an eclipsing system, either an eclipsing binary
    or a transiting exoplanet.
    """
    def in_primary_eclipse(self, time):
        """
        Returns `True` when ``time`` is during a primary eclipse.

        Parameters
        ----------
        time : `~astropy.time.Time`
            Time to evaluate

        Returns
        -------
        in_eclipse : `~numpy.ndarray` or bool
            ``True`` if ``time`` is during primary eclipse
        """
        phases = self.phase(time)
        return ((phases < float(self.duration/self.period)/2) |
                (phases > 1 - float(self.duration/self.period)/2))

    def in_secondary_eclipse(self, time, secondary_eclipse_phase=0.5):
        """
        Returns `True` when ``time`` is during a secondary eclipse.

        Parameters
        ----------
        time : `~astropy.time.Time`
            Time to evaluate
        secondary_eclipse_phase : float (optional)
            Defines the phase at secondary eclipse. For circular orbits, the
            secondary eclipse should happen at phase 0.5 -- in general, for
            non-circular orbits the secondary eclipse phase will be something
            else. Default is 0.5 (circular orbit).
        Returns
        -------
        in_eclipse : `~numpy.ndarray` or bool
            ``True`` if ``time`` is during secondary eclipse
        """
        phases = self.phase(time)
        return ((phases < secondary_eclipse_phase + float(self.duration/self.period)/2) &
                (phases > secondary_eclipse_phase - float(self.duration/self.period)/2))

    def out_of_eclipse(self, time):
        """
        Returns `True` when ``time`` is not during primary or secondary eclipse.

        Parameters
        ----------
        time : `~astropy.time.Time`
            Time to evaluate

        Returns
        -------
        in_eclipse : `~numpy.ndarray` or bool
            ``True`` if ``time`` is not during primary or secondary eclipse
        """
        return np.logical_not(np.logical_or(self.in_primary_eclipse(time),
                                            self.in_secondary_eclipse(time)))

    def next_primary_eclipse_time(self, time, n_eclipses=1):
        """
        Time of the next primary eclipse after ``time``.

        Parameters
        ----------
        time : `~astropy.time.Time`
            Find the next primary eclipse after ``time``
        n_eclipses : int (optional)
            Return the times of eclipse for the next ``n_eclipses`` after
            ``time``. Default is 1.

        Returns
        -------
        primary_eclipses : `~astropy.time.Time`
            Times of the next ``n_eclipses`` primary eclipses after ``time``
        """
        eclipse_times = ((1-self.phase(time)) * self.period + time +
                         np.arange(n_eclipses) * self.period)
        return eclipse_times

    def next_secondary_eclipse_time(self, time, n_eclipses=1):
        """
        Time of the next secondary eclipse after ``time``.

        Parameters
        ----------
        time : `~astropy.time.Time`
            Find the next secondary eclipse after ``time``
        n_eclipses : int (optional)
            Return the times of eclipse for the next ``n_eclipses`` after
            ``time``. Default is 1.

        Returns
        -------
        secondary_eclipses : `~astropy.time.Time`
            Times of the next ``n_eclipses`` secondary eclipses after ``time``
        """
        phase = self.phase(time)
        if phase >= 0.5:
            next_eclipse_phase = 1.5
        else:
            next_eclipse_phase = 0.5
        eclipse_times = ((next_eclipse_phase - phase) * self.period + time +
                         np.arange(n_eclipses) * self.period)
        return eclipse_times

    def primary_ingress_egress_time(self, time, n_eclipses=1):
        """
        Calculate the times of ingress and egress for the next ``n_eclipses``
        primary eclipses after ``time``

        Parameters
        ----------
        time : `~astropy.time.Time`
            Find the next primary ingress and egress after ``time``
        n_eclipses : int (optional)
            Return the times of eclipse for the next ``n_eclipses`` after
            ``time``. Default is 1.

        Returns
        -------
        primary_eclipses : `~astropy.time.Time` of shape (``n_eclipses``, 2)
            Times of ingress and egress for the next ``n_eclipses`` primary
            eclipses after ``time``
        """
        next_mid_eclipses = self.next_primary_eclipse_time(time, n_eclipses=n_eclipses)
        next_ingresses = next_mid_eclipses - self.duration/2
        next_egresses = next_mid_eclipses + self.duration/2

        ing_egr = np.vstack([next_ingresses.utc.jd, next_egresses.utc.jd]).T

        return Time(ing_egr, format='jd', scale='utc')

    def secondary_ingress_egress_time(self, time, n_eclipses=1):
        """
        Calculate the times of ingress and egress for the next ``n_eclipses``
        secondary eclipses after ``time``

        Parameters
        ----------
        time : `~astropy.time.Time`
            Find the next primary ingress and egress after ``time``
        n_eclipses : int (optional)
            Return the times of eclipse for the next ``n_eclipses`` after
            ``time``. Default is 1.

        Returns
        -------
        primary_eclipses : `~astropy.time.Time` of shape (``n_eclipses``, 2)
            Times of ingress and egress for the next ``n_eclipses`` primary
            eclipses after ``time``
        """
        next_mid_eclipses = self.next_secondary_eclipse_time(time, n_eclipses=n_eclipses)
        next_ingresses = next_mid_eclipses - self.duration/2
        next_egresses = next_mid_eclipses + self.duration/2

        ing_egr = np.vstack([next_ingresses.utc.jd, next_egresses.utc.jd]).T

        return Time(ing_egr, format='jd', scale='utc')


