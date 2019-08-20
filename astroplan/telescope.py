# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# Third-party
import numpy as np
import astropy.units as u
from astroplan.exceptions import UserInputError


__all__ = ["Telescope"]

AD_ERR_DEFAULT = np.sqrt(0.289) * (u.adu / u.pixel)


def _default_ccd_response():
    """
    A generic response function of a silicon detector from the
    table on the following page:
    https://www.apo.nmsu.edu/arc35m/Instruments/ARCTIC/#3p5
    """
    wl = [100., 309.45600215, 319.53901519, 329.99173871,
          339.70504127, 349.78805431, 379.53294279, 390.00376402,
          399.57293121, 409.3114413, 419.5961146, 429.14136695,
          439.52687038, 449.9795939, 459.69289647, 469.60785929,
          479.85892255, 489.5638226, 499.59835962, 509.99161922,
          519.5607864, 529.30446727, 539.85285015, 549.59976275,
          559.51472558, 569.94676599, 579.68075166, 589.59571448,
          599.9631202, 609.56008031, 619.65269622, 630.09581687,
          639.67467926, 649.50561697, 659.84070534, 669.7685951,
          679.50258077, 689.9637068, 699.78494931, 709.53555533,
          719.83463294, 729.72858948, 739.88431656, 749.9673296,
          759.66253445, 769.59042422, 780.37854306, 789.59647942,
          799.60677842, 810.07759966, 819.65646205, 829.52341053,
          839.82248813, 849.68943661, 859.77244965, 870.07152726,
          879.71760973, 889.81096433, 899.94245339, 909.73137855,
          919.69435573, 930.06545485, 939.64431724, 949.72733028,
          959.95862293, 969.6772918, 979.76030484, 990.05938245,
          1100.]
    response = [0., 70.30747425, 73.11442327, 75.26541144, 76.58538173,
                74.70200949, 77.10351031, 81.79762371, 84.528972,
                87.18136276, 89.0405892, 91.16038906, 92.49142617,
                93.66048522, 94.87582372, 95.64295585, 96.56602958,
                97.20551595, 98.01709918, 98.32588679, 98.26189462,
                98.79719419, 99.10133838, 99.09379281, 99.35788748,
                99.30100555, 99.14661175, 98.81209184, 98.74843825,
                98.30855134, 98.04495971, 97.98459522, 97.24513015,
                96.85779131, 96.40002722, 96.03435769, 95.87183789,
                95.09275863, 94.89671913, 94.49854563, 94.12126754,
                93.52139537, 92.57268607, 91.77633908, 90.30410094,
                88.68846298, 86.74856762, 84.85909033, 82.85400237,
                80.76562301, 78.18504085, 75.90628116, 72.63150731,
                68.93418199, 65.34249453, 61.79608045, 58.13799205,
                54.49429826, 49.87975185, 45.18326838, 41.57397462,
                36.82027063, 32.80603172, 28.7917928, 25.09446748,
                21.39714216, 18.4392819, 14.89286782, 0.]
    return np.array(wl) * u.nm, np.array(response) * u.Unit('') / 100


def _user_input_test(user_input, user_input_str):
    """
    Raise appropriate errors if the user doesn't input the
    correct array len or type.
    """

    if len(user_input) != 2:
        raise UserInputError(user_input_str +
                             " ccd_response must be "
                             "a numpy array of shape = (2, n)")

    elif (type(user_input[0]) is not u.Quantity
          or type(user_input[1]) is not u.Quantity):
        raise TypeError(user_input_str +
                        " must be a 2D-array of"
                        "`~astropy.units.Quantity` objects")


class Telescope(object):
    """A class to store telescope specifications"""
    @u.quantity_input(diameter=u.m, gain=u.ct / u.adu,
                      readnoise=u.ct / u.pixel, ad_err=u.adu / u.pixel)
    def __init__(self, diameter, bandpass, gain, ccd_response=False,
                 readnoise=0 * (u.ct / u.pixel), ad_err=AD_ERR_DEFAULT):
        """
        Parameters
        ----------
        diameter : `~astropy.units.Quantity`
            Diameter of the telescope aperture.
        bandpass : `synphot.spectrum.SpectralElement` or 2D-array
            The bandpass of the instrument/telescope. Can either be
            specified by a `~synphot.spectrum.SpectralElement` object
            or an array with dimensions 2 x n, where one row is the
            array of wavelengths associated with the bandpass, and the other
            row contains the fractional values of the bandpass.
        gain : `~astropy.units.Quantity`, optional
            Gain of the CCD with units of counts/ADU. Default is
            1 * astropy.units.ct / astropy.units.adu such that the
            contribution to the error due to the gain is assumed to be small.
        ccd_response : 'default' or 2D-array of `~astropy.units.Quantity`
            objects
            Note: the zeroth index of the given array must be the
            waveset of the response function, while the other
            must contain the values of the function.
        readnoise : `~astropy.units.Quantity`, optional
            Counts per pixel from the read noise with units of counts/pixel.
            Default is 0 * astropy.units.ct / astropy.units.pixel.
        ad_err : `~astropy.units.Quantity`, optional
            An estimate of the 1 sigma error within the A/D converter with
            units of adu/pixel. Default is set to
            sqrt(0.289) * astropy.units.adu / astropy.units.pixel, where
            sqrt(0.289) comes from the assumption that "for a charge
            level that is half way between two output ADU steps, there
            is an equal chance that it will be assigned to the lower
            or to the higher ADU value when converted to a digital number"
            (text from footnote on page 56 of [1]_, see also [2]_).

            References
            ----------
        .. [1] Howell, S. B. 2000, *Handbook of CCD Astronomy* (Cambridge,
            UK: Cambridge University Press)
        .. [2] Merline, W. & Howell, S. B. *A Realistic Model for
            Point-sources Imaged on Array Detectors: The Model and
            Initial Results*. ExA, 6:163 (1995)
        """
        self.diameter = diameter
        self.gain = gain
        self.ad_err = ad_err
        self.readnoise = readnoise
        self._bandpass = bandpass
        self._ccd_response = ccd_response

    @property
    def area(self):
        """Area of the telescope"""
        return np.pi * (self.diameter / 2) ** 2

    @property
    def bandpass(self):
        """Bandpass of the instrument"""
        from synphot.spectrum import SpectralElement

        if isinstance(self._bandpass, SpectralElement):
            waves = self._bandpass.waveset
            return waves, self._bandpass(waves)

        else:
            _user_input_test(self._bandpass, "bandpass")

        return self._bandpass

    @property
    def ccd_response(self):
        """Response function of the CCD, i.e. the quantum efficiency"""
        if self._ccd_response == 'default':
            return _default_ccd_response()

        elif self._ccd_response is False:
            return self._ccd_response

        else:
            _user_input_test(self._ccd_response, "ccd_response")

        return self._ccd_response[0], self._ccd_response[1]
