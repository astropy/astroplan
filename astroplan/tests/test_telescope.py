# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import astropy.units as u
from ..telescope import Telescope
import pytest
import numpy as np
from astroplan.exceptions import UserInputError
from astropy.tests.helper import raises


class TestTelescope(object):
    """Tests for telescope.py"""
    @pytest.mark.skipif(reason="synphot not installed")
    def setup_class(self):

        from synphot.spectrum import SpectralElement
        from astropy.utils.data import download_file

        svo_link = ('http://svo2.cab.inta-csic.es/theory/fps3/fps.php?ID=')
        filt_path = download_file(svo_link + 'SLOAN/SDSS.rprime_filter')

        self.diameter = 3.5 * u.m
        self.bandpass = SpectralElement.from_file(filt_path)
        self.gain = 1.9 * (u.ct / u.adu)
        self.default_readnoise = 0 * (u.ct / u.pixel)
        self.default_ad_err = np.sqrt(0.289) * (u.adu / u.pixel)
        self.telescope = Telescope(self.diameter, self.bandpass, self.gain)

    def test_diameter(self):
        """
        Tests that Telescope returns the correct diameter
        """
        assert(self.telescope.diameter == self.diameter)

    def test_bandpass(self):
        """
        Tests that Telescope returns the correct bandpass
        """
        waves = self.bandpass.waveset
        assert(all(self.telescope.bandpass[0] == waves) and
               all(self.telescope.bandpass[1] == self.bandpass(waves)))

    def test_gain(self):
        """
        Tests that Telescope returns the correct gain
        """
        assert(self.telescope.gain == self.gain)

    def test_ccd_response(self):
        """
        Tests that Telescope.ccd_response is False if a ccd response
        function isn't given.
        """
        assert(self.telescope.ccd_response is False)

    def test_readnoise(self):
        """Tests that Telescope returns the default readnoise"""
        assert(self.telescope.readnoise == self.default_readnoise)

    def test_ad_err(self):
        """Tests that Telescope returns the default ad_err"""
        assert(self.telescope.ad_err == self.default_ad_err)

    @raises(TypeError)
    def test_diameter_type(self):
        """
        Tests that Telescope raises an error if the user doesn't
        gives the bandpass as a unit quantity
        """
        return Telescope(1, self.bandpass, self.gain).diameter

    @raises(TypeError)
    def test_bandpass_type(self):
        """
        Tests that Telescope raises an error if the user doesn't
        gives the bandpass as a unit quantity
        """
        return Telescope(self.diameter, np.array([[], []]), self.gain).bandpass

    @raises(UserInputError)
    def test_bandpass_len(self):
        """
       Tests that Telescope raises an error if the user doesn't
       give a bandpass of the correct shape
        """
        return Telescope(self.diameter, np.array([]), self.gain).bandpass

    @raises(TypeError)
    def test_gain_type(self):
        """
        Tests that Telescope raises an error if gain isn't given
        as a `~astropy.units.quantity.Quantity`.
        """
        return Telescope(self.diameter, self.bandpass, 1).gain

    def test_default_ccd_response(self):
        """
        Tests that Telescope.ccd_response returns the following array
        if set to 'default'.
        """
        wl = np.array([100., 309.45600215, 319.53901519, 329.99173871,
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
                       1100.]) * u.nm
        response = np.array([0., 70.30747425, 73.11442327, 75.26541144,
                             76.58538173, 74.70200949, 77.10351031,
                             81.79762371, 84.528972, 87.18136276,
                             89.0405892, 91.16038906, 92.49142617,
                             93.66048522, 94.87582372, 95.64295585,
                             96.56602958, 97.20551595, 98.01709918,
                             98.32588679, 98.26189462, 98.79719419,
                             99.10133838, 99.09379281, 99.35788748,
                             99.30100555, 99.14661175, 98.81209184,
                             98.74843825, 98.30855134, 98.04495971,
                             97.98459522, 97.24513015, 96.85779131,
                             96.40002722, 96.03435769, 95.87183789,
                             95.09275863, 94.89671913, 94.49854563,
                             94.12126754, 93.52139537, 92.57268607,
                             91.77633908, 90.30410094, 88.68846298,
                             86.74856762, 84.85909033, 82.85400237,
                             80.76562301, 78.18504085, 75.90628116,
                             72.63150731, 68.93418199, 65.34249453,
                             61.79608045, 58.13799205, 54.49429826,
                             49.87975185, 45.18326838, 41.57397462,
                             36.82027063, 32.80603172, 28.7917928,
                             25.09446748, 21.39714216, 18.4392819,
                             14.89286782, 0.]) * u.Unit('') / 100
        tele = Telescope(self.diameter, self.bandpass, self.gain,
                         ccd_response='default')
        assert(all(wl == tele.ccd_response[0]) and
               all(response == tele.ccd_response[1]))

    def test_custom_ccd_response(self):
        """
        Tests that Telescope.ccd_response gives back the user specified
        ccd response function.
        """
        wl = np.arange(10) * u.angstrom
        response = np.arange(10) * u.Unit('')
        custom_ccd_response = (wl, response)

        tele = Telescope(self.diameter, self.bandpass, self.gain,
                         ccd_response=custom_ccd_response)
        assert(tele.ccd_response == custom_ccd_response)

    @raises(TypeError)
    def test_custom_ccd_response_type(self):
        """
        Tests that telescope raises an error if ccd_response isn't passed in
        as a `~astropy.units.quantity.Quantity`.
        """
        wl = np.arange(10)
        response = np.arange(10)
        custom_ccd_response = (wl, response)
        return Telescope(self.diameter, self.bandpass, self.gain,
                         ccd_response=custom_ccd_response).ccd_response

    @raises(UserInputError)
    def test_custom_ccd_response_len(self):
        """
        Tests that telescope raises an error if the given ccd_response isn't
        the right length
        """
        return Telescope(self.diameter, self.bandpass, self.gain,
                         ccd_response=np.array([])).ccd_response

    @raises(TypeError)
    def test_readnoise_type(self):
        """Tests that Telescope raises an error if readnoise isn't passed in
        as a `~astropy.units.quantity.Quantity`.
        """
        return Telescope(self.diameter, self.bandpass, self.gain,
                         readnoise=0)

    @raises(TypeError)
    def test_ad_err_type(self):
        """
        Tests that telescope raises an error if ad_err isn't given
        as a `~astropy.units.quantity.Quantity`.
        """
        return Telescope(self.diameter, self.bandpass, self.gain,
                         ad_err=1).ad_err
