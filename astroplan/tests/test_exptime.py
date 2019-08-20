# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# Third-party
import astropy.units as u
from astropy.tests.helper import assert_quantity_allclose
import numpy as np
from astroplan.exceptions import UserInputError
from astropy.tests.helper import raises
import astropy.units as u
import numpy as np
from astroplan import FixedTarget, Observer, Telescope
from astroquery.gaia import Gaia
from astropy.utils.data import download_file

# Local
from astroplan import exptime_from_ccd_snr

@raises(TypeError)
def _raise_type_error_test(testcase):
    return testcase

class TestExptime(object):
    """Tests for exptime.py"""
    @pytest.mark.skipif(reason="synphot not installed")
    def setup_class(self):

        from synphot.spectrum import SpectralElement
        from synphot.models import ConstFlux1D

    def test_synphot_specmodel(self):
        """
        Tests that exptime returns the expected result given a synphot
        spectral model, for both numeric and analytic solutions.
        This test is built off of the exptime tutorial in the
        documentation.
        """
        # numeric, i.e. when n_background != inf or gain_err > 1
        # note that gain_err > 1 by default, since gain_err = gain * ad_err, where
        # ad_err is >1 by default
        snr = 1

        waveset = np.array([1, 2]) * u.angstrom
        flux = SpectralElement(ConstFlux1D, 1 * (u.photon / u.s / u.cm ** 2 / u.angstrom))

        area = 1 * u.cm ** 2
        diameter = np.sqrt(area / np.pi) * 2
        gain = 1 * u.photon / u.adu
        throughput = np.ones(len(waveset)) * u.Unit('')
        bandpass = SpectralElement(Empirical1D, points=waveset, lookup_table=throughput)
        ad_err = 0 * u.adu / u.pixel
        telescope = Telescope(diameter, bandpass, gain, ad_err=ad_err)

        result = exptime_from_ccd_snr(snr, waveset, flux, observer, telescope)
        assert_quantity_allclose(result, 1 * u.s, rtol=1e-3)

        # numeric, i.e. when n_background != inf or gain_err > 1
        telescope = Telescope(diameter, bandpass, gain)
        # gain error doesn't have that big of an effect so the result shouldn't
        # be too different
        assert_quantity_allclose(result, 1 * u.s, rtol=1e-1)

    def test_user_specmodel(self):
        """
        Tests that exptime returns the expected result given a user
        inputed spectral model, for both numeric and analytic solutions
        """
        # numeric, i.e. when n_background != inf or gain_err > 1
        # note that gain_err > 1 by default, since gain_err = gain * ad_err, where
        # ad_err is >1 by default
        T_eff = 4800.0
        distance_scale = 6.701456436998426e-20

        # setting up the synphot spectrum
        flux_url = ('ftp://phoenix.astro.physik.uni-goettingen.de/'
                'v2.0/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/'
                'lte{T_eff:05d}-{log_g:1.2f}-0.0.PHOENIX-ACES-AGSS-'
                'COND-2011-HiRes.fits').format(T_eff=int(T_eff),
                log_g=4.5)
        wavelength_url = ('ftp://phoenix.astro.physik.uni-goettingen.de/'
                          'v2.0/HiResFITS/WAVE_PHOENIX-ACES-AGSS-COND-'
                          '2011.fits')
        flux_of_target = fits.getdata(flux_url) * distance_scale
        waveset_of_target = fits.getdata(wavelength_url)

        # setting up the telescope
        diameter = 3.5 * u.m
        svo_link = ('http://svo2.cab.inta-csic.es/theory/fps3/fps.php?ID=')
        filt_path = download_file(svo_link + 'SLOAN/SDSS.rprime_filter')
        bandpass = SpectralElement.from_file(filt_path)
        gain = 1.9 * (u.ct / u.adu)

        snr = 100
        flux = flux_of_target * (u.erg / u.s / u.cm ** 2 / u.cm)
        waveset = waveset_of_target * u.Angstrom
        observer = Observer.at_site('apo')
        telescope = Telescope(diameter, bandpass, gain, ccd_response='default')
        result = exptime_from_ccd_snr(snr, waveset, flux, observer, telescope)
        assert_quantity_allclose(result, 0.0022556432 * u.s, rtol=1e-3)

        # analytic
        waveset = np.array([1, 2]) * u.angstrom
        throughput = np.ones(len(waveset)) * u.Unit('')
        bandpass = SpectralElement(Empirical1D, points=waveset, lookup_table=throughput)
        flux = np.ones(len(waveset)) * (u.photon / u.s / u.cm ** 2 / u.angstrom)
        ad_err = 0 * u.adu / u.pixel

    @raises(TypeError)
    def test_user_specmodel_typeerror(self):
        """
        Tests that exptime raises an error when the user gives a spectral
        model of the wrong type (ie not an array of
        `~astropy.units.quantity.Quantity` objects)
        """

    @raises(UserInputError)
    def test_user_specmodel_lenerror(self):
        """
        Tests that exptime raises an error when the user gives a spectral
        model of the wrong length
        """

    def test_optional_param_type_errors_raised(self):
        """
        Tests that exptime raises type errors whenever any optional
        parameters is given in the incorrect unit
        """
        _raise_type_error_test(testcase1)
        _raise_type_error_test(testcase2)
        _raise_type_error_test(testcase3)


def test_t_exp_numeric():
    """
    A test to check that the numerical method in
    observation.exptime_from_ccd_snr (i.e. when the error in background
    noise or gain is non-negligible) is done correctly. Based on the worked
    example on pg. 56-57 of [1]_.

    References
    ----------
    .. [1] Howell, S. B. 2000, *Handbook of CCD Astronomy* (Cambridge, UK:
        Cambridge University Press)
    """
    t = 300 * u.s
    snr = 342
    gain = 5 * (u.ct / u.adu)
    countrate = 24013 * u.adu * gain / t
    npix = 1 * u.pixel
    n_background = 200 * u.pixel
    background_rate = 620 * (u.adu / u.pixel) * gain / t
    darkcurrent_rate = 22 * (u.ct / u.pixel / u.hr)
    readnoise = 5 * (u.ct / u.pixel)

    result = exptime_from_ccd_snr(snr, countrate, npix=npix,
                                  n_background=n_background,
                                  background_rate=background_rate,
                                  darkcurrent_rate=darkcurrent_rate,
                                  readnoise=readnoise, gain=gain)

    assert_quantity_allclose(result, t, atol=1 * u.s)


def test_t_exp_analytic():
    """
    A test to check that the analytic method in
    observation.exptime_from_ccd_snr is done correctly.
    """
    snr_set = 50
    countrate = 1000 * (u.ct / u.s)
    npix = 1 * u.pixel
    background_rate = 100 * (u.ct / u.pixel / u.s)
    darkcurrent_rate = 5 * (u.ct / u.pixel / u.s)
    readnoise = 1 * (u.ct / u.pixel)

    t = exptime_from_ccd_snr(snr_set, countrate, npix=npix,
                             background_rate=background_rate,
                             darkcurrent_rate=darkcurrent_rate,
                             readnoise=readnoise)

    # if t is correct, ccd_snr() should return snr_set:
    snr_calc = ccd_snr(countrate * t,
                       npix=npix,
                       background=background_rate * t,
                       darkcurrent=darkcurrent_rate * t,
                       readnoise=readnoise)

    assert_quantity_allclose(snr_calc, snr_set,
                             atol=0.5 * u.Unit(""))