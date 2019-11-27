# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# Third-party
import numpy as np
from astropy import units as u


__all__ = ['exptime_from_ccd_snr']


@u.quantity_input(waveset=u.angstrom, flux=u.erg / u.s / u.cm ** 2 / u.cm,
                  npix=u.pixel, n_background=u.pixel,
                  background_rate=u.ct / u.pixel / u.s,
                  darkcurrent_rate=u.ct / u.pixel / u.s)
def exptime_from_ccd_snr(snr, waveset, flux, observer, telescope,
                         npix=1 * u.pixel,
                         n_background=np.inf * u.pixel,
                         background_rate=0 * (u.ct / u.pixel / u.s),
                         darkcurrent_rate=0 * (u.ct / u.pixel / u.s),
                         force_overlap='taper'):
    """
    Returns the exposure time needed in units of seconds to achieve
    the specified (idealized theoretical) signal to noise ratio
    (from pg 57-58 of [1]_).

    Parameters
    ----------
    snr : float, int, or `~astropy.units.Quantity`
        The signal to noise ratio of the given observation in dimensionless
        units.
    waveset : array-like
        The wavelengths associated with the target's flux.
    flux : array-like
        The flux of the target.
    observer : `~astroplan.observer.Observer`
        The Observer object.
    telescope : `~astroplan.telescope.Telescope`
        The Telescope object.
    npix : `~astropy.units.Quantity`, optional
        Number of pixels under consideration for the signal with units of
        pixels. Default is 1 * astropy.units.pixel.
    n_background : `~astropy.units.Quantity`, optional
        Number of pixels used in the background estimation with units of
        pixels. Default is set to np.inf * astropy.units.pixel
        such that there is no contribution of error due to background
        estimation. This assumes that n_background will be >> npix.
    background_rate : `~astropy.units.Quantity`, optional
        Photons per pixel per second due to the backround/sky with units of
        counts/second/pixel.
        Default is 0 * (astropy.units.ct /
        astropy.units.second / astropy.units.pixel)
    darkcurrent_rate : `~astropy.units.Quantity`, optional
        Counts per pixel per second due to the dark current with units
        of counts/second/pixel.
        Default is 0 * (astropy.units.ct /
        astropy.units.second / astropy.units.pixel)
    force_overlap : {'taper', 'extrap', 'none', None}, optional
        Force the spectral element x source spectrum convolution by the
        specified method even when they don't fully overlap in wavelength.
        'taper' forces the incomplete spectral component to zero for its
        missing wavelengths, while 'extrap' attempts to extrapolate the
        incomplete part of the spectrum.
        Default is 'taper'

    References
    ----------
    .. [1] Howell, S. B. 2000, *Handbook of CCD Astronomy* (Cambridge, UK:
        Cambridge University Press)

    Returns
    -------
    t : `~astropy.units.Quantity`
        The exposure time needed (in seconds) to achieve the given signal
        to noise ratio.
    """
    # import synphot, which isn't a required package of astroplan
    from synphot.models import Empirical1D
    from synphot.spectrum import SourceSpectrum, SpectralElement
    from synphot.observation import Observation

    # set the source spectrum with synphot
    source_spec = SourceSpectrum(Empirical1D, points=waveset,
                                 lookup_table=flux)

    # set the spectral elements if given (quantum efficiency,skymodel,bandpass)
    qe = telescope.ccd_response
    skymodel = observer.skymodel

    spec_elements = _get_spectral_element(telescope.bandpass,
                                          SpectralElement, Empirical1D)
    if skymodel is not False:
        spec_elements *= _get_spectral_element(skymodel,
                                               SpectralElement, Empirical1D)
    if qe is not False:
        spec_elements *= _get_spectral_element(qe,
                                               SpectralElement, Empirical1D)

    # get the synphot observation object
    synphot_obs = Observation(source_spec, spec_elements, force=force_overlap)

    # make sure the gain is in the correct units
    gain = telescope.gain
    if gain.unit in (u.electron / u.adu, u.photon / u.adu):
        gain = gain.value * (u.ct / u.adu)
    elif gain.unit != u.ct / u.adu:
        raise u.UnitsError('gain must have units of (either '
                           'astropy.units.ct, astropy.units.electron, or '
                           'astropy.units.photon) / astropy.units.adu')

    # get the countrate from the synphot observation object
    countrate = synphot_obs.countrate(area=telescope.area) / telescope.gain

    # define counts to be in ADU, which are not technically convertible in
    # astropy.units:
    countrate = countrate.value * (u.ct / u.s)

    # necessary for units to work in countrate calculation:
    if not hasattr(snr, 'unit'):
        snr = snr * np.sqrt(1 * u.ct)
    readnoise = _get_shotnoise(telescope.readnoise)
    gain_err = _get_shotnoise(telescope.gain * telescope.ad_err)

    # solve t with the quadratic equation (pg. 57 of Howell 2000)
    A = countrate ** 2
    B = (-1) * snr ** 2 * (countrate + npix * (background_rate +
                                               darkcurrent_rate))
    C = (-1) * snr ** 2 * npix * readnoise ** 2

    t = (-B + np.sqrt(B ** 2 - 4 * A * C)) / (2 * A)

    if gain_err.value > 1 or np.isfinite(n_background.value):
        from scipy.optimize import fsolve
        # solve t numerically
        t = fsolve(_t_with_small_errs, t, args=(background_rate,
                                                darkcurrent_rate,
                                                gain_err, readnoise, countrate,
                                                npix, n_background))
        t = float(t) * u.s

    return t


def _get_spectral_element(spec_tuple, SpectralElement, Empirical1D):
    """Returns a synphot SpectralElement made from the given 2D array"""
    points, lookup_table = spec_tuple
    return SpectralElement(Empirical1D, points=points,
                           lookup_table=lookup_table)


def _get_shotnoise(detector_property):
    """
    Returns the shot noise (i.e. non-Poissonion noise) in the correct
    units. ``detector_property`` must be a Quantity.
    """
    # Ensure detector_property is in the correct units:
    detector_property = detector_property.to(u.ct / u.pixel)
    return detector_property.value * np.sqrt(1 * (u.ct / u.pixel))


def _t_with_small_errs(t, background_rate, darkcurrent_rate, gain_err,
                       readnoise, countrate, npix, n_background):
    """
    Returns the full expression for the exposure time including the
    contribution to the noise from the background and the gain.
    """
    if not hasattr(t, 'unit'):
        t = t * u.s

    detector_noise = (background_rate * t + darkcurrent_rate * t +
                      gain_err ** 2 + readnoise ** 2)
    radicand = countrate * t + (npix * (1 + npix / n_background) *
                                detector_noise)

    return countrate * t / np.sqrt(radicand)
