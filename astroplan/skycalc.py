# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# Standard library
import json
import os

# Third-party
from astropy.io import fits
import astropy.units as u

# Local
from .exceptions import SkyCalcError


__all__ = ['skymodel']


def skymodel(airmass=1.0, pwv_mode='pwv', season=0,
             time=0, pwv=3.5, msolflux=130.0,
             incl_moon='Y', moon_sun_sep=90.0,
             moon_target_sep=45.0, moon_alt=45.0,
             moon_earth_dist=1.0, incl_starlight='Y',
             incl_zodiacal='Y',
             ecl_lon=135.0, ecl_lat=90.0,
             incl_loweratm='Y', incl_upperatm='Y',
             incl_airglow='Y', incl_therm='N',
             therm_t1=0.0, therm_e1=0.0,
             therm_t2=0.0, therm_e2=0.0, therm_t3=0.0,
             therm_e3=0.0, vacair='vac', wmin=300.0,
             wmax=2000.0,
             wgrid_mode='fixed_wavelength_step',
             wdelta=0.1, wres=20000, lsf_type='none',
             lsf_gauss_fwhm=5.0, lsf_boxcar_fwhm=5.0,
             observatory='paranal'):
    """
    Returns the model atmospheric transmittance curve queried from the SkyCalc
    Sky Model Calculator. The default parameters used here are the default
    parameters provided by SkyCalc:
    https://www.eso.org/observing/etc/doc/skycalc/helpskycalccli.html

    Parameters
    ----------
    airmass : float within range [1.0, 3.0]
        Airmass. Alt and airmass are coupled through the plane parallel
        approximation airmass=sec(z), z being the zenith distance
        z=90°−Alt
    pwv_mode : {'pwv', 'season'}. default = 'pwv'
    season : {0,1,2,3,4,5,6}. default = 0
        Time of year if not in pwv mode.
        (0 = all year, 1 = dec/jan, 2 = feb/mar...)
    time : {0, 1, 2, 3}. default = 0
        0 = all year, 1, 2, 3 = third of night
    pwv : {-1.0, 0.5, 1.0, 1.5, 2.5, 3.5, 5.0, 7.5, 10.0, 20.0}. default = 3.5
        Precipitable Water Vapor.
    msolflux : float
        Monthly Averaged Solar Flux, s.f.u float > 0 (default is 130.0)
    incl_moon : {'Y', 'N'}. default = 'Y'
        Flag for inclusion of scattered moonlight.
        Moon coordinate constraints: |z – zmoon| ≤ ρ ≤ |z + zmoon| where
        ρ=moon/target separation, z=90°−target altitude and
        zmoon=90°−moon altitude.
    moon_sun_sep : float within range [0.0, 360.0]. default = 90.0
        Degrees of separation between Sun and Moon as seen from Earth
        (i.e. the "moon phase").
    moon_target_sep : float in range [0.0, 180.0]. defualt = 45.0
        Degrees of separation between the moon and the target.
    moon_alt : float in range [-90.0, 90.0] defualt = 45.0
        Moon Altitude over Horizon in degrees.
    moon_earth_dist : float in range [0.91, 1.08]. default = 1.0
        Moon-Earth Distance (mean=1).
    incl_starlight : {'Y', 'N'}. default = 'Y'
        Flag for inclusion of scattered starlight.
    incl_zodiacal : {'Y', 'N'}. default = 'Y'
        Flag for inclusion of zodiacal light.
    ecl_lon : float in range [-180.0, 180.0]. default = 135.0
        Heliocentric ecliptic in degrees.
    ecl_lat : float in range [-90.0, 90.0]. default = 90.0
        Ecliptic latitude in degrees.
    incl_loweratm : {'Y', 'N'}. default = 'Y'
        Flag for inclusion of molecular emission of lower atmosphere.
    incl_upperatm : {'Y', 'N'}. default = 'Y'
        Flag for inclusion of molecular emission of upper atmosphere.
    incl_airglow : {'Y', 'N'}. default = 'Y'
        Flag for inclusion of airglow continuum (residual continuum)
    incl_therm : {'Y', 'N'}. default = 'N'
        Flag for inclusion of instrumental thermal radiation.
        Note: This radiance component represents an instrumental effect.
        The emission is provided relative to the other model components.
        To obtain the correct absolute flux, an instrumental response curve
        must be applied to the resulting model spectrum.
        See section 6.2.4 in the SkyCalc documentation
        `here <http://localhost/observing/etc/doc/skycalc/'
        'The_Cerro_Paranal_Advanced_Sky_Model.pdf>`_.
    therm_t1, therm_t2, therm_t3 : float. default = 0.0
        Temperature in K
    therm_e1, therm_e2, therm_e3 : float in range [0.0, 1.0]. default = 0.0
    vacair : {'vac', 'air'}. default = 'vac'
        In regards to the wavelength grid.
    wmin : float in range [300.0, 30000.0]. default = 300.0
        Minimum wavelength (nm) in the wavelength grid.
        Must be < wmax.
    wmax : float in range [300.0, 30000.0]. default = 2000.0
        Maximum wavelength (nm) in the wavelength grid.
        Must be > wmin.
    wgrid_mode : {'fixed_spectral_resolution', 'fixed_wavelength_step',
        'user'}.
        default = 'fixed_wavelength_step'
        Mode of the wavelength grid.
    wdelta : float in range [0, 30000.0]. default = 0.1
        Wavelength sampling step dlam (i.e. nm/step)
    wres : int in range [0,1.0e6]. default = 20000
        lam/dlam where dlam is wavelength step.
    wgrid_user : list of floats
        default = [500.0, 510.0, 520.0, 530.0, 540.0, 550.0]
    lsf_type : {'none', 'Gaussian', 'Boxcar'}. default = 'none'
        Line spread function type for convolution.
    lsf_gauss_fwhm : float. default = 5.0
        Gaussian full-width half-max for line spread function
        wavelength bins.
    lsf_boxcar_fwhm : float. default = 5.0
        Boxcar full-width half-max for line spread function
        wavelength bins.
    observatory : {'paranal', 'lasilla', 'armazones'}.
        default = 'paranal'
        Observatory where observation takes place.

    Returns
    -------
    trans_waves, transmission: tuple of arrays of floats
        'trans_waves' is an array of wavelengths in angstroms (float),
        'transmission' is an array of fractional atmospheric
        transmittance (float).
    """
    params = locals()

    if params['observatory'] == 'lasilla':
        params['observatory'] = '2400'
    elif params['observatory'] == 'paranal':
        params['observatory'] = '2640'
    elif (params['observatory'] == '3060m' or
          params['observatory'] == 'armazones'):
        params['observatory'] = '3060'
    else:
        raise ValueError('Wrong Observatory name, please refer to the '
                         'skycalc_cli documentation.')

    # import packages not required by astroplan
    import requests

    # Use the bit from skycalc_cli which queries from the SkyCalc Sky Model
    server = 'http://etimecalret-001.eso.org'
    url = server + '/observing/etc/api/skycalc'

    response = requests.post(url, data=json.dumps(params))
    results = json.loads(response.text)

    status = results['status']
    tmpdir = results['tmpdir']
    tmpurl = server + '/observing/etc/tmp/' + tmpdir + '/skytable.fits'

    if status == 'success':
        try:
            response = requests.get(tmpurl, stream=True)
            data = response.content
        except requests.exceptions.RequestException as e:
            print(e, 'could not retrieve FITS data from server')
    else:
        raise SkyCalcError('HTML request failed. A custom parameter you ' +
                           'set is most likely not accepted by the ' +
                           'SkyCalc Calculator. json.loads() returns: ' +
                           '{}'.format(results))

    # Create a temporary file to write the binary results to
    tmp_data_file = './tmp_skycalc_data.fits'

    with open(tmp_data_file, 'wb') as f:
        f.write(data)

    hdu = fits.open(tmp_data_file)
    trans_waves = hdu[1].data["LAM"] * u.um  # wavelengths
    transmission = hdu[1].data["TRANS"] * u.Unit('')

    # Delete the file after reading from it
    os.remove(tmp_data_file)

    return trans_waves.to(u.angstrom), transmission
