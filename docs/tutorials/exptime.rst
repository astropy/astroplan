.. _exptime_tutorial:

.. doctest-skip-all

****************************
Predicting an Exposure Time
****************************

This module of astroplan requires that `scipy <https://www.scipy.org/>`_
and `synphot <https://synphot.readthedocs.io/en/latest/>`_ be installed.

Contents
========

 In this tutorial we will show how to predict the exposure time needed
 to obtain a given signal to noise ratio for a target star.

* :ref:`exptime-querying_target_properties`
* :ref:`exptime-obtaining_a_model_spectrum`
* :ref:`exptime-setting_the_location_of_the_observer`
* :ref:`exptime-creating_a_telescope_model`
* :ref:`exptime-getting_the_exposure_time`
* :ref:`exptime-adding_an_atmospheric_transmission_model`

.. code-block:: python

    >>> import astropy.units as u

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt

.. _exptime-querying_target_properties:

Querying Target Properties
==========================

In this example we will use HAT-P-11 as the target. We will model the
spectrum of HAT-P-11 with a PHOENIX stellar atmosphere model, which we
can obtain by specifying a stellar temperature and querying from `the
Gottingen Spectral Library. <http://phoenix.astro.physik.uni-goettingen.de/>`_

To get the stellar temperature, we can use :meth:`astroplan.FixedTarget.from_name`
to get the coordinates of HAT-P-11, then give these coordinates to
`astroquery.gaia`:

.. code-block:: python

    >>> from astroplan import FixedTarget
    >>> from astroquery.gaia import Gaia

    >>> hatp11 = FixedTarget.from_name('HAT-P-11')

    >>> # width / height of search:
    >>> width = u.Quantity(1, u.arcmin)
    >>> height = u.Quantity(1, u.arcmin)
    >>> search_results = Gaia.query_object_async(coordinate=hatp11.coord,
                                                 width=width, height=height)

    >>> # the queried star should be the one nearest to the given coordinates
    >>> search_results.add_index('dist', unique=True)
    >>> hatp11_info = search_results.loc['dist', min(search_results['dist'])]

    >>> T_eff = round(hatp11_info['teff_val'], -2)  # round to nearest 100 K

Since the spectrum is defined at the source, we have to scale for the
distance to the target to get the flux at the telescope:

.. code-block:: python

    >>> stellar_radius = hatp11_info['radius_val'] * u.R_sun
    >>> parallax = hatp11_info['parallax'] * u.mas
    >>> # parallax is given by Gaia in milliarcseconds, so convert to parsec:
    >>> distance = parallax.to(u.parsec, equivalencies=u.parallax())

    >>> distance_scale = float(stellar_radius / distance) ** 2 / np.pi

.. _exptime-obtaining_a_model_spectrum:

Obtaining a Model Spectrum
==========================

If you don't already have a model spectrum, you may query one (at least for
a star) from the Gottingen Spectral Library using `astropy.io.fits.getdata`:

.. code-block:: python

    >>> from astropy.io import fits

    >>> flux_url = ('ftp://phoenix.astro.physik.uni-goettingen.de/'
                    'v2.0/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/'
                    'lte{T_eff:05d}-{log_g:1.2f}-0.0.PHOENIX-'ACES-AGSS-'
                    'COND-2011-HiRes.fits').format(T_eff=int(T_eff),
                    log_g=4.5)
    >>> wavelength_url = ('ftp://phoenix.astro.physik.uni-goettingen.de/'
                          'v2.0/HiResFITS/WAVE_PHOENIX-ACES-AGSS-COND-'
                          '2011.fits')

    >>> flux_of_target = fits.getdata(flux_url) * distance_scale
    >>> waveset_of_target = fits.getdata(wavelength_url)

The units we use below are specified by the PHOENIX models.

.. note::
    It is good practice (and for some packages, essential) to attach units to
    quantities so that you can be sure the final results are scaled correctly.

.. code-block:: python

    >>> flux_of_target = flux_of_target * (u.erg / u.s / u.cm ** 2 / u.cm)
    >>> waveset_of_target = waveset_of_target * u.Angstrom

    >>> plt.plot(waveset_of_target, flux_of_target)

.. plot::

    import astropy.units as u

    import numpy as np
    import matplotlib.pyplot as plt

    from astroplan import FixedTarget
    from astroquery.gaia import Gaia

    hatp11 = FixedTarget.from_name('Hat-p-11')

    # width / height of search:
    width = u.Quantity(1, u.arcmin)
    height = u.Quantity(1, u.arcmin)
    search_results = Gaia.query_object_async(coordinate=hatp11.coord, width=width, height=height)

    # the queried star should be the one nearest to the given coordinates
    search_results.add_index('dist', unique=True)
    hatp11_info = search_results.loc['dist', min(search_results['dist'])]

    T_eff = round(hatp11_info['teff_val'], -2)  # round to nearest 100 K

    stellar_radius = hatp11_info['radius_val'] * u.R_sun
    parallax = hatp11_info['parallax'] * u.mas
    # parallax is given by Gaia in milliarcseconds, so convert to parsec:
    distance = parallax.to(u.parsec, equivalencies=u.parallax())

    distance_scale = float(stellar_radius / distance) ** 2 / np.pi

    from astropy.io import fits

    flux_url = ('ftp://phoenix.astro.physik.uni-goettingen.de/v2.0/HiResFITS/'
           'PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte{T_eff:05d}-{log_g:1.2f}-0.0.PHOENIX-'
           'ACES-AGSS-COND-2011-HiRes.fits').format(T_eff=int(T_eff), log_g=4.5)
    wavelength_url = ('ftp://phoenix.astro.physik.uni-goettingen.de/v2.0/HiResFITS/'
                      'WAVE_PHOENIX-ACES-AGSS-COND-2011.fits')

    flux_of_target = fits.getdata(flux_url) * distance_scale
    waveset_of_target = fits.getdata(wavelength_url)

    flux_of_target = flux_of_target * (u.erg / u.s / u.cm ** 2 / u.cm)
    waveset_of_target = waveset_of_target * u.Angstrom

    plt.plot(waveset_of_target, flux_of_target)
    plt.show()

.. _exptime-setting_the_location_of_the_observer:

Setting the Location of the Observer
====================================

Create an :class:`astroplan.observer.Observer` object via the
:meth:`~astroplan.observer.Observer.at_site` method to set the
location of the observing site:

.. code-block:: python

    >>> from astroplan import Observer

    >>> observer = Observer.at_site('apo')

.. _exptime-creating_a_telescope_model:

Creating a Telescope Model
==========================

Create a :class:`astroplan.telescope.Telescope` object, which models the telescope
you wish to use (we model ours off of `the ARCTIC instrument <https://www.apo.nmsu.
edu/arc35m/Instruments/ARCTIC/>`_ on APO's 3.5m telescope), and input:

- the diameter of the telescope aperture
- a model of your bandpass (we use `synphot <https://synphot.readthedocs.io/en/latest/>`_
  to generate a bandpass downloaded from `the Spanish Virtual Observatory <http://svo2.
  cab.inta-csic.es/theory/fps/index.php?mode=browse>`_)
- the gain of the instrument, which is usually found on the observatory's website

.. code-block:: python

    >>> from synphot.spectrum import SpectralElement
    >>> from astropy.utils.data import download_file

    >>> diameter = 3.5 * u.m
    >>> gain = 1.9 * (u.ct / u.adu)

    >>> # model the bandpass:
    >>> svo_link = ('http://svo2.cab.inta-csic.es/' +
                    'theory/fps3/fps.php?ID=')
    >>> filt_path = download_file(svo_link + 'SLOAN/SDSS.rprime_filter')
    >>> bandpass = SpectralElement.from_file(filt_path)
    >>> bandpass.plot()

.. plot::

    import astropy.units as u

    import numpy as np
    import matplotlib.pyplot as plt

    from astroplan import FixedTarget
    from astroquery.gaia import Gaia

    hatp11 = FixedTarget.from_name('Hat-p-11')

    # width / height of search:
    width = u.Quantity(1, u.arcmin)
    height = u.Quantity(1, u.arcmin)
    search_results = Gaia.query_object_async(coordinate=hatp11.coord, width=width, height=height)

    # the queried star should be the one nearest to the given coordinates
    search_results.add_index('dist', unique=True)
    hatp11_info = search_results.loc['dist', min(search_results['dist'])]

    T_eff = round(hatp11_info['teff_val'], -2)  # round to nearest 100 K

    stellar_radius = hatp11_info['radius_val'] * u.R_sun
    parallax = hatp11_info['parallax'] * u.mas
    # parallax is given by Gaia in milliarcseconds, so convert to parsec:
    distance = parallax.to(u.parsec, equivalencies=u.parallax())

    distance_scale = float(stellar_radius / distance) ** 2 / np.pi

    from astropy.io import fits

    flux_url = ('ftp://phoenix.astro.physik.uni-goettingen.de/v2.0/HiResFITS/'
           'PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte{T_eff:05d}-{log_g:1.2f}-0.0.PHOENIX-'
           'ACES-AGSS-COND-2011-HiRes.fits').format(T_eff=int(T_eff), log_g=4.5)
    wavelength_url = ('ftp://phoenix.astro.physik.uni-goettingen.de/v2.0/HiResFITS/'
                      'WAVE_PHOENIX-ACES-AGSS-COND-2011.fits')

    flux_of_target = fits.getdata(flux_url) * distance_scale
    waveset_of_target = fits.getdata(wavelength_url)

    flux_of_target = flux_of_target * (u.erg / u.s / u.cm ** 2 / u.cm)
    waveset_of_target = waveset_of_target * u.Angstrom

    from astroplan import Observer

    observer = Observer.at_site('apo')

    from synphot.spectrum import SpectralElement
    from astropy.utils.data import download_file
    from astroplan import Telescope

    diameter = 3.5 * u.m
    gain = 1.9 * (u.ct / u.adu)

    # model the bandpass:
    svo_link = ('http://svo2.cab.inta-csic.es/' +
               'theory/fps3/fps.php?ID=')
    filt_path = download_file(svo_link + 'SLOAN/SDSS.rprime_filter')
    bandpass = SpectralElement.from_file(filt_path)

    bandpass.plot()

You can also give it a CCD response function, either from your own model or
the default model given by :class:`~astroplan.telescope.Telescope` as we've done here:

.. code-block:: python

    >>> from astroplan import Telescope

    >>> telescope = Telescope(diameter, bandpass, gain, ccd_response='default')

.. code-block:: python

    >>> ccd_response_wls, ccd_response = telescope.ccd_response
    >>> plt.plot(ccd_response_wls, ccd_response)
    >>> plt.xlabel(ccd_response_wls.unit)

.. plot::

    import astropy.units as u

    import numpy as np
    import matplotlib.pyplot as plt

    from astroplan import FixedTarget
    from astroquery.gaia import Gaia

    hatp11 = FixedTarget.from_name('HAT-P-11')

    # width / height of search:
    width = u.Quantity(1, u.arcmin)
    height = u.Quantity(1, u.arcmin)
    search_results = Gaia.query_object_async(coordinate=hatp11.coord, width=width, height=height)

    # the queried star should be the one nearest to the given coordinates
    search_results.add_index('dist', unique=True)
    hatp11_info = search_results.loc['dist', min(search_results['dist'])]

    T_eff = round(hatp11_info['teff_val'], -2)  # round to nearest 100 K

    stellar_radius = hatp11_info['radius_val'] * u.R_sun
    parallax = hatp11_info['parallax'] * u.mas
    # parallax is given by Gaia in milliarcseconds, so convert to parsec:
    distance = parallax.to(u.parsec, equivalencies=u.parallax())

    distance_scale = float(stellar_radius / distance) ** 2 / np.pi

    from astropy.io import fits

    flux_url = ('ftp://phoenix.astro.physik.uni-goettingen.de/v2.0/HiResFITS/'
           'PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte{T_eff:05d}-{log_g:1.2f}-0.0.PHOENIX-'
           'ACES-AGSS-COND-2011-HiRes.fits').format(T_eff=int(T_eff), log_g=4.5)
    wavelength_url = ('ftp://phoenix.astro.physik.uni-goettingen.de/v2.0/HiResFITS/'
                      'WAVE_PHOENIX-ACES-AGSS-COND-2011.fits')

    flux_of_target = fits.getdata(flux_url) * distance_scale
    waveset_of_target = fits.getdata(wavelength_url)

    flux_of_target = flux_of_target * (u.erg / u.s / u.cm ** 2 / u.cm)
    waveset_of_target = waveset_of_target * u.Angstrom

    from astroplan import Observer

    observer = Observer.at_site('apo')

    from synphot.spectrum import SpectralElement
    from astropy.utils.data import download_file
    from astroplan import Telescope

    diameter = 3.5 * u.m
    gain = 1.9 * (u.ct / u.adu)

    # model the bandpass:
    svo_link = ('http://svo2.cab.inta-csic.es/' +
               'theory/fps3/fps.php?ID=')
    filt_path = download_file(svo_link + 'SLOAN/SDSS.rprime_filter')
    bandpass = SpectralElement.from_file(filt_path)

    telescope = Telescope(diameter, bandpass, gain, ccd_response='default')
    ccd_response_wls, ccd_response = telescope.ccd_response
    plt.plot(ccd_response_wls, ccd_response)
    plt.xlabel(ccd_response_wls.unit)
    plt.show()

.. _exptime-getting_the_exposure_time:

Getting the Exposure Time
=========================

Then specify the signal to noise ratio you wish to reach and
call the :func:`~astroplan.exptime_from_ccd_snr()` function:

.. code-block:: python

    >>> from astroplan import exptime_from_ccd_snr

    >>> snr = 100
    >>> exptime_from_ccd_snr(snr, waveset_of_target, flux_of_target, observer, telescope)
    0.0022556432s

.. note::
    SNR is derived from the equations on pages 57-58 of Howell (2000). See the
    docstring of :func:`~astroplan.exptime_from_ccd_snr()` for more details on
    the reference and optional input parameters.

.. _exptime-adding_an_atmospheric_transmission_model:

(optional) Adding an Atmospheric Transmission Model
===================================================

Say you want to add the effect of the atmosphere to make the exposure
time prediction more realistic. You can pass your own sky model as an
array into your observer object, or you can query from the `SKYCALC
Sky Model Calculator <https://www.eso.org/observing/etc/bin/gen/form?
INS.MODE=swspectr+INS.NAME=SKYCALC>`_ with :func:`astroplan.skymodel()`:

.. code-block:: python

    >>> from astroplan import skymodel

    >>> skymodel_1pt5airmass = skymodel(airmass=1.5)

`skymodel` has many different sky parameters that can be set, including
the precipitable water vapor, whether or not to include the moon, and more.
For options and defaults, see :func:`astroplan.skymodel()`.

The skymodel attribute of the `Observer` object must then be set:

.. code-block:: python

    >>> observer.skymodel = skymodel_1pt5airmass

    >>> skymodel_waveset, skymodel_flux = observer.skymodel
    >>> plt.plot(skymodel_waveset, skymodel_flux)

.. note::

    Alternatively, you could pass it in as an optional argument
    when initiating the observer object:

    .. code-block:: python

        >>> observer = Observer.at_site('apo', skymodel=skymodel_1pt5airmass)

.. plot::

    import astropy.units as u

    import numpy as np
    import matplotlib.pyplot as plt

    from astroplan import FixedTarget
    from astroquery.gaia import Gaia

    hatp11 = FixedTarget.from_name('Hat-p-11')

    # width / height of search:
    width = u.Quantity(1, u.arcmin)
    height = u.Quantity(1, u.arcmin)
    search_results = Gaia.query_object_async(coordinate=hatp11.coord, width=width, height=height)

    # the queried star should be the one nearest to the given coordinates
    search_results.add_index('dist', unique=True)
    hatp11_info = search_results.loc['dist', min(search_results['dist'])]

    T_eff = round(hatp11_info['teff_val'], -2)  # round to nearest 100 K

    stellar_radius = hatp11_info['radius_val'] * u.R_sun
    parallax = hatp11_info['parallax'] * u.mas
    # parallax is given by Gaia in milliarcseconds, so convert to parsec:
    distance = parallax.to(u.parsec, equivalencies=u.parallax())

    distance_scale = float(stellar_radius / distance) ** 2 / np.pi

    from astropy.io import fits

    flux_url = ('ftp://phoenix.astro.physik.uni-goettingen.de/v2.0/HiResFITS/'
           'PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte{T_eff:05d}-{log_g:1.2f}-0.0.PHOENIX-'
           'ACES-AGSS-COND-2011-HiRes.fits').format(T_eff=int(T_eff), log_g=4.5)
    wavelength_url = ('ftp://phoenix.astro.physik.uni-goettingen.de/v2.0/HiResFITS/'
                      'WAVE_PHOENIX-ACES-AGSS-COND-2011.fits')

    flux_of_target = fits.getdata(flux_url) * distance_scale
    waveset_of_target = fits.getdata(wavelength_url)

    flux_of_target = flux_of_target * (u.erg / u.s / u.cm ** 2 / u.cm)
    waveset_of_target = waveset_of_target * u.Angstrom

    from astroplan import Observer

    observer = Observer.at_site('apo')

    from synphot.spectrum import SpectralElement
    from astropy.utils.data import download_file
    from astroplan import Telescope

    diameter = 3.5 * u.m
    gain = 1.9 * (u.ct / u.adu)

    # model the bandpass:
    svo_link = ('http://svo2.cab.inta-csic.es/' +
               'theory/fps3/fps.php?ID=')
    filt_path = download_file(svo_link + 'SLOAN/SDSS.rprime_filter')
    bandpass = SpectralElement.from_file(filt_path)

    telescope = Telescope(diameter, bandpass, gain, ccd_response='default')
    ccd_response_wls, ccd_response = telescope.ccd_response

    from astroplan import exptime_from_ccd_snr

    snr = 100
    exptime_from_ccd_snr(snr, waveset_of_target, flux_of_target, observer, telescope)

    from astroplan import skymodel

    skymodel_1pt5airmass = skymodel(airmass=1.5)

    observer.skymodel = skymodel_1pt5airmass
    # alternatively: >>> observer = Observer.at_site('apo', skymodel=skymodel_1pt5airmass)
    skymodel_waveset, skymodel_flux = observer.skymodel
    plt.plot(skymodel_waveset, skymodel_flux)
    plt.show()

And get the exposure time:

.. code-block:: python

    >>> exptime_from_ccd_snr(snr, waveset_of_target, flux_of_target, observer, telescope)
    0.0026073144s

.. note::
    The exposure time is slightly longer when accounting for the atmosphere than it was without an atmosphere due to slight atmospheric extinction in r' the band.
