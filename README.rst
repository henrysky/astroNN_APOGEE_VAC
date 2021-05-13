astroNN APOGEE VAC codes
===========================

This repository contains codes required to generate astroNN APOGEE VAC

- astroNN DR16 VAC: https://www.sdss.org/dr16/data_access/value-added-catalogs/?vac_id=the-astronn-catalog-of-abundances,-distances,-and-ages-for-apogee-dr16-stars
- astroNN DR1 7VAC: To be released publicly

Requirement
---------------

``python>=3.6``, ``numpy``, ``scipy``, ``astropy``, ``astroquery``, ``galpy``, ``astroNN``, ``pyyaml``, ``statsmodels``, ``h5py``, ``gaia_tools``

``zero_point``: https://gitlab.com/icc-ub/public/gaiadr3_zeropoint

Scripts
---------

Please review the configurations stored in the following following files before running any script.

-   | `config.py`_
    | Global configuration for the scripts
-   | `gaia_credential`_
    | Gaia login if you have any
-   | `utils.py`_
    | Contains a few useful utility functions
-   | `aspcap_norm.py`_
    | A simple script to retrieve all ASPCAP normalized spectra and save to a single file

.. _config.py: config.py
.. _gaia_credential: gaia_credential
.. _utils.py: utils.py
.. _aspcap_norm.py: aspcap_norm.py

Please run the scripts in ASCENDING ORDER

-   | `0_dr14_reference_contspec.py`_
    | This script continuum normalizing the original data set that the models trained on (aka dr14)
-   | `1_continuum_norm.py`_
    | This script continuum normalizing all the APOGEE spectra and save them to a single fits
-   | `2_gaia_xmatch.py`_
    | This script cross matching APOGEE-Gaia to get all columns and save them to a single fits
-   | `3_astroNN_chem_dist_ages.py`_
    | This script uses neural network models to get abundances, distance, ages and save them to fits
-   | `4_orbital_parameters.py`_
    | This script uses galpy to generate orbital parameters and save them to fits
-   | `5_compile_vac.py`_
    | This script uses all the fits files we have generated and compile VAC
-   | `6_diagnostic_plots.py`_
    | This script make some summary plots for the VAC

.. _0_dr14_reference_contspec.py: 0_dr14_reference_contspec.py
.. _1_continuum_norm.py: 1_continuum_norm.py
.. _2_gaia_xmatch.py: 2_gaia_xmatch.py
.. _3_astroNN_chem_dist_ages.py: 3_astroNN_chem_dist_ages.py
.. _4_orbital_parameters.py: 4_orbital_parameters.py
.. _5_compile_vac.py: 5_compile_vac.py
.. _6_diagnostic_plots.py: 6_diagnostic_plots.py

External Data
---------------

- ``APOKASC_cat_v6.6.1.fits.zip``: https://trac.sdss.org/attachment/wiki/APOGEE2/APOKASC/Catalog/APOKASC_cat_v6.6.1.fits.zip
- ``kepler_low_metallicity_with_samples.fits``: Internal use only

DR17 VAC
----------

We have retrained our models for DR17 VAC, the training scripts and sanity checks notebook are under the folder ``dr17-VAC-notebooks``

The folder ``astroNN_0512_run002``, ``astroNN_gaia_dr17_model_3``, ``APOKASC2_BCNN_age_combined_dr17_4`` are new neural network models used in DR17 VAC

Major Authors
---------------

-  | **Henry Leung** - henrysky_
   | Contact Henry: henrysky.leung [at] utoronto.ca
-  | **Jo Bovy** - jobovy_
-  | **Ted Mackereth** - jmackereth_


.. _henrysky: https://github.com/henrysky
.. _jobovy: https://github.com/jobovy
.. _jmackereth: https://github.com/jmackereth
