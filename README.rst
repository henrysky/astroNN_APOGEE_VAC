astroNN APOGEE VAC codes
===========================

This repository contains codes required to generate astroNN APOGEE VAC

astroNN DR16 VAC: https://www.sdss.org/dr16/data_access/value-added-catalogs/?vac_id=the-astronn-catalog-of-abundances,-distances,-and-ages-for-apogee-dr16-stars

Requirement
---------------

``python>=3.6``, ``numpy``, ``scipy``, ``astropy``, ``astroquery``, ``galpy``, ``astroNN``, ``pyyaml``, ``statsmodels``, ``h5py``

Scripts
---------

Please review the configurations stored in the following following files before running any script.

-   | `config.py`_
    | Global configuration for the scripts
-   | `gaia_credential`_
    | Gaia login if you have any

.. _config.py: config.py
.. _gaia_credential: gaia_credential

Please run the scripts in ASCENDING ORDER

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

.. _1_continuum_norm.py: 1_continuum_norm.py
.. _2_gaia_xmatch.py: 2_gaia_xmatch.py
.. _3_astroNN_chem_dist_ages.py: 3_astroNN_chem_dist_ages.py
.. _4_orbital_parameters.py: 4_orbital_parameters.py
.. _5_compile_vac.py: 5_compile_vac.py

Major Authors
---------------

-  | **Henry Leung** - henrysky_
   | Contact Henry: henrysky.leung [at] utoronto.ca
-  | **Jo Bovy** - jobovy_
-  | **Ted Mackereth** - jmackereth_


.. _henrysky: https://github.com/henrysky
.. _jobovy: https://github.com/jobovy
.. _jmackereth: https://github.com/jmackereth
