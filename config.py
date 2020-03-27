# This script defines global configuration/constant

# gaia table to use to cross matching
gaia_table_name = "gaiadr2.gaia_source"

# path to allstar
allstar_path = "/yngve_data/sdss/dr16/apogee/spectro/aspcap/r12/l33/allStar-r12-l33.fits"
# allstar_path = "C:/sdss_local/dr16/apogee/spectro/aspcap/r12/l33/allStar-r12-l33.fits"  # path on my windows

# base path to apStar/asStar files
base_path = "/yngve_data/sdss/apogee/dr16/apogee/spectro/redux/r12/stars"

# file name for master continuum normalized spectra, can be a path if file is not stored in current dir
contspac_file_name = "contspec_dr16.fits"

# file name for gaia all column, can be a path if file is not stored in current dir
gaia_allcolumns_f = "apogeedr16_gaiadr2_xmatch_allcolumns.fits"

# file name for gaia selected columns but row matched to allstar, can be a path if file is not stored in current dir
gaia_rowmatch_f = "apogeedr16_gaiadr2_xmatch.fits"

# file name for astroNN abundances, can be a path if file is not stored in current dir
astronn_chem_f = "astroNN_apogee_dr16post_catalog.fits"

# file name for astroNN abundances, can be a path if file is not stored in current dir
astronn_dist_f = "apogee_dr16post_nn_dist.fits"

# neural network models folders name, can be a path if file is not stored in current dir
astronn_chem_model = "astroNN_0617_run001"  # https://github.com/henrysky/astroNN_spectra_paper_figures/
astronn_dist_model = "astroNN_constant_model_reduced"  # https://github.com/henrysky/astroNN_gaia_dr2_paper/
astronn_age_model = "ted"

# file name for astroNN APGOEE VAC, can be a path if file is not stored in current dir
astronn_apogee_vac_f = "apogee_astroNN-DR16.fits"
