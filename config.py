# This script defines global configuration/constant

from astroNN.apogee import allstar

# milkyway potential parameters
_R0 = 8.125  # distance to galactic center in kpc
_v0 = 242.  # sun tangential velocity in km/x
_z0 = 0.0208  # sun above galactic plane in kpc

# gaia table to use to cross matching
gaia_table_name = "gaiadr3.gaia_source"

# path to the reference allstar (which is dr14 the models are trained on)
allstar14_path = allstar(dr=14)
# file name for dr14 continuum normalized spectra, can be a path if file is not stored in current dir
contspac14_file_name = "contspec_dr14.fits"
# base path to dr14apStar files
base14_path = "/yngve_data/sdss/apogee/dr14/apogee/spectro/redux/r8/stars"


# diff correction (calculate median diff to dr14, and subtract it from current dataset to make it looks more like dr14)
# will do north/south separately
corr_flag = True

# path to your current allstar
# allstar_path = "/yngve_data/sdss/apogee/apogeework/apogee/spectro/aspcap/r13/l33/allStar-r13-l33-58932beta.fits"
allstar_path = allstar(dr=17)  # path on my windows

# base path to apStar/asStar files
base_path = "/epsen_data/data_overflow/sdss/apogee/dr17/apogee/spectro/redux/dr17/stars/"
aspcap_base_path = "/epsen_data/data_overflow/sdss/apogee/dr17/apogee/spectro/aspcap/dr17/synspec/"
# base_path = "C:/sdss_local/dr16/apogee/spectro/redux/r12/stars"  # path on my windows

# file name for master continuum normalized spectra, can be a path if file is not stored in current dir
contspac_file_name = "contspec_dr17_synspec.fits"

# file name for gaia all column, can be a path if file is not stored in current dir
gaia_allcolumns_f = "apogeedr17_syncspec_gaiadr3_xmatch_allcolumns.fits"
# file name for gaia selected columns but row matched to allstar, can be a path if file is not stored in current dir
gaia_rowmatch_f = "apogeedr17_syncspec_gaiadr3_xmatch.fits"

# file name for astroNN abundances, can be a path if file is not stored in current dir
astronn_chem_f = "astroNN_dr17_synspe_nn_chem.fits"

# file name for astroNN distances, can be a path if file is not stored in current dir
astronn_dist_f = "apogee_dr17_synspe_dr3_nn_dist.fits"

# file name for astroNN ages, can be a path if file is not stored in current dir
astronn_ages_f = "apogee_dr17_synspe_nn_ages.fits"

# neural network models folders name, can be a path if file is not stored in current dir
astronn_chem_model = "astroNN_0512_run002"
astronn_dist_model = "astroNN_gaia_dr17_model_3"
astronn_age_model = "APOKASC2_BCNN_age_combined_dr17_4"

# file name for astroNN APGOEE VAC, can be a path if file is not stored in current dir
astronn_apogee_vac_f = "apogee_astroNN-DR17_syncspec_dr3.fits"