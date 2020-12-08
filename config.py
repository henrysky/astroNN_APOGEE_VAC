# This script defines global configuration/constant

from astroNN.apogee import allstar

# milkyway potential parameters
_R0 = 8.125  # distance to galactic center in kpc
_v0 = 220.  # sun tangential velocity in km/x
_z0 = 0.0208  # sun above galactic plane in kpc

# gaia table to use to cross matching
gaia_table_name = "gaiaedr3.gaia_source"

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
allstar_path = "allStar-r13-l33-58932beta.fits"  # path on my windows

# base path to apStar/asStar files
base_path = "/yngve_data/sdss/apogee/dr16/apogee/spectro/redux/r12/stars"
# base_path = "C:/sdss_local/dr16/apogee/spectro/redux/r12/stars"  # path on my windows

# file name for master continuum normalized spectra, can be a path if file is not stored in current dir
contspac_file_name = "contspec_dr16.fits"

# file name for gaia all column, can be a path if file is not stored in current dir
gaia_allcolumns_f = "apogeedr17_gaiaedr3_xmatch_allcolumns.fits"
# file name for gaia selected columns but row matched to allstar, can be a path if file is not stored in current dir
gaia_rowmatch_f = "apogeedr17_gaiaedr3_xmatch.fits"

# file name for astroNN abundances, can be a path if file is not stored in current dir
astronn_chem_f = "astroNN_dr17_nn_chem.fits"

# file name for astroNN distances, can be a path if file is not stored in current dir
astronn_dist_f = "apogee_dr17_edr3_nn_dist.fits"

# file name for astroNN ages, can be a path if file is not stored in current dir
astronn_ages_f = "apogee_dr17_nn_ages.fits"

# file name for astroNN ages, can be a path if file is not stored in current dir
galpy_orbitparams_f = "apogee_dr17_orbitparams_edr3.fits"

# neural network models folders name, can be a path if file is not stored in current dir
astronn_chem_model = "astroNN_0617_run001"  # https://github.com/henrysky/astroNN_spectra_paper_figures/
astronn_dist_model = "astroNN_constant_model_reduced"  # https://github.com/henrysky/astroNN_gaia_dr2_paper/
astronn_age_model = "APOKASC2_BCNN_age_only_corrections0.1"

# file name for astroNN APGOEE VAC, can be a path if file is not stored in current dir
astronn_apogee_vac_f = "apogee_astroNN-DR17_edr3.fits"
