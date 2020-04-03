# This script uses neural network models to get abundances, distance, ages and save them to fits

import numpy as np
from astropy.io import fits
from astroNN.models import load_folder
from astroNN.gaia import extinction_correction, fakemag_to_pc, fakemag_to_parallax

from config import allstar_path, contspac_file_name, gaia_rowmatch_f, astronn_chem_model, astronn_dist_model, \
    astronn_age_model, astronn_chem_f, astronn_dist_f

allstar_data = fits.getdata(allstar_path)

allspec_f = fits.open(contspac_file_name)
all_spec = allspec_f[0].data
bad_spec_idx = np.all(all_spec == 0., axis=1)

# ====================================== Abundances ====================================== #

net = load_folder(astronn_chem_model)

pred, pred_error = net.test(all_spec)

# some spectra are all zeros, set prediction for those spectra to np.nan
pred[bad_spec_idx] = np.nan
pred_error['total'][bad_spec_idx] = np.nan

# deal with infinite error issue, set them to np.nan
inf_err_idx = np.array([pred_error['total'] == np.inf])[0]
pred[inf_err_idx] = np.nan
pred_error['total'][inf_err_idx] = np.nan

# save a fits
columns_list = [fits.Column(name='APOGEE_ID', array=allstar_data['APOGEE_ID'], format="18A"),
                fits.Column(name='LOCATION_ID', array=allstar_data['LOCATION_ID'], format="J"),
                fits.Column(name='RA', array=allstar_data['RA'], format='D'),
                fits.Column(name='DEC', array=allstar_data['DEC'], format='D'),
                fits.Column(name='astroNN', array=pred, format='22E'),
                fits.Column(name='astroNN_error', array=pred_error['total'], format='22E')]

t = fits.BinTableHDU.from_columns(columns_list)
t.writeto(astronn_chem_f)

# ====================================== Distance ====================================== #

corrected_K = extinction_correction(allstar_data['K'], allstar_data['AK_TARG'])

# cross matched APOGEE-Gaia DR2
apogeegaia_file = fits.getdata(gaia_rowmatch_f)
# add the offset we found for inv var weighting parallax
parallax = apogeegaia_file["parallax"] + 0.052
parallax_error = apogeegaia_file["parallax_error"]

# set negative parallax after constant offset correction to np.nan
parallax[(parallax < 0)] = np.nan
parallax_error[(parallax < 0)] = np.nan

# inference
net = load_folder(astronn_dist_model)
pred, pred_err = net.test(all_spec)
pred[:, 0][pred[:, 0] == 0.] = np.nan

# unit conversion
nn_dist, nn_dist_err = fakemag_to_pc(pred[:, 0], corrected_K, pred_err['total'][:, 0])
_, nn_dist_model_err = fakemag_to_pc(pred[:, 0], corrected_K, pred_err['model'][:, 0])
nn_parallax, nn_parallax_err = fakemag_to_parallax(pred[:, 0], corrected_K, pred_err['total'][:, 0])
_, nn_parallax_model_err = fakemag_to_parallax(pred[:, 0], corrected_K, pred_err['model'][:, 0])

# remove astropy units
nn_dist = nn_dist.value
nn_dist_err = nn_dist_err.value
nn_dist_model_err = nn_dist_model_err.value
nn_parallax = nn_parallax.value
nn_parallax_err = nn_parallax_err.value
nn_parallax_model_err = nn_parallax_model_err.value

# set bad value to np.nan
bad_idx = np.all(all_spec == 0., axis=1) | (pred[:, 0] < 0.)
nn_dist[bad_idx] = np.nan
nn_dist_err[bad_idx] = np.nan
nn_dist_model_err[bad_idx] = np.nan
nn_parallax[bad_idx] = np.nan
nn_parallax_err[bad_idx] = np.nan
nn_parallax_model_err[bad_idx] = np.nan

deno = ((1. / nn_parallax_model_err ** 2) + (1. / apogeegaia_file['parallax_error'] ** 2))
weighted_parallax = ((nn_parallax / nn_parallax_model_err ** 2) + (
            parallax / apogeegaia_file['parallax_error'] ** 2)) / deno
weighted_parallax_err = 1 / deno

# if one of them is -9999, use the value from the other one
idx = ((parallax == np.nan) & (nn_parallax != np.nan))
weighted_parallax[idx] = nn_parallax[idx]
weighted_parallax_err[idx] = nn_parallax_model_err[idx] ** 2  # still variance
# the other case
idx = ((nn_parallax == np.nan) & (parallax != np.nan))
weighted_parallax[idx] = parallax[idx]
weighted_parallax_err[idx] = parallax_error[idx] ** 2  # still variance

# if both of them is -9999, then np.nan
weighted_parallax[(parallax == np.nan) & (nn_parallax == np.nan)] = np.nan
weighted_parallax_err[(parallax == np.nan) & (nn_parallax == np.nan)] = np.nan

# need to take sqrt
weighted_parallax_err = np.ma.sqrt(np.ma.array(weighted_parallax_err, mask=(weighted_parallax_err == np.nan),
                                               fill_value=np.nan))

# change to weighted_dist
weighted_dist = 1000 / weighted_parallax
weighted_dist_err = weighted_dist * (weighted_parallax_err / weighted_parallax)

# if both of them is -9999, then np.nan
weighted_dist[(parallax == np.nan) & (nn_parallax == np.nan)] = np.nan
weighted_dist_err[(parallax == np.nan) & (nn_parallax == np.nan)] = np.nan

# prepare astropy fits columns
columns_list = [fits.Column(name='apogee_id', array=allstar_data['APOGEE_ID'], format="18A"),
                fits.Column(name='location_id', array=allstar_data['LOCATION_ID'], format="J"),
                fits.Column(name='ra_apogee', array=allstar_data['RA'], format='D'),
                fits.Column(name='dec_apogee', array=allstar_data['DEC'], format='D'),
                fits.Column(name='dist', array=nn_dist, format='D'),
                fits.Column(name='dist_error', array=nn_dist_err, format='D'),
                fits.Column(name='dist_model_error', array=nn_dist_model_err, format='D'),
                fits.Column(name='nn_parallax', array=nn_parallax, format='D'),
                fits.Column(name='nn_parallax_model_error', array=nn_parallax_model_err, format='D'),
                fits.Column(name='nn_parallax_error', array=nn_parallax_err, format='D'),
                fits.Column(name='fakemag', array=pred[:, 0], format='D'),
                fits.Column(name='fakemag_error', array=pred_err['total'][:, 0], format='D'),
                fits.Column(name='weighted_dist', array=weighted_dist, format='D'),
                fits.Column(name='weighted_dist_error', array=weighted_dist_err, format='D'),
                fits.Column(name='ra', array=apogeegaia_file['RA'], format='D'),
                fits.Column(name='dec', array=apogeegaia_file['DEC'], format='D'),
                fits.Column(name='pmra', array=apogeegaia_file['pmra'], format='D'),
                fits.Column(name='pmra_error', array=apogeegaia_file['pmra_error'], format='D'),
                fits.Column(name='pmdec', array=apogeegaia_file['pmdec'], format='D'),
                fits.Column(name='pmdec_error', array=apogeegaia_file['pmdec_error'], format='D'),
                fits.Column(name='phot_g_mean_mag', array=apogeegaia_file['phot_g_mean_mag'], format='D'),
                fits.Column(name='bp_rp', array=apogeegaia_file['bp_rp'], format='D')]

t = fits.BinTableHDU.from_columns(columns_list)
t.writeto(astronn_dist_f)

# ====================================== Ages ====================================== #

net = load_folder(astronn_age_model)

pred, pred_error = net.test(all_spec)