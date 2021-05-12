# This script uses neural network models to get abundances, distance, ages and save them to fits

import os
import h5py
import numpy as np
from astropy.io import fits
from astroNN.models import load_folder
from astroNN.datasets import xmatch
from astroNN.gaia import extinction_correction, fakemag_to_pc, fakemag_to_parallax
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.interpolate import interp1d
import scipy.optimize as op
from astroNN.apogee import bitmask_boolean

from config import (
    allstar_path,
    contspac_file_name,
    gaia_rowmatch_f,
    astronn_chem_model,
    astronn_dist_model,
    astronn_age_model,
    astronn_chem_f,
    astronn_dist_f,
    astronn_ages_f,
)

if os.path.exists(astronn_chem_f):
    raise FileExistsError(f"{astronn_chem_f} already existed")
if os.path.exists(astronn_dist_f):
    raise FileExistsError(f"{astronn_dist_f} already existed")
if os.path.exists(astronn_ages_f):
    raise FileExistsError(f"{astronn_ages_f} already existed")

allstar_data = fits.getdata(allstar_path)

allspec_f = fits.open(contspac_file_name)
all_spec = allspec_f[0].data
bad_spec_idx = (
    ~np.array(allspec_f[1].data, bool)
    | bitmask_boolean(allstar_data["STARFLAG"], target_bit=[0, 3, 4])
)[0]

# ====================================== Abundances ====================================== #

net = load_folder(astronn_chem_model)

pred, pred_error = net.predict_dataset(fits.getdata(contspac_file_name))

# some spectra are all zeros, set prediction for those spectra to np.nan
pred[bad_spec_idx] = np.nan
pred_error["total"][bad_spec_idx] = np.nan

# deal with infinite error issue, set them to np.nan
inf_err_idx = np.array([pred_error["total"] == np.inf])[0]
pred[inf_err_idx] = np.nan
pred[bad_spec_idx] = np.nan
pred_error["total"][inf_err_idx] = np.nan
pred_error["total"][bad_spec_idx] = np.nan

# save a fits
columns_list = [
    fits.Column(name="APOGEE_ID", array=allstar_data["APOGEE_ID"], format="18A"),
    fits.Column(name="LOCATION_ID", array=allstar_data["LOCATION_ID"], format="J"),
    fits.Column(name="RA", array=allstar_data["RA"], format="D"),
    fits.Column(name="DEC", array=allstar_data["DEC"], format="D"),
    fits.Column(name="astroNN", array=pred, format="22E"),
    fits.Column(name="astroNN_error", array=pred_error["total"], format="22E"),
]

t = fits.BinTableHDU.from_columns(columns_list)
t.writeto(astronn_chem_f)

# ====================================== Distance ====================================== #

corrected_K = extinction_correction(allstar_data["K"], allstar_data["AK_TARG"])

# cross matched APOGEE-Gaia
apogeegaia_file = fits.getdata(gaia_rowmatch_f)
# add the offset we found for inv var weighting parallax
try:
    parallax = apogeegaia_file["parallax_w_zp"]
except KeyError:
    parallax = apogeegaia_file["parallax"]

parallax_error = apogeegaia_file["parallax_error"]

# set negative parallax after constant offset correction to np.nan
parallax_error[(parallax < 0) & (parallax > 1e10)] = np.nan
parallax[(parallax < 0) & (parallax > 1e10)] = np.nan

# inference
net = load_folder(astronn_dist_model)
# pred, pred_err = net.predict(all_spec)
pred, pred_err = net.predict_dataset(fits.getdata(contspac_file_name))

# unit conversion
nn_dist, nn_dist_err = fakemag_to_pc(pred[:, 0], corrected_K, pred_err["total"][:, 0])
_, nn_dist_model_err = fakemag_to_pc(pred[:, 0], corrected_K, pred_err["model"][:, 0])
nn_parallax, nn_parallax_err = fakemag_to_parallax(
    pred[:, 0], corrected_K, pred_err["total"][:, 0]
)
_, nn_parallax_model_err = fakemag_to_parallax(
    pred[:, 0], corrected_K, pred_err["model"][:, 0]
)

# remove astropy units
nn_dist = nn_dist.value
nn_dist_err = nn_dist_err.value
nn_dist_model_err = nn_dist_model_err.value
nn_parallax = nn_parallax.value
nn_parallax_err = nn_parallax_err.value
nn_parallax_model_err = nn_parallax_model_err.value

# set bad value to np.nan
bad_idx = (
    bad_spec_idx
    | (pred[:, 0] < 1e-4)
    | (nn_dist == -9999.0)
    | np.isnan(nn_dist)
    | np.isinf(nn_dist)
)
nn_dist[bad_idx] = np.nan
nn_dist_err[bad_idx] = np.nan
nn_dist_model_err[bad_idx] = np.nan
nn_parallax[bad_idx] = np.nan
nn_parallax_err[bad_idx] = np.nan
nn_parallax_model_err[bad_idx] = np.nan
pred[:, 0][bad_idx] = np.nan
pred_err["total"][bad_idx] = np.nan
pred_err["model"][bad_idx] = np.nan

deno = (1.0 / nn_parallax_model_err ** 2) + (
    1.0 / apogeegaia_file["parallax_error"] ** 2
)
weighted_parallax = (
    (nn_parallax / nn_parallax_model_err ** 2)
    + (parallax / apogeegaia_file["parallax_error"] ** 2)
) / deno
weighted_parallax_err = 1 / deno

# if one of them is -9999 or NaN, use the value from the other one
idx = (parallax == np.nan) & (nn_parallax != np.nan)
weighted_parallax[idx] = nn_parallax[idx]
weighted_parallax_err[idx] = nn_parallax_model_err[idx] ** 2  # still variance
# the other case
idx = (nn_parallax == np.nan) & (parallax != np.nan)
weighted_parallax[idx] = parallax[idx]
weighted_parallax_err[idx] = parallax_error[idx] ** 2  # still variance

# if both of them is -9999 or NaN, then np.nan
weighted_parallax[(parallax == np.nan) & (nn_parallax == np.nan)] = np.nan
weighted_parallax_err[(parallax == np.nan) & (nn_parallax == np.nan)] = np.nan

# need to take sqrt
weighted_parallax_err = np.ma.sqrt(
    np.ma.array(
        weighted_parallax_err, mask=(weighted_parallax_err == np.nan), fill_value=np.nan
    )
)

# change to weighted_dist
weighted_dist = 1000 / weighted_parallax
weighted_dist_err = weighted_dist * (weighted_parallax_err / weighted_parallax)

# if both of them is -9999 or NaN, then np.nan
weighted_dist[(parallax == np.nan) & (nn_parallax == np.nan)] = np.nan
weighted_dist_err[(parallax == np.nan) & (nn_parallax == np.nan)] = np.nan

# prepare astropy fits columns
columns_list = [
    fits.Column(name="apogee_id", array=allstar_data["APOGEE_ID"], format="18A"),
    fits.Column(name="location_id", array=allstar_data["LOCATION_ID"], format="J"),
    fits.Column(name="ra_apogee", array=allstar_data["RA"], format="D"),
    fits.Column(name="dec_apogee", array=allstar_data["DEC"], format="D"),
    fits.Column(name="dist", array=nn_dist, format="D"),
    fits.Column(name="dist_error", array=nn_dist_err, format="D"),
    fits.Column(name="dist_model_error", array=nn_dist_model_err, format="D"),
    fits.Column(name="nn_parallax", array=nn_parallax, format="D"),
    fits.Column(
        name="nn_parallax_model_error", array=nn_parallax_model_err, format="D"
    ),
    fits.Column(name="nn_parallax_error", array=nn_parallax_err, format="D"),
    fits.Column(name="fakemag", array=pred[:, 0], format="D"),
    fits.Column(name="fakemag_error", array=pred_err["total"][:, 0], format="D"),
    fits.Column(name="weighted_dist", array=weighted_dist, format="D"),
    fits.Column(name="weighted_dist_error", array=weighted_dist_err, format="D"),
    fits.Column(name="ra", array=apogeegaia_file["ra"], format="D"),
    fits.Column(name="dec", array=apogeegaia_file["dec"], format="D"),
    fits.Column(name="ra_error", array=apogeegaia_file["ra_error"], format="D"),
    fits.Column(name="dec_error", array=apogeegaia_file["dec_error"], format="D"),
    fits.Column(name="pmra", array=apogeegaia_file["pmra"], format="D"),
    fits.Column(name="pmra_error", array=apogeegaia_file["pmra_error"], format="D"),
    fits.Column(name="pmdec", array=apogeegaia_file["pmdec"], format="D"),
    fits.Column(name="pmdec_error", array=apogeegaia_file["pmdec_error"], format="D"),
    fits.Column(
        name="phot_g_mean_mag", array=apogeegaia_file["phot_g_mean_mag"], format="D"
    ),
    fits.Column(name="bp_rp", array=apogeegaia_file["bp_rp"], format="D"),
]

t = fits.BinTableHDU.from_columns(columns_list)
t.writeto(astronn_dist_f)

# ====================================== Ages ====================================== #

# APOKASC processing
apokasc3 = fits.getdata("APOKASC_cat_v6.6.1.fits.zip")
good_ages = apokasc3["APOKASC2_AGE"] != -9999.0
apokasc3 = apokasc3[good_ages]

# ages for low metallicity
f_age_low_M = fits.getdata("kepler_low_metallicity_with_samples.fits")

allstar_f = fits.getdata(allstar_path)
ra = allstar_f["ra"]
dec = allstar_f["dec"]
ra[0] = 0
dec[0] = 0
idx_1, idx_2, sep = xmatch(apokasc3["RA"], apokasc3["DEC"], ra, dec)
idx_3, idx_4, sep = xmatch(f_age_low_M["RA"], f_age_low_M["DEC"], ra, dec)

idx_combined, unique_indices = np.unique(
    np.concatenate([idx_4, idx_2]), return_index=True
)

all_age = np.concatenate([f_age_low_M["Age_med"] / 1e9, apokasc3["APOKASC2_AGE"]])[
    unique_indices
]
all_age_err = np.concatenate(
    [f_age_low_M["Age_Sd"] / 1e9, apokasc3["APOKASC2_AGE_MERR"]]
)[unique_indices]
all_mass = np.concatenate([f_age_low_M["Mass_med"], apokasc3["APOKASC2_MASS"]])[
    unique_indices
]
all_mass_err = np.concatenate(
    [f_age_low_M["Mass_Sd"], apokasc3["APOKASC2_MASS_RANERR"]]
)[unique_indices]

net = load_folder(astronn_age_model)

pred, pred_error = net.predict_dataset(fits.getdata(contspac_file_name))

# some spectra are all zeros, set prediction for those spectra to NaN
pred[bad_spec_idx] = np.nan
pred_error["total"][bad_spec_idx] = np.nan
pred_error["model"][bad_spec_idx] = np.nan

# If age too small or too large set them to Nan
bad_age = (pred < 0.1) | (pred > 20.0)
pred[bad_age] = np.nan
pred_error["total"][bad_age] = np.nan
pred_error["model"][bad_age] = np.nan

# deal with infinite error issue if it exists, set them to NaN
inf_err_idx = pred_error["total"] == np.inf
pred[inf_err_idx] = np.nan
pred_error["total"][inf_err_idx] = np.nan
pred_error["model"][inf_err_idx] = np.nan

out = lowess(all_age, pred[:, 0][idx_combined], frac=0.8, delta=0.1)
correction = interp1d(
    out[:, 0], out[:, 1], bounds_error=False, fill_value="extrapolate"
)

columns_list = [
    fits.Column(name="apogee_id", array=allstar_data["APOGEE_ID"], format="18A"),
    fits.Column(name="location_id", array=allstar_data["LOCATION_ID"], format="J"),
    fits.Column(name="ra_apogee", array=allstar_data["RA"], format="D"),
    fits.Column(name="dec_apogee", array=allstar_data["DEC"], format="D"),
    fits.Column(name="age", array=pred[:, 0], format="D"),
    fits.Column(
        name="age_linear_correct",
        array=(pred[:, 0] - 0.23834204) / 0.88723217,
        format="D",
    ),
    fits.Column(name="age_lowess_correct", array=correction(pred[:, 0]), format="D"),
    fits.Column(name="age_total_error", array=pred_error["total"][:, 0], format="D"),
    fits.Column(name="age_model_error", array=pred_error["model"][:, 0], format="D"),
    fits.Column(name="mass", array=pred[:, 1], format="D"),
    fits.Column(name="mass_total_error", array=pred_error["total"][:, 1], format="D"),
    fits.Column(name="mass_model_error", array=pred_error["model"][:, 1], format="D"),
]

t = fits.BinTableHDU.from_columns(columns_list)
t.writeto(astronn_ages_f)
