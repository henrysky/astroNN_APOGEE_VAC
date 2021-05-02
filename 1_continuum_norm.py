# This script continuum normalizing all the APOGEE spectra and save them to a single fits

import os
import tqdm
import numpy as np
from astropy.io import fits
from astroNN.apogee.chips import apogee_continuum
from gaia_tools.xmatch import xmatch

from config import allstar_path, base_path, contspac_file_name, corr_flag

if os.path.exists(contspac_file_name):
    raise FileExistsError(f"{contspac_file_name} already existed")

# read allstar
allstar_data = fits.getdata(allstar_path)

# target bit in spectra to be set to zero
target_bit = [0, 1, 2, 3, 4, 5, 6, 7, 12]


def apstar_normalization(spectra, spectra_err, _spec_mask):
    return apogee_continuum(spectra=spectra, spectra_err=spectra_err, cont_mask=None, deg=2, dr=16, bitmask=_spec_mask,
                            target_bit=target_bit)

total_num = allstar_data['RA'].shape[0]
spec = np.zeros((total_num, 7514), dtype=np.float32)
good_flag = np.zeros(total_num, dtype=int)  # flag to indicate spectra exist (not all zeros)

for counter in tqdm.tqdm(range(0, total_num)):
    if allstar_data['LOCATION_ID'][counter] == 1:
        continue
    ap_path = f"{base_path}/{allstar_data['TELESCOPE'][counter]}/{allstar_data['FIELD'][counter]}/{allstar_data['FILE'][counter]}"
    if os.path.isfile(ap_path) is False:
        pass
    else:
        apstar_file = fits.open(ap_path)
        nvisits = apstar_file[0].header['NVISITS']
        _spec = apstar_file[1].data[0]
        _spec_err = apstar_file[2].data[0]
        _spec_mask = apstar_file[3].data[0]

        if not np.all(_spec == 0.) and not np.all(np.isnan(_spec)):
            _spec_err[np.isnan(_spec)] = 1e20
            _spec_mask[np.isnan(_spec)] = 1
            _spec[np.isnan(_spec)] = 1.
            _spec, _spec_err = apstar_normalization(_spec, _spec_err, _spec_mask)
            spec[counter] = _spec
            good_flag[counter] = 1

# save a fits
hdu = fits.PrimaryHDU(spec)
flag_hdu = fits.ImageHDU(good_flag)
hdul = fits.HDUList([hdu, flag_hdu])
hdul.writeto(contspac_file_name)

good_flag = good_flag.astype(bool)

if corr_flag:  # correct north-south
    north_spec_idx = np.where(allstar_data['TELESCOPE'][good_flag] == 'apo25m')[0]
    south_spec_idx = np.where(allstar_data['TELESCOPE'][good_flag] == 'lco25m')[0]

    idx1, idx2, _ = xmatch(allstar_data[good_flag][north_spec_idx], allstar_data[good_flag][south_spec_idx])

    spec = fits.getdata(contspac_file_name)

    lco_apo_mediff = np.median(spec[good_flag][south_spec_idx][idx2] - spec[good_flag][north_spec_idx][idx1], axis=0)
    np.save("lco_apo_median_dr17sync.npy", lco_apo_mediff)

    spec[good_flag][south_spec_idx] -= lco_apo_mediff

    os.rename(contspac_file_name, "{0}_{2}{1}".format(*os.path.splitext(contspac_file_name) + ("uncorrected",)))

    # save a fits
    hdu = fits.PrimaryHDU(spec)
    flag_hdu = fits.ImageHDU(good_flag.astype(int))
    hdul = fits.HDUList([hdu, flag_hdu])
    hdul.writeto(contspac_file_name)