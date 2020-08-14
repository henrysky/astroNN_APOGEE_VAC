# This script continuum normalizing all the APOGEE spectra and save them to a single fits

import os
import tqdm
import numpy as np
from astropy.io import fits
from astroNN.apogee.chips import apogee_continuum
from gaia_tools.xmatch import xmatch

from config import allstar_path, allstar14_path, base_path, contspac_file_name, corr_flag, contspac14_file_name

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

for counter in tqdm.tqdm(range(0, total_num)):
    if allstar_data['LOCATION_ID'][counter] == 1:
        continue
    ap_path = f"{base_path}/{allstar_data['TELESCOPE'][counter]}/{allstar_data['FIELD'][counter]}/{allstar_data['FILE'][counter]}"
    if os.path.exists(ap_path) is False:
        pass
    else:
        apstar_file = fits.open(ap_path)
        nvisits = apstar_file[0].header['NVISITS']
        if nvisits == 1:
            _spec = apstar_file[1].data
            _spec_err = apstar_file[2].data
            _spec_mask = apstar_file[3].data
        else:
            _spec = apstar_file[1].data[1]
            _spec_err = apstar_file[2].data[1]
            _spec_mask = apstar_file[3].data[1]

        if not np.all(_spec == 0.):
            _spec, _spec_err = apstar_normalization(_spec, _spec_err, _spec_mask)
            spec[counter] = _spec

# save a fits
hdu = fits.PrimaryHDU(spec)
hdu.writeto(contspac_file_name)

if corr_flag:
    allstar_dr14_data = fits.getdata(allstar14_path)

    north_spec = np.where(allstar_data['TELESCOPE'] == 'apo25m')[0]
    south_spec = np.where(allstar_data['TELESCOPE'] == 'lco25m')[0]

    idx1, idx2, _ = xmatch(allstar_data[north_spec], allstar_dr14_data)
    idx3, idx4, _ = xmatch(allstar_data[north_spec], allstar_data[south_spec])

    dr14_spec = fits.getdata(contspac14_file_name)

    apo16_apo14_mediff = np.median(spec[north_spec][idx1] - dr14_spec[idx2], axis=0)
    lco16_apo14_mediff = np.median(spec[south_spec][idx4] - spec[north_spec][idx3], axis=0) + apo16_apo14_mediff
    np.save("apo16_apo14_median_dr17.npy", apo16_apo14_mediff)
    np.save("lco16_apo14_median_dr17.npy", lco16_apo14_mediff)

    for counter in tqdm.tqdm(range(0, total_num)):
        if not np.all(spec[counter] == 0.):
            telescope = allstar_data['TELESCOPE'][counter]

            if telescope == 'lco25m':
                corr = np.array(lco16_apo14_mediff)
            elif telescope == 'apo25m':
                corr = np.array(apo16_apo14_mediff)
            else:
                corr = np.zeros_like(spec[0])

            masked = (spec[counter] == 1.)

            spec[counter] = spec[counter] - corr
            spec[counter][masked] = 1.

    os.rename(contspac_file_name, "{0}_{2}{1}".format(*os.path.splitext(contspac_file_name) + ("uncorrected",)))

    # save a fits
    hdu = fits.PrimaryHDU(spec)
    hdu.writeto(contspac_file_name)
