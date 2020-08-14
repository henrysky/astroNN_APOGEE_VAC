# This script continuum normalizing all the APOGEE DR14 spectra and save them to a single fits

import os
import tqdm
import numpy as np
from astropy.io import fits
from astroNN.apogee.chips import apogee_continuum

from config import allstar14_path, base14_path, contspac14_file_name

if os.path.exists(contspac14_file_name):
    raise FileExistsError(f"{contspac14_file_name} already existed")

# read allstar
allstar_data = fits.getdata(allstar14_path)

# target bit in spectra to be set to zero
target_bit = [0, 1, 2, 3, 4, 5, 6, 7, 12]


def apstar_normalization(spectra, spectra_err, _spec_mask):
    return apogee_continuum(spectra=spectra, spectra_err=spectra_err, cont_mask=None, deg=2, dr=14, bitmask=_spec_mask,
                            target_bit=target_bit)

total_num = allstar_data['RA'].shape[0]
spec = np.zeros((total_num, 7514), dtype=np.float32)

for counter in tqdm.tqdm(range(0, total_num)):
    if allstar_data['LOCATION_ID'][counter] == 1:
        continue
    ap_path = f"{base14_path}/{allstar_data['TELESCOPE'][counter]}/{allstar_data['LOCATION_ID'][counter]}/{allstar_data['FILE'][counter]}"
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
hdu.writeto(contspac14_file_name)
