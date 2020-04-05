# This script continuum normalizing all the APOGEE spectra and save them to a single fits

import os
import time
import numpy as np
from astropy.io import fits
from astroNN.apogee.chips import apogee_continuum

from config import allstar_path, base_path, contspac_file_name

if os.path.exists(contspac_file_name):
    raise FileExistsError(f"{contspac_file_name} already existed")

# read allstar
allstar_data = fits.getdata(allstar_path)

# target bit in spectra to be set to zero
target_bit = [0, 1, 2, 3, 4, 5, 6, 7, 12]


def apstar_normalization(spectra, spectra_err, _spec_mask):
    return apogee_continuum(spectra=spectra, spectra_err=spectra_err, cont_mask=None, deg=2, dr=16, bitmask=_spec_mask,
                            target_bit=target_bit)

start_time = time.time()

total_num = allstar_data['RA'].shape[0]
spec = np.zeros((total_num, 7514), dtype=np.float32)

for counter in range(0, total_num):
    if allstar_data['LOCATION_ID'][counter] == 1:
        continue
    ap_path = f"{base_path}/{allstar_data['TELESCOPE'][counter]}/{allstar_data['FIELD'][counter]}/{allstar_data['FILE'][counter]}"
    if ap_path is False:
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

    if counter % 1000 == 0:
        print(f'Completed {counter} of {total_num}, {(time.time() - start_time):.2f}s elapsed')

# save a fits
hdu = fits.PrimaryHDU(spec)
hdu.writeto(contspac_file_name)
