# This script gets all the APOGEE spectra normalized by ASPCAP and save them to a single fits

import os
import tqdm
import numpy as np
from astropy.io import fits
from astroNN.apogee.chips import gap_delete

from config import allstar_path, aspcap_base_path

# read allstar
allstar_data = fits.getdata(allstar_path)

# target bit in spectra to be set to zero
target_bit = [0, 1, 2, 3, 4, 5, 6, 7, 12]

total_num = allstar_data['RA'].shape[0]
spec = np.zeros((total_num, 7514), dtype=np.float32)
good_flag = np.zeros(total_num, dtype=int)  # flag to indicate spectra exist (not all zeros)

for counter in tqdm.tqdm(range(0, total_num)):
    if allstar_data['LOCATION_ID'][counter] == 1:
        continue
    ap_path = f"{aspcap_base_path}/{allstar_data['TELESCOPE'][counter]}/{allstar_data['FIELD'][counter]}/aspcapStar-dr17-{allstar_data['APOGEE_ID'][counter]}.fits"
    if os.path.isfile(ap_path) is False:
        pass
    else:
        apstar_file = fits.open(ap_path)
        _spec = gap_delete(apstar_file[1].data)

        if not np.all(_spec == 0.) and not np.all(np.isnan(_spec)):
            _spec[np.isnan(_spec)] = 1.
            spec[counter] = _spec
            good_flag[counter] = 1

# save a fits
hdu = fits.PrimaryHDU(spec)
flag_hdu = fits.ImageHDU(good_flag)
hdul = fits.HDUList([hdu, flag_hdu])
hdul.writeto("aspcap_dr17_synspec.fits")
