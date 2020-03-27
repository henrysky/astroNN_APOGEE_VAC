# This script uses galpy to generate orbital parameters and save them to fits

import galpy
from astropy.io import fits

from config import allstar_path, gaia_rowmatch_f, astronn_dist_f

allstar_data = fits.getdata(allstar_path)
f_gaia = fits.getdata(gaia_rowmatch_f)

f_dist = fits.getdata(astronn_dist_f)
