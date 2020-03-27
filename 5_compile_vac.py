# This script uses all the fits files we have generated and compile VAC

from astropy.io import fits

from config import allstar_path, gaia_rowmatch_f, astronn_chem_f, astronn_dist_f, astronn_apogee_vac_f

allstar_data = fits.getdata(allstar_path)
f_gaia = fits.getdata(gaia_rowmatch_f)

f_chem= fits.getdata(astronn_chem_f)
f_dist = fits.getdata(astronn_dist_f)

c = [fits.Column(name='APOGEE_ID', array=allstar_data['APOGEE_ID'], format="18A"),
     fits.Column(name='LOCATION_ID', array=allstar_data['LOCATION_ID'], format="J"),
     fits.Column(name='TELESCOPE', array=allstar_data['TELESCOPE'], format="8A"),
     fits.Column(name='FIELD', array=allstar_data['FIELD'], format="19A"),
     fits.Column(name='RA_APOGEE', array=allstar_data['RA'], format='D'),
     fits.Column(name='DEC_APOGEE', array=allstar_data['DEC'], format='D'),
     fits.Column(name='VHELIO_AVG', array=allstar_data['VHELIO_AVG'], format='E'),
     fits.Column(name='source_id', array=f_gaia['source_id'], format='K'),
     fits.Column(name='TEFF', array=f_chem['astroNN'][:, 0], format='E'),
     fits.Column(name='TEFF_ERR', array=f_chem['astroNN_error'][:, 0], format='E'),
     fits.Column(name='LOGG', array=f_chem['astroNN'][:, 1], format='E'),
     fits.Column(name='LOGG_ERR', array=f_chem['astroNN_error'][:, 1], format='E'),
     fits.Column(name='C_H', array=f_chem['astroNN'][:, 2], format='E'),
     fits.Column(name='C_H_ERR', array=f_chem['astroNN_error'][:, 2], format='E'),
     fits.Column(name='CI_H', array=f_chem['astroNN'][:, 3], format='E'),
     fits.Column(name='CI_H_ERR', array=f_chem['astroNN_error'][:, 3], format='E'),
     fits.Column(name='N_H', array=f_chem['astroNN'][:, 4], format='E'),
     fits.Column(name='N_H_ERR', array=f_chem['astroNN_error'][:, 4], format='E'),
     fits.Column(name='O_H', array=f_chem['astroNN'][:, 5], format='E'),
     fits.Column(name='O_H_ERR', array=f_chem['astroNN_error'][:, 5], format='E'),
     fits.Column(name='NA_H', array=f_chem['astroNN'][:, 6], format='E'),
     fits.Column(name='NA_H_ERR', array=f_chem['astroNN_error'][:, 6], format='E'),
     fits.Column(name='MG_H', array=f_chem['astroNN'][:, 7], format='E'),
     fits.Column(name='MG_H_ERR', array=f_chem['astroNN_error'][:, 7], format='E'),
     fits.Column(name='AL_H', array=f_chem['astroNN'][:, 8], format='E'),
     fits.Column(name='AL_H_ERR', array=f_chem['astroNN_error'][:, 8], format='E'),
     fits.Column(name='SI_H', array=f_chem['astroNN'][:, 9], format='E'),
     fits.Column(name='SI_H_ERR', array=f_chem['astroNN_error'][:, 9], format='E'),
     fits.Column(name='P_H', array=f_chem['astroNN'][:, 10], format='E'),
     fits.Column(name='P_H_ERR', array=f_chem['astroNN_error'][:, 10], format='E'),
     fits.Column(name='S_H', array=f_chem['astroNN'][:, 11], format='E'),
     fits.Column(name='S_H_ERR', array=f_chem['astroNN_error'][:, 11], format='E'),
     fits.Column(name='K_H', array=f_chem['astroNN'][:, 12], format='E'),
     fits.Column(name='K_H_ERR', array=f_chem['astroNN_error'][:, 12], format='E'),
     fits.Column(name='CA_H', array=f_chem['astroNN'][:, 13], format='E'),
     fits.Column(name='CA_H_ERR', array=f_chem['astroNN_error'][:, 13], format='E'),
     fits.Column(name='TI_H', array=f_chem['astroNN'][:, 14], format='E'),
     fits.Column(name='TI_H_ERR', array=f_chem['astroNN_error'][:, 14], format='E'),
     fits.Column(name='TIII_H', array=f_chem['astroNN'][:, 15], format='E'),
     fits.Column(name='TIII_H_ERR', array=f_chem['astroNN_error'][:, 15], format='E'),
     fits.Column(name='V_H', array=f_chem['astroNN'][:, 16], format='E'),
     fits.Column(name='V_H_ERR', array=f_chem['astroNN_error'][:, 16], format='E'),
     fits.Column(name='CR_H', array=f_chem['astroNN'][:, 17], format='E'),
     fits.Column(name='CR_H_ERR', array=f_chem['astroNN_error'][:, 17], format='E'),
     fits.Column(name='MN_H', array=f_chem['astroNN'][:, 18], format='E'),
     fits.Column(name='MN_H_ERR', array=f_chem['astroNN_error'][:, 18], format='E'),
     fits.Column(name='FE_H', array=f_chem['astroNN'][:, 19], format='E'),
     fits.Column(name='FE_H_ERR', array=f_chem['astroNN_error'][:, 19], format='E'),
     fits.Column(name='CO_H', array=f_chem['astroNN'][:, 20], format='E'),
     fits.Column(name='CO_H_ERR', array=f_chem['astroNN_error'][:, 20], format='E'),
     fits.Column(name='NI_H', array=f_chem['astroNN'][:, 21], format='E'),
     fits.Column(name='NI_H_ERR', array=f_chem['astroNN_error'][:, 21], format='E'),
     fits.Column(name='dist', array=f_dist['dist'], format='E'),
     fits.Column(name='dist_error', array=f_dist['dist_error'], format='E'),
     fits.Column(name='dist_model_error', array=f_dist['dist_model_error'], format='E'),
     fits.Column(name='nn_parallax', array=f_dist['nn_parallax'], format='E'),
     fits.Column(name='nn_parallax_error', array=f_dist['nn_parallax_error'], format='E'),
     fits.Column(name='nn_parallax_model_error', array=f_dist['nn_parallax_model_error'], format='E'),
     fits.Column(name='fakemag', array=f_dist['fakemag'], format='E'),
     fits.Column(name='fakemag_error', array=f_dist['fakemag_error'], format='E'),
     fits.Column(name='RA', array=f_gaia['RA'], format='D'),
     fits.Column(name='DEC', array=f_gaia['DEC'], format='D'),
     fits.Column(name='pmra', array=f_gaia['pmra'], format='D'),
     fits.Column(name='pmra_error', array=f_gaia['pmra_error'], format='D'),
     fits.Column(name='pmdec', array=f_gaia['pmdec'], format='D'),
     fits.Column(name='pmdec_error', array=f_gaia['pmdec_error'], format='D'),
     fits.Column(name='phot_g_mean_mag', array=f_gaia['phot_g_mean_mag'], format='D'),
     fits.Column(name='bp_rp', array=f_gaia['bp_rp'], format='D')]

# save a fits
t = fits.BinTableHDU.from_columns(c)
t.writeto(astronn_apogee_vac_f)