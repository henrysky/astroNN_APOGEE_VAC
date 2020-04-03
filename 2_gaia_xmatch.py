# This script cross matching APOGEE-Gaia to get all columns and save them to a single fits

import os
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astroquery.gaia import Gaia
from config import allstar_path, gaia_allcolumns_f, gaia_rowmatch_f, gaia_table_name

# open apogee allstar
allstar_data = fits.getdata(allstar_path)
ra_apogee = allstar_data['ra']
gaia_matched_idx = np.where(allstar_data['GAIA_SOURCE_ID'] != 0)

# try to login if any
if os.stat("gaia_credential").st_size != 0:
       Gaia.login(credentials_file='gaia_credential')

t = Table({'source_id':allstar_data['GAIA_SOURCE_ID']})
t.write('temptable.xml', format='votable')

# launching job at https://gea.esac.esa.int/archive/
job = Gaia.launch_job_async(
    f"""
    select g.*
    from {gaia_table_name} as g 
    inner join tap_upload.my_table as m on m.source_id = g.source_id
    """,
    upload_resource='temptable.xml',
    upload_table_name="my_table")

os.remove('temptable.xml')


# parse job result and save
gaia2_matches = job.results
gaia2_matches.remove_columns(['datalink_url', 'epoch_photometry_url', 'designation', 'phot_variable_flag'])
gaia2_matches.write(gaia_allcolumns_f)

# prepare the row matching allstar-gaia
xmatched_allcolumns = fits.getdata(gaia_allcolumns_f)

ra = np.ones(ra_apogee.shape[0]) * np.nan
dec = np.ones(ra_apogee.shape[0]) * np.nan
ref_epoch = np.ones(ra_apogee.shape[0]) * np.nan
parallax = np.ones(ra_apogee.shape[0]) * np.nan
parallax_error = np.ones(ra_apogee.shape[0]) * np.nan
visibility_periods_used  = np.ones(ra_apogee.shape[0]) * np.nan
astrometric_chi2_al  = np.ones(ra_apogee.shape[0]) * np.nan
astrometric_n_good_obs_al  = np.ones(ra_apogee.shape[0]) * np.nan
pmra = np.ones(ra_apogee.shape[0]) * np.nan
pmra_error = np.ones(ra_apogee.shape[0]) * np.nan
pmdec = np.ones(ra_apogee.shape[0]) * np.nan
pmdec_error = np.ones(ra_apogee.shape[0]) * np.nan
phot_g_mean_mag = np.ones(ra_apogee.shape[0]) * np.nan
bp_rp = np.ones(ra_apogee.shape[0]) * np.nan
source_id = np.zeros(ra_apogee.shape[0], dtype=np.int64) - 1

ra[gaia_matched_idx] = xmatched_allcolumns['ra']
dec[gaia_matched_idx] = xmatched_allcolumns['dec']
ref_epoch[gaia_matched_idx] = xmatched_allcolumns['ref_epoch']
parallax[gaia_matched_idx] = xmatched_allcolumns['parallax']
parallax_error[gaia_matched_idx] = xmatched_allcolumns['parallax_error']
visibility_periods_used[gaia_matched_idx] = xmatched_allcolumns['visibility_periods_used']
astrometric_chi2_al[gaia_matched_idx] = xmatched_allcolumns['astrometric_chi2_al']
astrometric_n_good_obs_al[gaia_matched_idx] = xmatched_allcolumns['astrometric_n_good_obs_al']
pmra[gaia_matched_idx] = xmatched_allcolumns['pmra']
pmra_error[gaia_matched_idx] = xmatched_allcolumns['pmra_error']
pmdec[gaia_matched_idx] = xmatched_allcolumns['pmdec']
pmdec_error[gaia_matched_idx] = xmatched_allcolumns['pmdec_error']
phot_g_mean_mag[gaia_matched_idx] = xmatched_allcolumns['phot_g_mean_mag']
bp_rp[gaia_matched_idx] = xmatched_allcolumns['bp_rp']
source_id[gaia_matched_idx] = xmatched_allcolumns['source_id']

col = [fits.Column(name='APOGEE_ID', array=allstar_data['APOGEE_ID'], format="18A"),
       fits.Column(name='LOCATION_ID', array=allstar_data['LOCATION_ID'], format="J"),
       fits.Column(name='RA_APOGEE', array=allstar_data['RA'], format='D'),
       fits.Column(name='DEC_APOGEE', array=allstar_data['DEC'], format='D'),
       fits.Column(name='RA', array=ra, format='D'),
       fits.Column(name='DEC', array=dec, format='D'),
       fits.Column(name='ref_epoch', array=ref_epoch, format='D'),
       fits.Column(name='parallax', array=parallax, format='D'),
       fits.Column(name='parallax_error', array=parallax_error, format='D'),
       fits.Column(name='visibility_periods_used', array=visibility_periods_used, format='I'),
       fits.Column(name='astrometric_chi2_al', array=astrometric_chi2_al, format='E'),
       fits.Column(name='astrometric_n_good_obs_al', array=astrometric_n_good_obs_al, format='J'),
       fits.Column(name='pmra', array=pmra, format='D'),
       fits.Column(name='pmra_error', array=pmra_error, format='D'),
       fits.Column(name='pmdec', array=pmdec, format='D'),
       fits.Column(name='pmdec_error', array=pmdec_error, format='D'),
       fits.Column(name='phot_g_mean_mag', array=phot_g_mean_mag, format='D'),
       fits.Column(name='bp_rp', array=bp_rp, format='D'),
       fits.Column(name='source_id', array=source_id, format='K')]

t = fits.BinTableHDU.from_columns(col)
t.writeto(gaia_rowmatch_f)