# This script cross matching APOGEE-Gaia to get all columns and save them to a single fits

import os
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astroquery.gaia import Gaia
from config import allstar_path, gaia_allcolumns_f, gaia_rowmatch_f, gaia_table_name
from utils import dr2source_dr3source

try:
       from zero_point import zpt
       zpt.load_tables()
       ZPT = True
except ImportError:
       ZPT = False

# open apogee allstar
allstar_data = fits.getdata(allstar_path)
ra_apogee = allstar_data['ra']
# dr3_source_id = dr2source_dr3source(allstar_data['GAIAEDR3_SOURCE_ID'])
dr3_source_id = allstar_data['GAIAEDR3_SOURCE_ID']

gaia_matched_idx = np.where(dr3_source_id > 1)

# try to login if any
if os.stat("gaia_credential").st_size != 0:
       Gaia.login(credentials_file='gaia_credential')

t = Table({'source_id': dr3_source_id})

# launching job at https://gea.esac.esa.int/archive/
job = Gaia.launch_job_async(
    f"""
    select g.*
    from {gaia_table_name} as g 
    inner join tap_upload.my_table as m on m.source_id = g.source_id
    """,
    upload_resource=t,
    upload_table_name="my_table")

# parse job result and save
gaia2_matches = job.results
gaia2_matches.remove_columns(['designation'])
gaia2_matches.write(gaia_allcolumns_f)

# prepare the row matching allstar-gaia
xmatched_allcolumns = fits.getdata(gaia_allcolumns_f)

ra = np.ones(ra_apogee.shape[0]) * np.nan
dec = np.ones(ra_apogee.shape[0]) * np.nan
ra_error = np.ones(ra_apogee.shape[0]) * np.nan
dec_error = np.ones(ra_apogee.shape[0]) * np.nan
ra_dec_corr = np.ones(ra_apogee.shape[0]) * np.nan
ra_parallax_corr = np.ones(ra_apogee.shape[0]) * np.nan
ra_pmra_corr = np.ones(ra_apogee.shape[0]) * np.nan
ra_pmdec_corr = np.ones(ra_apogee.shape[0]) * np.nan
dec_parallax_corr = np.ones(ra_apogee.shape[0]) * np.nan
dec_pmra_corr = np.ones(ra_apogee.shape[0]) * np.nan
dec_pmdec_corr = np.ones(ra_apogee.shape[0]) * np.nan
parallax_pmra_corr = np.ones(ra_apogee.shape[0]) * np.nan
parallax_pmdec_corr = np.ones(ra_apogee.shape[0]) * np.nan
pmra_pmdec_corr = np.ones(ra_apogee.shape[0]) * np.nan
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
# dr2_source_id = np.zeros(ra_apogee.shape[0], dtype=np.int64) - 1
# new in EDR3
bp_g = np.ones(ra_apogee.shape[0]) * np.nan
g_rp = np.ones(ra_apogee.shape[0]) * np.nan
pseudocolour = np.ones(ra_apogee.shape[0]) * np.nan
pseudocolour_error = np.ones(ra_apogee.shape[0]) * np.nan
nu_eff_used_in_astrometry = np.ones(ra_apogee.shape[0]) * np.nan
astrometric_params_solved = np.ones(ra_apogee.shape[0]) * np.nan
ecl_lat = np.ones(ra_apogee.shape[0]) * np.nan
ruwe = np.ones(ra_apogee.shape[0]) * np.nan
ipd_gof_harmonic_amplitude = np.ones(ra_apogee.shape[0]) * np.nan
ipd_frac_multi_peak = np.ones(ra_apogee.shape[0]) * np.nan
phot_bp_rp_excess_factor = np.ones(ra_apogee.shape[0]) * np.nan

ra[gaia_matched_idx] = xmatched_allcolumns['ra']
dec[gaia_matched_idx] = xmatched_allcolumns['dec']
ra_error[gaia_matched_idx] = xmatched_allcolumns['ra_error']
dec_error[gaia_matched_idx] = xmatched_allcolumns['dec_error']
ra_dec_corr[gaia_matched_idx] = xmatched_allcolumns['ra_dec_corr']
ra_parallax_corr[gaia_matched_idx] = xmatched_allcolumns['ra_parallax_corr']
ra_pmra_corr[gaia_matched_idx] = xmatched_allcolumns['ra_pmra_corr']
ra_pmdec_corr[gaia_matched_idx] = xmatched_allcolumns['ra_pmdec_corr']
dec_parallax_corr[gaia_matched_idx] = xmatched_allcolumns['dec_parallax_corr']
dec_pmra_corr[gaia_matched_idx] = xmatched_allcolumns['dec_pmra_corr']
dec_pmdec_corr[gaia_matched_idx] = xmatched_allcolumns['dec_pmdec_corr']
parallax_pmra_corr[gaia_matched_idx] = xmatched_allcolumns['parallax_pmra_corr']
parallax_pmdec_corr[gaia_matched_idx] = xmatched_allcolumns['parallax_pmdec_corr']
pmra_pmdec_corr[gaia_matched_idx] = xmatched_allcolumns['pmra_pmdec_corr']
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
# dr2_source_id[allstar_data["GAIA_SOURCE_ID"] > 1] = \
#        allstar_data["GAIA_SOURCE_ID"][allstar_data["GAIA_SOURCE_ID"] > 1]
# new in EDR3
bp_g[gaia_matched_idx] = xmatched_allcolumns['bp_g']
g_rp[gaia_matched_idx] = xmatched_allcolumns['g_rp']
pseudocolour[gaia_matched_idx] = xmatched_allcolumns['pseudocolour']
pseudocolour_error[gaia_matched_idx] = xmatched_allcolumns['pseudocolour_error']
nu_eff_used_in_astrometry[gaia_matched_idx] = xmatched_allcolumns['nu_eff_used_in_astrometry']
astrometric_params_solved[gaia_matched_idx] = xmatched_allcolumns['astrometric_params_solved']
ecl_lat[gaia_matched_idx] = xmatched_allcolumns['ecl_lat']
ruwe[gaia_matched_idx] = xmatched_allcolumns['ruwe']
ipd_gof_harmonic_amplitude[gaia_matched_idx] = xmatched_allcolumns['ipd_gof_harmonic_amplitude']
ipd_frac_multi_peak[gaia_matched_idx] = xmatched_allcolumns['ipd_frac_multi_peak']
phot_bp_rp_excess_factor[gaia_matched_idx] = xmatched_allcolumns['ipd_frac_multi_peak']

if ZPT:
       good_idx = np.where(((astrometric_params_solved==31) | (astrometric_params_solved==95)) & 
                           (phot_g_mean_mag<21) & (6<phot_g_mean_mag))[0]
       zp = zpt.get_zpt(phot_g_mean_mag[good_idx], 
                        nu_eff_used_in_astrometry[good_idx], 
                        pseudocolour[good_idx], 
                        ecl_lat[good_idx], 
                        astrometric_params_solved[good_idx])
       # use median zero-point for all stars by default
       zp_row_matched = np.ones(len(ra)) * np.median(zp)

       zp_row_matched[good_idx] = zp


col = [fits.Column(name='APOGEE_ID', array=allstar_data['APOGEE_ID'], format="18A"),
       fits.Column(name='LOCATION_ID', array=allstar_data['LOCATION_ID'], format="J"),
       fits.Column(name='RA_APOGEE', array=allstar_data['RA'], format='D'),
       fits.Column(name='DEC_APOGEE', array=allstar_data['DEC'], format='D'),
       fits.Column(name='ra', array=ra, format='D'),
       fits.Column(name='dec', array=dec, format='D'),
       fits.Column(name='ra_error', array=ra_error, format='D'),
       fits.Column(name='dec_error', array=dec_error, format='D'),
       fits.Column(name='ra_dec_corr', array=ra_dec_corr, format='D'),
       fits.Column(name='ra_parallax_corr', array=ra_parallax_corr, format='D'),
       fits.Column(name='ra_pmra_corr', array=ra_pmra_corr, format='D'),
       fits.Column(name='ra_pmdec_corr', array=ra_pmdec_corr, format='D'),
       fits.Column(name='dec_parallax_corr', array=dec_parallax_corr, format='D'),
       fits.Column(name='dec_pmra_corr', array=dec_pmra_corr, format='D'),
       fits.Column(name='dec_pmdec_corr', array=dec_pmdec_corr, format='D'),
       fits.Column(name='parallax_pmra_corr', array=parallax_pmra_corr, format='D'),
       fits.Column(name='parallax_pmdec_corr', array=parallax_pmdec_corr, format='D'),
       fits.Column(name='pmra_pmdec_corr', array=pmra_pmdec_corr, format='D'),
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
       fits.Column(name='bp_g', array=bp_g, format='D'),
       fits.Column(name='g_rp', array=g_rp, format='D'),
       fits.Column(name='source_id', array=source_id, format='K'),
       # fits.Column(name='dr2_source_id', array=allstar_data["GAIA_SOURCE_ID"], format='K'),
       fits.Column(name='pseudocolour', array=pseudocolour, format='D'),
       fits.Column(name='pseudocolour_error', array=pseudocolour_error, format='D'),
       fits.Column(name='nu_eff_used_in_astrometry', array=nu_eff_used_in_astrometry, format='D'),
       fits.Column(name='astrometric_params_solved', array=astrometric_params_solved, format='D'),
       fits.Column(name='ecl_lat', array=ecl_lat, format='D'),
       fits.Column(name='ruwe', array=ruwe, format='D'),
       fits.Column(name='ipd_gof_harmonic_amplitude', array=ipd_gof_harmonic_amplitude, format='D'),
       fits.Column(name='ipd_frac_multi_peak', array=ipd_frac_multi_peak, format='D'),
       fits.Column(name='phot_bp_rp_excess_factor', array=phot_bp_rp_excess_factor, format='D')]

if ZPT:  # official gaia zero-point
       col.append(fits.Column(name='parallax_w_zp', array=parallax-zp_row_matched, format='D'))


t = fits.BinTableHDU.from_columns(col)
t.writeto(gaia_rowmatch_f)