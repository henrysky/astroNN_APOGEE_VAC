# This script uses galpy to generate orbital parameters and save them to fits

import os
import tqdm
import numpy as np
import multiprocessing
from astropy.io import fits
from galpy.orbit import Orbit
from galpy.util import bovy_conversion
from galpy.actionAngle import actionAngleStaeckel, estimateDeltaStaeckel
from galpy.potential import MWPotential2014, evaluatePotentials, rl

from config import allstar_path, gaia_rowmatch_f, astronn_dist_f, galpy_orbitparams_f, _R0, _v0, _z0

if os.path.exists(galpy_orbitparams_f):
    raise FileExistsError(f"{galpy_orbitparams_f} already existed")

allstar_data = fits.getdata(allstar_path)
f_gaia = fits.getdata(gaia_rowmatch_f)
f_dist = fits.getdata(astronn_dist_f)
f_astronn = fits.getdata(astronn_dist_f)

_freq = bovy_conversion.freq_in_Gyr(_v0, _R0)

# extract gaia info
source_ids = f_gaia['source_id']
apogee_ids = f_astronn['APOGEE_ID']
ra = f_gaia['ra']
dec = f_gaia['dec']
parallax = f_gaia['parallax']
distance = f_astronn['weighted_dist'] / 1e3  # pc to kpc
pmra = f_gaia['pmra']
pmdec = f_gaia['pmdec']
rv = allstar_data['vhelio_avg']

# extract relevant fields
ra_err = f_gaia['ra_error']
dec_err = f_gaia['dec_error']
parallax_err = f_gaia['parallax_error']
distance_err = f_astronn['weighted_dist_error'] / 1e3  # pc to kpc
pmra_err = f_gaia['pmra_error']
pmdec_err = f_gaia['pmdec_error']
rv_err = allstar_data['verr']
ra_dec_corr = f_gaia['ra_dec_corr']
ra_parallax_corr = f_gaia['ra_parallax_corr']
ra_pmra_corr = f_gaia['ra_pmra_corr']
ra_pmdec_corr = f_gaia['ra_pmdec_corr']
dec_parallax_corr = f_gaia['dec_parallax_corr']
dec_pmra_corr = f_gaia['dec_pmra_corr']
dec_pmdec_corr = f_gaia['dec_pmdec_corr']
parallax_pmra_corr = f_gaia['parallax_pmra_corr']
parallax_pmdec_corr = f_gaia['parallax_pmdec_corr']
pmra_pmdec_corr = f_gaia['pmra_pmdec_corr']

# build covariance matrices
covariance = np.zeros([len(ra), 6, 6])
covariance[:, 0, 0] = ra_err ** 2
covariance[:, 1, 1] = dec_err ** 2
covariance[:, 2, 2] = distance_err ** 2
covariance[:, 3, 3] = pmra_err ** 2
covariance[:, 4, 4] = pmdec_err ** 2
covariance[:, 5, 5] = rv_err ** 2
covariance[:, 0, 1] = ra_dec_corr * ra_err * dec_err
covariance[:, 0, 2] = 0.
covariance[:, 0, 3] = ra_pmra_corr * ra_err * pmra_err
covariance[:, 0, 4] = ra_pmdec_corr * ra_err * pmdec_err
covariance[:, 1, 0] = covariance[:, 0, 1]
covariance[:, 1, 2] = 0.
covariance[:, 1, 3] = dec_pmra_corr * dec_err * pmra_err
covariance[:, 1, 4] = dec_pmdec_corr * dec_err * pmdec_err
covariance[:, 2, 0] = covariance[:, 0, 2]
covariance[:, 2, 1] = covariance[:, 1, 2]
covariance[:, 2, 3] = 0.
covariance[:, 2, 4] = 0.
covariance[:, 3, 0] = covariance[:, 0, 3]
covariance[:, 3, 1] = covariance[:, 1, 3]
covariance[:, 3, 2] = covariance[:, 2, 3]
covariance[:, 3, 4] = pmra_pmdec_corr * pmra_err * pmdec_err
covariance[:, 4, 0] = covariance[:, 0, 4]
covariance[:, 4, 1] = covariance[:, 1, 4]
covariance[:, 4, 2] = covariance[:, 2, 4]
covariance[:, 4, 3] = covariance[:, 3, 4]


def process_single(i):
    vxvv = [ra[i], dec[i], distance[i], pmra[i], pmdec[i], rv[i]]
    if not np.all(np.isfinite(vxvv)) or not np.all(np.isfinite(covariance[i])):
        return np.ones(56) * np.nan
    samp = np.random.multivariate_normal(vxvv, covariance[i], size=1000)
    if not (np.all(samp[:, 1] < 90.) & np.all(samp[:, 1] > -90.)):
        return np.ones(56) * np.nan

    os = Orbit(np.array([samp[:, 0], samp[:, 1], samp[:, 2], samp[:, 3], samp[:, 4], samp[:, 5]]).T, radec=True,ro=_R0,
               vo=_v0, zo=_z0, solarmotion=[-11.1, 25.7, 7.25])

    sXYZ = np.dstack([os.x(), os.y(), os.z()])[0] / _R0
    sRpz = np.dstack([os.R() / _R0, os.phi(), os.z() / _R0])[0]
    svRvTvz = np.dstack([os.vR(), os.vT(), os.vz()])[0] / _v0

    deltas = estimateDeltaStaeckel(MWPotential2014, np.median(sRpz[:, 0]), np.median(sRpz[:, 2]), no_median=True)
    aAS = actionAngleStaeckel(pot=MWPotential2014, delta=np.mean(deltas))
    e, zmax, rperi, rap = aAS.EccZmaxRperiRap(sRpz[:, 0],
                                              svRvTvz[:, 0],
                                              svRvTvz[:, 1],
                                              sRpz[:, 2],
                                              svRvTvz[:, 2],
                                              sRpz[:, 1], delta=deltas)
    tcov = np.cov(np.dstack([e, zmax, rperi, rap])[0].T)
    errs = np.sqrt(np.diag(tcov))
    e_zmax_corr = tcov[0, 1] / (errs[0] * errs[1])
    e_rperi_corr = tcov[0, 2] / (errs[0] * errs[2])
    e_rap_corr = tcov[0, 3] / (errs[0] * errs[3])
    zmax_rperi_corr = tcov[1, 2] / (errs[1] * errs[2])
    zmax_rap_corr = tcov[1, 3] / (errs[1] * errs[3])
    rperi_rap_corr = tcov[2, 3] / (errs[2] * errs[3])
    e_err, zmax_err, rperi_err, rap_err = errs
    action = np.array(
        aAS.actionsFreqsAngles(sRpz[:, 0], svRvTvz[:, 0], svRvTvz[:, 1], sRpz[:, 2], svRvTvz[:, 2], sRpz[:, 1]))
    tcov_after_action = np.cov(action)
    errs = np.sqrt(np.diag(tcov_after_action))
    jr_lz_corr = tcov_after_action[0, 1] / (errs[0] * errs[1])
    jr_jz_corr = tcov_after_action[0, 2] / (errs[0] * errs[2])
    lz_jz_corr = tcov_after_action[1, 2] / (errs[1] * errs[2])
    jr_err, lz_err, jz_err, or_err, op_err, oz_err, tr_err, tphi_err, tz_err = errs

    Rc = np.array([rl(MWPotential2014, lz) for lz in action[1]]) * _R0
    Ec = (evaluatePotentials(MWPotential2014, Rc, 0.) + 0.5 * (action[1]) ** 2. / Rc ** 2.) * _v0 ** 2
    E = os.E(pot=MWPotential2014)

    # galactocentric coord and vel uncertainty
    galcen_tcov = np.cov(np.dstack([sRpz[:, 0], sRpz[:, 1], sRpz[:, 2]])[0].T)
    galcen_errs = np.sqrt(np.diag(galcen_tcov))
    galr_err, galphi_err, galz_err = galcen_errs

    galcenv_tcov = np.cov(np.dstack([svRvTvz[:, 0], svRvTvz[:, 1], svRvTvz[:, 2]])[0].T)
    galcenv_errs = np.sqrt(np.diag(galcenv_tcov))
    galvr_err, galvt_err, galvz_err = galcenv_errs

    galvr_galvt_corr = galcenv_tcov[0, 1] / (galcenv_errs[0] * galcenv_errs[1])
    galvr_galvz_corr = galcenv_tcov[0, 2] / (galcenv_errs[0] * galcenv_errs[2])
    galvt_galvz_corr = galcenv_tcov[1, 2] / (galcenv_errs[1] * galcenv_errs[2])

    # galr mean to avoid issue near GC, error propagation
    x_mean = np.nanmean(sXYZ[:, 0])
    y_mean = np.nanmean(sXYZ[:, 1])
    x_err = np.nanstd(sXYZ[:, 0])
    y_err = np.nanstd(sXYZ[:, 1])
    x_err_percentage = x_err / np.nanmean(sXYZ[:, 0])
    y_err_percentage = y_err / np.nanmean(sXYZ[:, 1])
    x2_err = (x_mean ** 2) * (2 * x_err_percentage)
    y2_err = (y_mean ** 2) * (2 * y_err_percentage)
    galr = np.sqrt(x_mean ** 2 + y_mean ** 2)
    galr2_err = np.sqrt(x2_err**2 + y2_err**2)
    galr2_err_percentage = galr2_err / (galr**2)
    galr_err = galr * (0.5 * galr2_err_percentage)

    # galphi mean to avoid issue near GC, error propagation
    galphi = np.arctan(y_mean/x_mean)
    galphi_err_x = (y_mean / (x_mean ** 2 + y_mean ** 2)) * x_err  # error propagation from x_mean
    galphi_err_y = (-x_mean / (x_mean ** 2 + y_mean ** 2)) * y_err  # error propagation from y_mean
    galphi_err = np.sqrt(galphi_err_x**2 + galphi_err_y**2)  # add them up

    return np.nanmean(e), \
           e_err, \
           np.nanmean(zmax) * _R0, \
           zmax_err * _R0, \
           np.nanmean(rperi) * _R0, \
           rperi_err * _R0, \
           np.nanmean(rap) * _R0, \
           rap_err * _R0, \
           e_zmax_corr, \
           e_rperi_corr, \
           e_rap_corr, \
           zmax_rperi_corr, \
           zmax_rap_corr, \
           rperi_rap_corr, \
           np.nanmean(action[0]) * _R0 * _v0, \
           jr_err * _R0 * _v0, \
           np.nanmean(action[1]) * _R0 * _v0, \
           lz_err * _R0 * _v0, \
           np.nanmean(action[2]) * _R0 * _v0, \
           jz_err * _R0 * _v0, \
           jr_lz_corr, \
           jr_jz_corr, \
           lz_jz_corr, \
           np.nanmean(action[3]) * _freq, or_err * _freq, \
           np.nanmean(action[4]) * _freq, op_err * _freq, \
           np.nanmean(action[5]) * _freq, oz_err * _freq, \
           np.nanmean(action[6]), \
           tr_err, \
           np.nanmean(action[7]), \
           tphi_err, \
           np.nanmean(action[8]), \
           tz_err, \
           np.nanmean(Rc), \
           np.nanstd(Rc), \
           np.nanmean(E), \
           np.nanstd(E), \
           np.nanmean(E - Ec), \
           np.nanstd(E - Ec), \
           galr * _R0, \
           galphi, \
           np.nanmean(sRpz[:, 2]) * _R0, \
           np.nanmean(svRvTvz[:, 0]) * _v0, \
           np.nanmean(svRvTvz[0, 1]) * _v0, \
           np.nanmean(svRvTvz[:, 2]) * _v0, \
           galr_err * _R0, \
           galphi_err, \
           galz_err * _R0, \
           galvr_err * _v0, \
           galvt_err * _v0, \
           galvz_err * _v0, \
           galvr_galvt_corr, \
           galvr_galvz_corr, \
           galvt_galvz_corr


if __name__ == '__main__':  # needed for multiprocessing on Windows
    print('starting MP run...')
    with multiprocessing.Pool(int(multiprocessing.cpu_count() / 2)) as p:
        output = list(tqdm.tqdm(p.imap(process_single, range(len(ra))), total=len(ra)))

    output = np.array(output)

    rec = np.recarray([len(output), ], dtype=[('source_id', np.int64),
                                              ('APOGEE_ID', '<U18'),
                                              ('ra', float),
                                              ('dec', float),
                                              ('e', float),
                                              ('e_err', float),
                                              ('zmax', float),
                                              ('zmax_err', float),
                                              ('rperi', float),
                                              ('rperi_err', float),
                                              ('rap', float),
                                              ('rap_err', float),
                                              ('e_zmax_corr', float),
                                              ('e_rperi_corr', float),
                                              ('e_rap_corr', float),
                                              ('zmax_rperi_corr', float),
                                              ('zmax_rap_corr', float),
                                              ('rperi_rap_corr', float),
                                              ('jr', float),
                                              ('jr_err', float),
                                              ('Lz', float),
                                              ('Lz_err', float),
                                              ('jz', float),
                                              ('jz_err', float),
                                              ('jr_Lz_corr', float),
                                              ('jr_jz_corr', float),
                                              ('lz_jz_corr', float),
                                              ('omega_r', float),
                                              ('omega_r_err', float),
                                              ('omega_phi', float),
                                              ('omega_phi_err', float),
                                              ('omega_z', float),
                                              ('omega_z_err', float),
                                              ('theta_r', float),
                                              ('theta_r_err', float),
                                              ('theta_phi', float),
                                              ('theta_phi_err', float),
                                              ('theta_z', float),
                                              ('theta_z_err', float),
                                              ('rl', float),
                                              ('rl_err', float),
                                              ('Energy', float),
                                              ('Energy_err', float),
                                              ('EminusEc', float),
                                              ('EminusEc_err', float),
                                              ('galr', float),
                                              ('galphi', float),
                                              ('galz', float),
                                              ('galvr', float),
                                              ('galvt', float),
                                              ('galvz', float),
                                              ('galr_err', float),
                                              ('galphi_err', float),
                                              ('galz_err', float),
                                              ('galvr_err', float),
                                              ('galvt_err', float),
                                              ('galvz_err', float),
                                              ('galvr_galvt_corr', float),
                                              ('galvr_galvz_corr', float),
                                              ('galvt_galvz_corr', float)])

    keys = ['e',
            'e_err',
            'zmax',
            'zmax_err',
            'rperi',
            'rperi_err',
            'rap',
            'rap_err',
            'e_zmax_corr',
            'e_rperi_corr',
            'e_rap_corr',
            'zmax_rperi_corr',
            'zmax_rap_corr',
            'rperi_rap_corr',
            'jr',
            'jr_err',
            'Lz',
            'Lz_err',
            'jz',
            'jz_err',
            'jr_Lz_corr',
            'jr_jz_corr',
            'lz_jz_corr',
            'omega_r',
            'omega_r_err',
            'omega_phi',
            'omega_phi_err',
            'omega_z',
            'omega_z_err',
            'theta_r',
            'theta_r_err',
            'theta_phi',
            'theta_phi_err',
            'theta_z',
            'theta_z_err',
            'rl',
            'rl_err',
            'Energy',
            'Energy_err',
            'EminusEc',
            'EminusEc_err',
            'galr',
            'galphi',
            'galz',
            'galvr',
            'galvt',
            'galvz',
            'galr_err',
            'galphi_err',
            'galz_err',
            'galvr_err',
            'galvt_err',
            'galvz_err',
            'galvr_galvt_corr',
            'galvr_galvz_corr',
            'galvt_galvz_corr']

    rec['source_id'] = source_ids
    rec['APOGEE_ID'] = apogee_ids
    rec['ra'] = ra
    rec['dec'] = dec
    for i in range(len(keys)):
        rec[keys[i]] = output[:, i]

    hdu = fits.BinTableHDU.from_columns(rec)
    hdu.writeto(galpy_orbitparams_f)
