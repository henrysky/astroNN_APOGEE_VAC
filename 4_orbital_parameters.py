# This script uses galpy to generate orbital parameters and save them to fits

import tqdm
import numpy as np
import multiprocessing
from astropy.io import fits
from galpy.util import bovy_coords, bovy_conversion
from galpy.actionAngle import actionAngleStaeckel, estimateDeltaStaeckel
from galpy.potential import MWPotential2014, evaluatePotentials, rl

from config import allstar_path, gaia_rowmatch_f, astronn_dist_f, galpy_orbitparams_f, _R0, _v0, _z0

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
distance = f_astronn['weighted_dist'] / 1e3
pmra = f_gaia['pmra']
pmdec = f_gaia['pmdec']
rv = allstar_data['vhelio_avg']

# extract relevant fields
ra_err = f_gaia['ra_error']
dec_err = f_gaia['dec_error']
parallax_err = f_gaia['parallax_error']
distance_err = f_astronn['weighted_dist_error'] / 1e3
pmra_err = f_gaia['pmra_error']
pmdec_err = f_gaia['pmdec_error']
rv_err = allstar_data['vscatter']
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


def obs_to_galcen(ra, dec, dist, pmra, pmdec, rv, ro=_R0, vo=_v0, zo=_z0):
    vxvv = np.dstack([ra, dec, dist, pmra, pmdec, rv])[0]
    ra, dec = vxvv[:, 0], vxvv[:, 1]
    lb = bovy_coords.radec_to_lb(ra, dec, degree=True)
    pmra, pmdec = vxvv[:, 3], vxvv[:, 4]
    pmllpmbb = bovy_coords.pmrapmdec_to_pmllpmbb(pmra, pmdec, ra, dec, degree=True)
    d, vlos = vxvv[:, 2], vxvv[:, 5]
    rectgal = bovy_coords.sphergal_to_rectgal(lb[:, 0], lb[:, 1], d, vlos, pmllpmbb[:, 0], pmllpmbb[:, 1], degree=True)
    vsolar = np.array([-11.1, 245.7, 7.25])
    vsun = vsolar / vo
    X = rectgal[:, 0] / ro
    Y = rectgal[:, 1] / ro
    Z = rectgal[:, 2] / ro
    vx = rectgal[:, 3] / vo
    vy = rectgal[:, 4] / vo
    vz = rectgal[:, 5] / vo
    XYZ = np.dstack([X, Y, Z])[0]
    vxyz = np.dstack([vx, vy, vz])[0]
    Rpz = bovy_coords.XYZ_to_galcencyl(XYZ[:, 0], XYZ[:, 1], XYZ[:, 2], Zsun=zo / ro)
    vRvTvz = bovy_coords.vxvyvz_to_galcencyl(vxyz[:, 0], vxyz[:, 1], vxyz[:, 2], Rpz[:, 0], Rpz[:, 1], Rpz[:, 2],
                                             vsun=vsun,
                                             Xsun=1.,
                                             Zsun=zo / ro,
                                             galcen=True)

    return XYZ, vxyz, Rpz, vRvTvz


def process_single(i):
    vxvv = [ra[i], dec[i], distance[i], pmra[i], pmdec[i], rv[i]]
    if not np.all(np.isfinite(vxvv)) or not np.all(np.isfinite(covariance[i])):
        return np.ones(41) * np.nan
    samp = np.random.multivariate_normal(vxvv, covariance[i], size=100)
    if not (np.all(samp[:, 1] < 90.) & np.all(samp[:, 1] > -90.)):
        return np.ones(41) * np.nan
    sXYZ, svxyz, sRpz, svRvTvz = obs_to_galcen(samp[:, 0],
                                               samp[:, 1],
                                               samp[:, 2],
                                               samp[:, 3],
                                               samp[:, 4],
                                               samp[:, 5],
                                               ro=_R0,
                                               vo=_v0,
                                               zo=_z0)
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
    E = (evaluatePotentials(MWPotential2014, sRpz[:, 0], sRpz[:, 2], phi=sRpz[:, 1]) + np.sum(svRvTvz ** 2 / 2.,
                                                                                              axis=1)) * _v0 ** 2

    return np.nanmean(e), e_err, np.nanmean(zmax) * _R0, zmax_err * _R0, np.nanmean(
        rperi) * _R0, rperi_err * _R0, np.nanmean(rap) * _R0, rap_err * _R0, \
           e_zmax_corr, e_rperi_corr, e_rap_corr, zmax_rperi_corr, zmax_rap_corr, rperi_rap_corr, \
           np.nanmean(action[0]) * _R0 * _v0, jr_err * _R0 * _v0, np.nanmean(
        action[1]) * _R0 * _v0, lz_err * _R0 * _v0, np.nanmean(action[2]) * _R0 * _v0, jz_err * _R0 * _v0, \
           jr_lz_corr, jr_jz_corr, lz_jz_corr, \
           np.nanmean(action[3]) * _freq, or_err * _freq, np.nanmean(action[4]) * _freq, op_err * _freq, np.nanmean(
        action[5]) * _freq, oz_err * _freq, \
           np.nanmean(action[6]), tr_err, np.nanmean(action[7]), tphi_err, np.nanmean(action[8]), tz_err, \
           np.nanmean(Rc), np.nanstd(Rc), np.nanmean(E), np.nanstd(E), np.nanmean(E - Ec), np.nanstd(E - Ec)


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
                                              ('EminusEc_err', float)])

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
            'EminusEc_err']
    rec['source_id'] = source_ids
    rec['APOGEE_ID'] = apogee_ids
    rec['ra'] = ra
    rec['dec'] = dec
    for i in range(len(keys)):
        rec[keys[i]] = output[:, i]

    hdu = fits.BinTableHDU.from_columns(rec)
    hdu.writeto(galpy_orbitparams_f)
