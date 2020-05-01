import os
import numpy as np
import pylab as plt
from astropy.io import fits
from config import allstar_path, astronn_apogee_vac_f
from galpy.orbit import Orbit
from matplotlib.colors import LogNorm


f_allstar = fits.getdata(allstar_path)
f_new_vac = fits.getdata(astronn_apogee_vac_f)

if os.path.exists("diagnostic_plots"):
    raise FileExistsError("Folder 'diagnostic_plots' already existed!")
else:
    os.mkdir("diagnostic_plots")

elements_ls = ['C', 'CI', 'N', 'O', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'K', 'Ca', 'Ti',
               'TiII', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni']


for elem in elements_ls:
    fig = plt.figure(dpi=200)
    ax = fig.gca()
    good_id = ((f_new_vac[f"Fe_H_ERR"] < 0.2) & (f_new_vac[f"{elem}_H_ERR"] < 0.2))
    elem_fe = f_new_vac[f"{elem}_H"]-f_new_vac[f"Fe_H"]
    ax.hexbin(f_new_vac[f"Fe_H"][good_id], elem_fe[good_id], extent=(-1.5, 0.5, -0.5, 0.5), bins='log', gridsize=(25, 25))
    ax.set_xlim(-1.5, 0.5)
    ax.set_ylim(-0.5, 0.5)
    ax.set_xlabel(f'[Fe/H]')
    ax.set_ylabel(f'[{elem}/Fe]')
    fig.tight_layout()
    fig.savefig(f"diagnostic_plots/{elem}_Fe.png")
    plt.close("all")

fig = plt.figure(dpi=200)
ax = fig.gca()
good_id = ((f_new_vac[f"TEFF_ERR"] < 100) & (f_new_vac[f"LOGG_ERR"] < 0.2) & (f_new_vac[f"FE_H_ERR"] < 0.2))
ax.scatter(f_new_vac[f"TEFF"][good_id], f_new_vac[f"LOGG"][good_id], c=f_new_vac[f"FE_H"][good_id], s=1e-2)
ax.set_xlim(6000, 3000)
ax.set_ylim(5, 0)
ax.set_xlabel(r'$T_\mathrm{eff}$')
ax.set_ylabel(r'$\log{g}$')
fig.tight_layout()
fig.savefig(f"diagnostic_plots/teff_logg_feh.png")
plt.close("all")

fig = plt.figure(dpi=200)
ax = fig.gca()
good_id = ((f_new_vac[f"Fe_H_ERR"] < 0.2) & (f_new_vac[f"MG_H_ERR"] < 0.2) &
           (f_new_vac[f"AGE_TOTAL_ERROR"]/f_new_vac[f"AGE"] < 0.4))
ax.scatter(f_new_vac[f"Fe_H"][good_id], f_new_vac[f"MG_H"][good_id]-f_new_vac[f"Fe_H"][good_id],
           c=f_new_vac[f"AGE"][good_id], s=1e-2, norm=LogNorm())
ax.set_xlim(-1.5, 0.5)
ax.set_ylim(-0.5, 0.5)
ax.set_xlabel(f'[Fe/H]')
ax.set_ylabel(f'[MG/Fe]')
fig.tight_layout()
fig.savefig(f"diagnostic_plots/alpha_age.png")
plt.close("all")


o= Orbit(np.array([f_new_vac['ra'],f_new_vac['dec'],f_new_vac['weighted_dist']/1000.,
                   f_new_vac['pmra'],f_new_vac['pmdec'],f_new_vac['vhelio_avg']]).T,
         ro=8.125,vo=220.,zo=0.0208,solarmotion=[-11.1, 25.7, 7.25],
         radec=True)

fig = plt.figure(dpi=200)
ax = fig.gca()
ax.scatter(o.R(), f_new_vac['GALR'], s=1e-2)
ax.plot([0, 20], [0, 20], c='k', ls='--')
ax.set_xlim(0.,20.)
ax.set_ylim(0.,20.)
ax.set_xlabel(r'$R\,\mathrm{Orbit}$')
ax.set_ylabel(r'$R\,\mathrm{catalog}$')
fig.tight_layout()
fig.savefig(f"diagnostic_plots/galr_check.png")
plt.close("all")


X = f_new_vac['GALR'] * np.cos(f_new_vac['GALPHI'])
Y = f_new_vac['GALR'] * np.sin(f_new_vac['GALPHI'])

fig = plt.figure(dpi=200)
ax = fig.gca()
ax.hist2d(X, Y, range=((-1,15),(-12.5,12.5)), bins=40, norm=LogNorm())
ax.set_xlabel('X')
ax.set_ylabel('Y')
fig.tight_layout()
fig.savefig(f"diagnostic_plots/galaxy_hist.png")
plt.close("all")
