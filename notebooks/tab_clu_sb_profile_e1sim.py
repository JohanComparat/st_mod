"""
-rwxrwx--- 1 root vboxsf 1,1G Feb  9 09:24 J_A+A_665_A78_agnev.dat.gz.fits.gz
-rwxrwx--- 1 root vboxsf  33G Feb  9 09:54 J_A+A_665_A78_agnin.dat.gz.fits.gz
-rwxrwx--- 1 root vboxsf 7,1G Feb  9 09:39 J_A+A_665_A78_bkgev.dat.gz.fits.gz
-rwxrwx--- 1 root vboxsf 155M Feb  9 09:23 J_A+A_665_A78_clustev.dat.gz.fits.gz
-rwxrwx--- 1 root vboxsf 1,5G Feb  9 09:24 J_A+A_665_A78_clustin.dat.gz.fits.gz
-rwxrwx--- 1 root vboxsf 216M Feb  9 09:23 J_A+A_665_A78_esass1b.dat.gz.fits.gz
-rwxrwx--- 1 root vboxsf 119M Feb  9 09:23 J_A+A_665_A78_esass3b.dat.gz.fits.gz
-rwxrwx--- 1 root vboxsf  85M Feb  9 09:23 J_A+A_665_A78_starsev.dat.gz.fits.gz
-rwxrwx--- 1 root vboxsf  70M Feb  9 09:23 J_A+A_665_A78_starsin.dat.gz.fits.gz

"""
import numpy as np
import numpy
import h5py
import healpy
import os, sys, glob
#dir_2pcf = os.path.join( os.environ['GIT_ER4W'], 'data', 'xcorr' )
from astropy.table import Table, Column, vstack, hstack
import astropy.units as u
from astropy.coordinates import SkyCoord
import time
t0 = time.time()
from scipy.stats import scoreatpercentile
from sklearn.neighbors import BallTree
deg_to_rad = np.pi/180.

nl = lambda selection : len(selection.nonzero()[0])
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoUNIT = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
h = 0.6774
L_box = 1000.0 / h
cosmo = cosmoUNIT

NSIDE = 128

z_min = 0.05
z_max = 0.6

# events
events = Table.read('/home/comparat/sf_Shared/data/UNIT_fA1i/twin_seppi_2022/J_A+A_665_A78_clustev.dat.gz.fits.gz')
print(len(events))
#area_pix = healpy.nside2pixarea(NSIDE, degrees = True)
selection = ( events['Signal'] >= 0.5 ) & ( events['Signal']<=2.0 )
events = events[selection]
events['RA'], events['DEC'] = events['RAdeg'], events['DEdeg']
coords = SkyCoord(events['RA'], events['DEC'], unit='deg', frame='icrs')
events['g_lat']  = coords.galactic.b.value
events['g_lon'] = coords.galactic.l.value
sel_events = (events['g_lon'] >180) & (abs(events['g_lat']) > 20)
#events['HPX'] = healpy.ang2pix(NSIDE, np.pi/2. - events['DEC']*np.pi/180. , events['RA']*np.pi/180. , nest=True)
events = events[sel_events]
print(len(events))


t = Table.read('/home/comparat/sf_Shared/data/UNIT_fA1i/twin_seppi_2022/J_A+A_665_A78_clustin.dat.gz.fits.gz')
t['RA'], t['DEC'] = t['RAdeg'], t['DEdeg']
coords = SkyCoord(t['RA'], t['DEC'], unit='deg', frame='icrs')
t['g_lat']  = coords.galactic.b.value
t['g_lon'] = coords.galactic.l.value
t['HPX'] = healpy.ang2pix(NSIDE, np.pi/2. - t['DEC']*np.pi/180. , t['RA']*np.pi/180. , nest=True)
sel_clu = ( t['redshift_S'] >= 0.01 ) & ( t['redshift_S'] <= 0.3 ) & ( t['LX_soft'] >= 42) & (t['g_lon'] >180) & (abs(t['g_lat']) > 20)
halo = t[sel_clu]

M500_sort = np.argsort(np.log10(t['HALO500c']))
t = t[M500_sort[::-1] ]

texp = t['TexpModel:4/20']

coord_CLU = deg_to_rad * np.transpose([t['DEC'], t['RA'] ])
Tree_CLU = BallTree(coord_CLU, metric='haversine')
coord_EVT = deg_to_rad * np.transpose([events['DEC'], events['RA'] ])
Tree_EVT = BallTree(coord_EVT, metric='haversine')

r_values = np.hstack(( 10**np.arange(1, 2, 0.5), 10**np.arange(2, 3, 0.2), 10**np.arange(3, 3.51, 0.1) ))
N_in_ap2 = [np.zeros_like(r_max_rad.value)]
#for vv in N_in_ap:
	#N_in_ap2.append(vv)
for r_value in r_values:
	r_max_rad = r_value/(t['R500ckpc']/t['R500carcmin']) /60 * np.pi/180
	N_in_ap2.append( Tree_EVT.query_radius(coord_CLU, r = r_max_rad.value, count_only = True) )

N_in_ap_lo = N_in_ap2[:-1]
N_in_ap_hi = N_in_ap2[1:]
N_max_apertuire = N_in_ap_hi[-1]
area = np.pi* ( r_values**2. - np.hstack((0., r_values[:-1]))**2. )

erg_per_ct = 1.602177e-12 * 1000 # at 1 keV
ARF_1keV = 2200. #cm2
#average_CTS_05_20 = LX_05_20 * texp * ARF_1keV / ( dL**2 * erg_per_ct )
dl2 = 4*np.pi * (cosmo.luminosity_distance(t['zS']).to(u.cm)**2).value
#N_in_ap * erg_per_ct / texp # erg/s/kpc2

SB_in_ap = []
for N_0, N_1, area_val in zip(N_in_ap_lo, N_in_ap_hi, area):
	SB_in_ap.append( (N_1-N_0) * erg_per_ct *dl2 / ( texp.value *ARF_1keV* area_val) )

profiles = np.transpose(SB_in_ap)

x_lo = np.hstack((0., r_values[:-1]))
x_hi = r_values
x_m = (x_hi+x_lo)/2.

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as plt

fig_dir = os.path.join(os.environ['GIT_STMOD'], 'data', 'models', 'model_GAS')


figure_name = os.path.join(fig_dir, 'e1CLU_profiles_allM_mean.png')
plt.figure(3, (5.5,5))
m_mins = np.array([13.0, 13.5, 14.0, 14.5])
m_maxs = np.array([13.5, 14.0, 14.5, 15.0])
for m0, m1 in zip(m_mins, m_maxs):
	selection = (np.log10(t['HALO500c'])>=m0) & ( np.log10(t['HALO500c'])<m1) & (t['zS']<0.3)& (t['zS']>0.05)&(texp>200)&(texp<500)&(N_max_apertuire>10)
	prf = profiles[selection]
	y_val = prf.mean(axis=0)
	#y_50 = scoreatpercentile(prf, 50, axis=0)
	#plt.errorbar(x_m, y_50, yerr=[y_50-y_32, y_68-y_50], xerr=[x_m-x_lo, x_hi-x_m], lw=3, color='darkblue', label=r'median $\pm1\sigma$')
	plt.plot(x_m, y_val, label=str(m0)+'-'+str(m1))
plt.legend(loc=0, fontsize=12)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$r$ [kpc]")
plt.ylabel(r"$S_X(r)$ [erg kpc$^{-2}$ s$^{-1}$]")
plt.title(r'Mean')
plt.tight_layout()
plt.savefig( figure_name )
plt.clf()
print(figure_name, 'written')



figure_name = os.path.join(fig_dir, 'e1CLU_profiles_allM_median.png')
plt.figure(3, (5.5,5))
m_mins = np.array([13.0, 13.5, 14.0, 14.5])
m_maxs = np.array([13.5, 14.0, 14.5, 15.0])
for m0, m1 in zip(m_mins, m_maxs):
	selection = (np.log10(t['HALO500c'])>=m0) & ( np.log10(t['HALO500c'])<m1) & (t['zS']<0.3)& (t['zS']>0.05)&(texp>200)&(texp<500)&(N_max_apertuire>10)
	prf = profiles[selection]
	#y_val = prf.mean(axis=0)
	y_50 = scoreatpercentile(prf, 50, axis=0)
	plt.plot(x_m, y_50, label=str(m0)+'-'+str(m1))
plt.legend(loc=0, fontsize=12)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$r$ [kpc]")
plt.ylabel(r"$S_X(r)$ [erg kpc$^{-2}$ s$^{-1}$]")
plt.title(r'Median')
plt.tight_layout()
plt.savefig( figure_name )
plt.clf()
print(figure_name, 'written')


selection = (np.log10(t['HALO500c'])>14.5) & ( np.log10(t['HALO500c'])<15) & (t['zS']<0.3)& (t['zS']>0.05)&(texp>200)&(texp<500)&(N_max_apertuire>10)

prf = profiles[selection]
y_val = prf.mean(axis=0)
y_std = prf.std(axis=0)
y_68 = scoreatpercentile(prf, 68, axis=0)
y_50 = scoreatpercentile(prf, 50, axis=0)
y_32 = scoreatpercentile(prf, 32, axis=0)
y_95 = scoreatpercentile(prf, 95, axis=0)
y_05 = scoreatpercentile(prf,  5, axis=0)

figure_name = os.path.join(fig_dir, 'e1CLU_profiles_145M150.png')
plt.figure(1, (5.5,5))
plt.errorbar(x_m, y_50, yerr=[y_50-y_32, y_68-y_50], xerr=[x_m-x_lo, x_hi-x_m], lw=3, color='darkblue', label=r'median $\pm1\sigma$')
plt.plot(x_m, y_val, color='blue', label='mean')
plt.legend(loc=0)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$r$ [kpc]")
plt.ylabel(r"$S_X(r)$ [erg kpc$^{-2}$ s$^{-1}$]")
plt.title(r'14.5$<\log_{10}(M_{500c})<15$')
plt.tight_layout()

plt.savefig( figure_name )
plt.clf()
print(figure_name, 'written')



selection = (np.log10(t['HALO500c'])>14) & ( np.log10(t['HALO500c'])<14.5) & (t['zS']<0.3)& (t['zS']>0.05)&(texp>200)&(texp<500)&(N_max_apertuire>10)

prf = profiles[selection]
y_val = prf.mean(axis=0)
y_std = prf.std(axis=0)
y_68 = scoreatpercentile(prf, 68, axis=0)
y_50 = scoreatpercentile(prf, 50, axis=0)
y_32 = scoreatpercentile(prf, 32, axis=0)
y_95 = scoreatpercentile(prf, 95, axis=0)
y_05 = scoreatpercentile(prf,  5, axis=0)

figure_name = os.path.join(fig_dir, 'e1CLU_profiles_140M145.png')
plt.figure(1, (5.5,5))
plt.errorbar(x_m, y_50, yerr=[y_50-y_32, y_68-y_50], xerr=[x_m-x_lo, x_hi-x_m], lw=3, color='darkblue', label=r'median $\pm1\sigma$')
plt.plot(x_m, y_val, color='blue', label='mean')
plt.legend(loc=0)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$r$ [kpc]")
plt.ylabel(r"$S_X(r)$ [erg kpc$^{-2}$ s$^{-1}$]")
plt.title(r'14$<\log_{10}(M_{500c})<14.5$')
plt.tight_layout()

plt.savefig( figure_name )
plt.clf()
print(figure_name, 'written')


selection = (np.log10(t['HALO500c'])>13.5) & ( np.log10(t['HALO500c'])<14.) & (t['zS']<0.3)& (t['zS']>0.05)&(texp>200)&(texp<500)&(N_max_apertuire>10)

prf = profiles[selection]
y_val = prf.mean(axis=0)
y_std = prf.std(axis=0)
y_68 = scoreatpercentile(prf, 68, axis=0)
y_50 = scoreatpercentile(prf, 50, axis=0)
y_32 = scoreatpercentile(prf, 32, axis=0)
y_95 = scoreatpercentile(prf, 95, axis=0)
y_05 = scoreatpercentile(prf,  5, axis=0)

figure_name = os.path.join(fig_dir, 'e1CLU_profiles_135M140.png')
plt.figure(1, (5.5,5))
plt.errorbar(x_m, y_50, yerr=[y_50-y_32, y_68-y_50], xerr=[x_m-x_lo, x_hi-x_m], lw=3, color='darkblue', label=r'median $\pm1\sigma$')
plt.plot(x_m, y_val, color='blue', label='mean')
plt.legend(loc=0)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$r$ [kpc]")
plt.ylabel(r"$S_X(r)$ [erg kpc$^{-2}$ s$^{-1}$]")
plt.title(r'13.5$<\log_{10}(M_{500c})<14$')
plt.tight_layout()

plt.savefig( figure_name )
plt.clf()
print(figure_name, 'written')


selection = (np.log10(t['HALO500c'])>13) & ( np.log10(t['HALO500c'])<13.5) & (t['zS']<0.3)& (t['zS']>0.05)&(texp>200)&(texp<500)&(N_max_apertuire>10)

prf = profiles[selection]
y_val = prf.mean(axis=0)
y_std = prf.std(axis=0)
y_68 = scoreatpercentile(prf, 68, axis=0)
y_50 = scoreatpercentile(prf, 50, axis=0)
y_32 = scoreatpercentile(prf, 32, axis=0)
y_95 = scoreatpercentile(prf, 95, axis=0)
y_05 = scoreatpercentile(prf,  5, axis=0)

figure_name = os.path.join(fig_dir, 'e1CLU_profiles_130M135.png')
plt.figure(1, (5.5,5))
plt.errorbar(x_m, y_50, yerr=[y_50-y_32, y_68-y_50], xerr=[x_m-x_lo, x_hi-x_m], lw=3, color='darkblue', label=r'median $\pm1\sigma$')
plt.plot(x_m, y_val, color='blue', label='mean')
plt.legend(loc=0)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$r$ [kpc]")
plt.ylabel(r"$S_X(r)$ [erg kpc$^{-2}$ s$^{-1}$]")
plt.title(r'13$<\log_{10}(M_{500c})<13.5$')
plt.tight_layout()

plt.savefig( figure_name )
plt.clf()
print(figure_name, 'written')
