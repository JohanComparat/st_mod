
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as plt
import os, sys
import glob
import numpy as np
from astropy.table import Table, vstack
from scipy.interpolate import interp1d
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoUNIT = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
h = 0.6774
L_box = 2000.0 / h
cosmo = cosmoUNIT


t_GAMA = Table.read('/home/comparat/sf_Shared/data/GAMA/forJohan.fits')
kmag = 8.9 - 2.5*np.log10(t_GAMA['flux_Kt'])
zmag = 8.9 - 2.5*np.log10(t_GAMA['flux_Zt'])
rmag = 8.9 - 2.5*np.log10(t_GAMA['flux_rt'])
imag = 8.9 - 2.5*np.log10(t_GAMA['flux_it'])
keep_GAMA = (kmag>0) & (t_GAMA['Z']>0.01) & (t_GAMA['Z']<1.99) & (zmag>0) & (rmag>0)  & (imag>0)
GAMA = t_GAMA[keep_GAMA]
area_GAMA = 60.
completeness_GAMA = 0.98

validation_dir       = os.path.join(os.environ['GIT_STMOD'], 'data', 'validation','validation_GAL')
DDD = np.loadtxt( os.path.join(validation_dir, 'GAMA', 'driver-2022-SMF.txt'), unpack = True )
M_dr22,  phi_All, phi_All_err = DDD.T[(DDD[0]>9)].T
L15_ngal,   L15_M,     L15_phi,    L15_Err = np.loadtxt(os.path.join(validation_dir, 'GAMA', 'loveday_2015', 'smfs.txt'), unpack=True)

mdex = 0.1
mbins = np.arange(8, 12.1, mdex)

z_min, z_max = 0.1, 0.2
print('z_min, z_max=', z_min, z_max)

fig_dir = os.path.join(os.environ['GIT_STMOD'], 'data', 'benchmark', )
BM_dir = os.path.join(os.environ['GIT_STMOD'], 'data', 'benchmark', )


z_GAMA = (GAMA['Z']>z_min) & (GAMA['Z']<z_max)
volume_GAMA =  (( cosmo.comoving_volume(z_max) - cosmo.comoving_volume(z_min) ) * area_GAMA * np.pi / 129600.).value
SM_arr_GAMA = np.log10(GAMA['StellarMass_50'][z_GAMA])
sSFR_arr_GAMA = np.log10(GAMA['SFR_50'][z_GAMA])-SM_arr_GAMA
QU_GAMA = (sSFR_arr_GAMA<-11)
SF_GAMA = (sSFR_arr_GAMA>=-11)
#
SDSS = Table.read(os.path.join('/home/comparat/sf_Shared/data/SDSS/DR7.2','Ti20_full_MERGED.fits'))
z_SDSS = (SDSS['Z']>z_min) & (SDSS['Z']<z_max)
area_SDSS = 8030.
completeness_SDSS = 0.95
volume_SDSS =  (( cosmo.comoving_volume(z_max) - cosmo.comoving_volume(z_min) ) * area_SDSS * np.pi / 129600.).value
SM_arr_SDSS = SDSS['log10M_star'][z_SDSS]
QU_SDSS = (SDSS['color_flag'][z_SDSS]==1)
SF_SDSS = (SDSS['color_flag'][z_SDSS]==0)

N_all = np.histogram(SM_arr_GAMA, mbins)[0]
N_qu = np.histogram(SM_arr_GAMA[QU_GAMA], mbins)[0]
N_sf = np.histogram(SM_arr_GAMA[SF_GAMA], mbins)[0]
f_qu = interp1d(mbins[:-1]+mdex/2., 1.0 * N_qu / N_all)
f_sf = interp1d(mbins[:-1]+mdex/2., 1.0 * N_sf / N_all)
ER_qu = interp1d(mbins[:-1]+mdex/2., N_qu**(-0.5))
ER_sf = interp1d(mbins[:-1]+mdex/2., N_sf**(-0.5))
#
# SMF
#
fig_out = os.path.join(fig_dir, 'SMF_SDSS_GAMA_SF_QU.png' )
plt.figure(1, (4.5, 4.))
# literature SMF
plt.fill_between(M_dr22,  10**(phi_All-phi_All_err),  10**(phi_All+phi_All_err),
                 color='k', alpha=0.5)
plt.fill_between(M_dr22,
                 f_qu(M_dr22) * 10**(phi_All-phi_All_err)*(1-ER_qu(M_dr22)),
                 f_qu(M_dr22) * 10**(phi_All+phi_All_err)*(1+ER_qu(M_dr22)),
                 color='darkred', alpha=0.5)
plt.fill_between(M_dr22,
                 f_sf(M_dr22) * 10**(phi_All-phi_All_err)*(1-ER_sf(M_dr22)),
                 f_sf(M_dr22) * 10**(phi_All+phi_All_err)*(1+ER_sf(M_dr22)),
                 color='darkblue', alpha=0.5)

plt.plot(M_dr22,  10**phi_All , lw=2, color='k', ls='dashed')

#plt.plot(M_dr22,  10**phi_All * N_qu / N_all, lw=2, color='darkblue')# , label='Driver 2021, z=0') # +5*np.log10(h)
#plt.errorbar( L15_M, y =  L15_phi*0.7**3, yerr=   L15_Err*0.7**3, label='Loveday 15, z<0.65')
# literature DATA
# GAMA G09
#plt.hist(SM_arr_GAMA, lw=1, weights=np.ones_like(SM_arr_GAMA) / ( mdex * volume_GAMA * completeness_GAMA ) ,
                #bins = mbins, histtype = 'step', rasterized = True, label = 'GAMA data')
# SDSS MGS
#plt.hist(SM_arr_SDSS, lw=1, weights=np.ones_like(SM_arr_SDSS) / ( mdex * volume_SDSS * completeness_SDSS ) ,
                #bins = mbins, histtype = 'step', rasterized = True, label = 'SDSS data')

plt.ylabel(r'$\Phi$ [dex$^{-1}$ Mpc$^{-3}$]')
plt.xlabel('$\log_{10}$(stellar mass[M$_\odot$])')
plt.yscale('log')
plt.xlim(( 9.5, 11.5))
plt.ylim((1e-4, 0.01))
#plt.legend(loc=0, fontsize=12)#, frameon=False)
plt.tight_layout()
plt.savefig(fig_out)
plt.clf()
print(fig_out, 'written')


##figure_name = os.path.join(fig_dir, 'wprpFIG_M_106_110_BCRS.png')
##plt.figure(0, (4.5,4))
##x, y2, y1 = np.loadtxt( os.path.join(BM_dir, 'wprp_M_106_110_shift_4_BC.txt'), unpack = True)
##plt.plot(x, (y1/1e4+ y2/1e4)/2., color='darkblue', ls='dashed')
##plt.fill_between(x, y1/1e4, y2/1e4, color='darkblue', alpha=0.5)
##x, y2, y1 = np.loadtxt( os.path.join(BM_dir, 'wprp_M_106_110_shift_4_RS.txt'), unpack = True)
##plt.plot(x, (y1/1e4+ y2/1e4)/2., color='darkred', ls='dashed')
##plt.fill_between(x, y1/1e4, y2/1e4, color='darkred', alpha=0.5)
##plt.xscale('log')
##plt.yscale('log')
##plt.xlabel(r"$r_p$ [$h^{-1}$Mpc]")
##plt.ylabel(r"$w_p(r_p)$ [$h^{-1}$Mpc]")
##plt.title('Milky Way analog galaxies')
##plt.tight_layout()

##plt.savefig( figure_name )
##plt.clf()
##print(figure_name, 'written')
