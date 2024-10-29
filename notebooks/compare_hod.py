
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


mdex = 0.1
mbins = np.arange(8, 16.1, mdex)
xm = mbins[:-1]+mdex/2.
z_min, z_max = 0.02, 0.2
print('z_min, z_max=', z_min, z_max)

fig_dir = os.path.join(os.environ['GIT_STMOD'], 'data', 'benchmark', )
BM_dir = os.path.join(os.environ['GIT_STMOD'], 'data', 'benchmark', )
#
SDSS = Table.read(os.path.join('/home/comparat/sf_Shared/data/SDSS/DR7.2','Ti20_full_MERGED.fits'))
area_SDSS = 8030.
completeness_SDSS = 0.95
volume_SDSS =  (( cosmo.comoving_volume(z_max) - cosmo.comoving_volume(z_min) ) * area_SDSS * np.pi / 129600.).value
SM_arr_SDSS = SDSS['log10M_star']
Mhalo_arr_SDSS = np.log10(SDSS['M_halo'])
z_SEL = (SDSS['Z']>z_min) & (SDSS['Z']<z_max)
MWsel = (SM_arr_SDSS>10.5)&(SM_arr_SDSS<11)
sat = (SDSS['P_sat']>=0.5)
cen = (SDSS['P_sat']<0.5)
QU = (SDSS['color_flag']==1)
SF = (SDSS['color_flag']==0)

ALL_cen = np.histogram(Mhalo_arr_SDSS[cen & z_SEL ], mbins)[0]
MW_cen  = np.histogram(Mhalo_arr_SDSS[cen & z_SEL & MWsel], mbins)[0]
MW_cen_SF  = np.histogram(Mhalo_arr_SDSS[cen & z_SEL & MWsel & SF], mbins)[0]
MW_cen_QU  = np.histogram(Mhalo_arr_SDSS[cen & z_SEL & MWsel & QU], mbins)[0]
MW_sat_SF  = np.histogram(Mhalo_arr_SDSS[sat & z_SEL & MWsel & SF], mbins)[0]
MW_sat_QU  = np.histogram(Mhalo_arr_SDSS[sat & z_SEL & MWsel & QU], mbins)[0]


f_MW_cen    = MW_cen *1.    / (ALL_cen + 1 )
f_MW_cen_SF = MW_cen_SF *1. / (ALL_cen + 1 )
f_MW_cen_QU = MW_cen_QU *1. / (ALL_cen + 1 )
f_MW_sat_SF = MW_sat_SF *1. / (ALL_cen + 1 )
f_MW_sat_QU = MW_sat_QU *1. / (ALL_cen + 1 )

err_MW_cen    = (1./(MW_cen   +1)  + 1./(ALL_cen+1))**(0.5)
err_MW_cen_SF = (1./(MW_cen_SF+1)  + 1./(ALL_cen+1))**(0.5)
err_MW_cen_QU = (1./(MW_cen_QU+1)  + 1./(ALL_cen+1))**(0.5)
err_MW_sat_SF = (1./(MW_sat_SF+1)  + 1./(ALL_cen+1))**(0.5)
err_MW_sat_QU = (1./(MW_sat_QU+1)  + 1./(ALL_cen+1))**(0.5)

#
fig_out = os.path.join(fig_dir, 'Ti22_HOD_SF_QU.png' )
plt.figure(1, (4.5, 4.))
# literature SMF
#plt.fill_between(xm,
                 #f_MW_cen * (1-err_MW_cen),
                 #f_MW_cen*(1+err_MW_cen),
                 #color='k', alpha=0.5)
plt.fill_between(xm,
                 f_MW_cen_QU * (1-err_MW_cen_QU)+f_MW_sat_QU * (1-err_MW_sat_QU),
                 f_MW_cen_QU*(1+err_MW_cen_QU)+f_MW_sat_QU*(1+err_MW_sat_QU),
                 color='darkred', alpha=0.5)
plt.plot(xm, f_MW_cen_QU, color='darkred', ls='dashed')
plt.plot(xm, f_MW_sat_QU, color='darkred', ls='dotted')
plt.fill_between(xm,
                 f_MW_cen_SF * (1-err_MW_cen_SF)+f_MW_sat_SF * (1-err_MW_sat_SF),
                 f_MW_cen_SF*(1+err_MW_cen_SF)+f_MW_sat_SF*(1+err_MW_sat_SF),
                 color='darkblue', alpha=0.5)
plt.plot(xm, f_MW_cen_SF, color='darkblue', ls='dashed')
plt.plot(xm, f_MW_sat_SF, color='darkblue', ls='dotted')

plt.plot(xm-1000, f_MW_cen_SF, color='k', ls='dashed', label='central')
plt.plot(xm-1000, f_MW_sat_SF, color='k', ls='dotted', label='satellites')
plt.ylabel(r'$\langle N(M) \rangle$')
plt.xlabel('$\log_{10}$(halo mass [M$_\odot$])')
plt.yscale('log')
plt.xlim(( 11.1, 14.5))
plt.ylim((0.01,10))
plt.title('Halo occupation distribution')
plt.legend()
plt.tight_layout()
plt.savefig(fig_out)
plt.clf()
print(fig_out, 'written')

