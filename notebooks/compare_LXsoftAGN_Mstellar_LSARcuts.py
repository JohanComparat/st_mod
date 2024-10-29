"""
python compare_LXsoftAGN_Mstellar_LSARcuts.py z0p05 sigma_1.0_fsat_10.0
python compare_LXsoftAGN_Mstellar_LSARcuts.py z0p09 sigma_1.0_fsat_10.0
python compare_LXsoftAGN_Mstellar_LSARcuts.py z0p14 sigma_1.0_fsat_10.0
python compare_LXsoftAGN_Mstellar_LSARcuts.py z0p19 sigma_1.0_fsat_10.0

python compare_LXsoftAGN_Mstellar_LSARcuts.py z0p05 sigma_0.8_fsat_8.0
python compare_LXsoftAGN_Mstellar_LSARcuts.py z0p09 sigma_0.8_fsat_8.0
python compare_LXsoftAGN_Mstellar_LSARcuts.py z0p14 sigma_0.8_fsat_8.0
python compare_LXsoftAGN_Mstellar_LSARcuts.py z0p19 sigma_0.8_fsat_8.0

python compare_LXsoftAGN_Mstellar_LSARcuts.py z0p05 sigma_0.4_fsat_0.0
python compare_LXsoftAGN_Mstellar_LSARcuts.py z0p09 sigma_0.4_fsat_0.0
python compare_LXsoftAGN_Mstellar_LSARcuts.py z0p14 sigma_0.4_fsat_0.0
python compare_LXsoftAGN_Mstellar_LSARcuts.py z0p19 sigma_0.4_fsat_0.0

python compare_LXsoftAGN_Mstellar_LSARcuts.py z0p25 sigma_1.0_fsat_10.0
python compare_LXsoftAGN_Mstellar_LSARcuts.py z0p25 sigma_0.8_fsat_8.0
python compare_LXsoftAGN_Mstellar_LSARcuts.py z0p25 sigma_0.4_fsat_0.0


"""
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

z_dir = sys.argv[1]
agn_mod = sys.argv[2]

uchuu_dir       = os.path.join(os.environ['UCHUU'], 'GPX8')
p2_glists = np.array( glob.glob( os.path.join( uchuu_dir, z_dir,'replication_0.0_0.0_0.0/glist.fits' ) ) )
z_str = np.array([ el.split('/')[-3] for el in p2_glists])
z_val = np.array([ float(el[1:].split('p')[0]+'.'+el[1:].split('p')[1]) for el in z_str])
p2_agnlists = np.array([ os.path.join( os.path.dirname( el ), 'AGN_list_'+agn_mod+'.fits' ) for el in p2_glists ])

validation_dir       = os.path.join(os.environ['GIT_STMOD'], 'data', 'validation','validation_AGN', 'LX_vs_Ms_fullpopulation')
fig_dir = validation_dir

mdex = 0.25
mbins = np.arange(8, 12.1, mdex)
x_m = mbins+mdex/2.

#z_min, z_max = 0.20, 0.28
#selection = (z_val>=z_min)&(z_val<=z_max)
#print('z_min, z_max=', z_min, z_max)

GAL = Table.read(p2_glists[0])
AGN = Table.read(p2_agnlists[0])
AGN['lx39']=10**(AGN['LX_soft']-39)
AGN['EddRatio']=25*10**(AGN['LX_hard']-38)/( 1.26 * 0.002 * AGN['obs_sm'] )

N_GAL = np.histogram(np.log10(GAL['obs_sm']), bins = mbins)[0]
N_AGN = np.histogram(np.log10(AGN['obs_sm']), bins = mbins)[0]
AGN['lx39']=10**(AGN['LX_soft']-39)
frac_AGN = N_AGN*1./N_GAL

agn_t = {}
for m_min in mbins[:-1]:
    s1 = (np.log10(AGN['obs_sm'])>=m_min) & (np.log10(AGN['obs_sm'])<m_min+mdex)
    agn_t[np.round(m_min,1)] = AGN[s1]

mean_lx39 = np.zeros_like(mbins)
mean_lx39_ERgtm4 = np.zeros_like(mbins)
mean_lx39_ERgtm3 = np.zeros_like(mbins)
mean_lx39_ERgtm2 = np.zeros_like(mbins)
for jj, m_min in enumerate(mbins[:-1]):
    mean_lx39[jj] = np.sum(agn_t[np.round(m_min,1)]['lx39']) / N_GAL[jj]
    mean_lx39_ERgtm4 [jj] = np.sum(agn_t[np.round(m_min,1)]['lx39'][agn_t[np.round(m_min,1)]['EddRatio']>2e-3]) / N_GAL[jj]
    mean_lx39_ERgtm3 [jj] = np.sum(agn_t[np.round(m_min,1)]['lx39'][agn_t[np.round(m_min,1)]['EddRatio']>4e-3]) / N_GAL[jj]
    mean_lx39_ERgtm2 [jj] = np.sum(agn_t[np.round(m_min,1)]['lx39'][agn_t[np.round(m_min,1)]['EddRatio']>1e-2]) / N_GAL[jj]


mean_lx39_LH40 = np.zeros_like(mbins)
mean_lx39_LH41 = np.zeros_like(mbins)
mean_lx39_LH42 = np.zeros_like(mbins)
for jj, m_min in enumerate(mbins[:-1]):
    mean_lx39_LH40 [jj] = np.sum(agn_t[np.round(m_min,1)]['lx39'][agn_t[np.round(m_min,1)]['LX_hard']>40]) / N_GAL[jj]
    mean_lx39_LH41 [jj] = np.sum(agn_t[np.round(m_min,1)]['lx39'][agn_t[np.round(m_min,1)]['LX_hard']>41]) / N_GAL[jj]
    mean_lx39_LH42 [jj] = np.sum(agn_t[np.round(m_min,1)]['lx39'][agn_t[np.round(m_min,1)]['LX_hard']>42]) / N_GAL[jj]

figure_name = os.path.join(fig_dir, z_dir+agn_mod+'-LXsoft-Mstar-LsarSelections.png')
plt.figure(1, (6,5))
plt.plot(x_m, np.log10(mean_lx39       )+39, label='all', ls='solid', lw=2)
plt.plot(x_m, np.log10(mean_lx39_ERgtm4)+39, label=r'$\lambda_{Edd}$>2e-3', ls='dotted')
plt.plot(x_m, np.log10(mean_lx39_ERgtm3)+39, label=r'$\lambda_{Edd}$>4e-3', ls='dotted')
plt.plot(x_m, np.log10(mean_lx39_ERgtm2)+39, label=r'$\lambda_{Edd}$>1e-2', ls='dotted')
plt.plot(x_m, np.log10(mean_lx39_LH40)+39, label=r'$L^{2-10}_X$>40', ls='dashed')
plt.plot(x_m, np.log10(mean_lx39_LH41)+39, label=r'$L^{2-10}_X$>41', ls='dashed')
plt.plot(x_m, np.log10(mean_lx39_LH42)+39, label=r'$L^{2-10}_X$>42', ls='dashed')

#plt.xscale('log')
#plt.yscale('log')
plt.xlabel(r"stellar mass [$M_\odot$]")
plt.ylabel(r"$L^{0.5-2}_{X, AGN}$ [erg s$^{-1}$]")
plt.title(agn_mod)
plt.legend(ncol=2, title='z='+z_dir[1:], fontsize=12)
plt.tight_layout()
plt.savefig( figure_name )
plt.clf()
print(figure_name, 'written')


figure_name = os.path.join(fig_dir, z_dir+agn_mod+'-ratio-LXsoft-Mstar-LsarSelections.png')
plt.figure(1, (6,5))
plt.plot(x_m, mean_lx39_ERgtm4/mean_lx39, label=r'$\lambda_{Edd}$>2e-3', ls='dotted', lw=3)
plt.plot(x_m, mean_lx39_ERgtm3/mean_lx39, label=r'$\lambda_{Edd}$>4e-3', ls='dotted', lw=3)
plt.plot(x_m, mean_lx39_ERgtm2/mean_lx39, label=r'$\lambda_{Edd}$>1e-2', ls='dotted', lw=3)
plt.plot(x_m, mean_lx39_LH40/mean_lx39, label=r'$L^{2-10}_X$>40', ls='dashed', lw=3)
plt.plot(x_m, mean_lx39_LH41/mean_lx39, label=r'$L^{2-10}_X$>41', ls='dashed', lw=3)
plt.plot(x_m, mean_lx39_LH42/mean_lx39, label=r'$L^{2-10}_X$>42', ls='dashed', lw=3)

#plt.xscale('log')
#plt.yscale('log')
plt.xlabel(r"stellar mass [$M_\odot$]")
plt.ylabel(r"frac = $L^{0.5-2}_{X, AGN}$/ Full population [erg s$^{-1}$]")
plt.title(agn_mod)
plt.legend(ncol=2, title='z='+z_dir[1:], fontsize=12)
plt.tight_layout()
plt.savefig( figure_name )
plt.clf()
print(figure_name, 'written')
