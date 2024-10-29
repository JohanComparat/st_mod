"""
Validation of the galaxy model
Figures :
 * stellar mass function
 * luminosity functions
 * color histogram (bimodality)
"""
import time
t0 = time.time()
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as plt
import sys, os, glob

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

#import extinction

from scipy.interpolate import interp1d

from astropy.table import Table

import numpy as np

print('Plots LF of all bands')
print('------------------------------------------------')
print('------------------------------------------------')

cosmoUNIT = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
h = 0.6774
L_box = 2000.0 / h
cosmo = cosmoUNIT

nl = lambda sel : len(sel.nonzero()[0])
zs = np.arange(0.0000001, 7.1, 0.001)
dm_itp = interp1d(zs, cosmo.distmod(zs).value)

validation_dir       = os.path.join(os.environ['GIT_STMOD_DATA'], 'data', 'validation','validation_GAL')
#validation_dir_Stell = os.path.join(validation_dir, 'StellarMassFunction')
#validation_dir_Uband = os.path.join(validation_dir, 'UbandLuminosityFunction')
#validation_dir_Gband = os.path.join(validation_dir, 'GbandLuminosityFunction')
validation_dir_Rband = os.path.join(validation_dir, 'GP_RbandLuminosityFunction')
#validation_dir_Iband = os.path.join(validation_dir, 'IbandLuminosityFunction')
#validation_dir_Zband = os.path.join(validation_dir, 'ZbandLuminosityFunction')
#validation_dir_Kband = os.path.join(validation_dir, 'KbandLuminosityFunction')

os.system('mkdir -p ' + validation_dir       )
#os.system('mkdir -p ' + validation_dir_Stell )
#os.system('mkdir -p ' + validation_dir_Uband )
#os.system('mkdir -p ' + validation_dir_Gband )
os.system('mkdir -p ' + validation_dir_Rband )
#os.system('mkdir -p ' + validation_dir_Iband )
#os.system('mkdir -p ' + validation_dir_Zband )
#os.system('mkdir -p ' + validation_dir_Kband )

#
# GAL MODEL
#
sys.path.append( os.path.join(os.environ['GIT_STMOD'], 'src') )
from models import GAL as GG

z_dir = sys.argv[1]
LC_dir = sys.argv[2] # 'FullSky'
C_GAL = GG.GAL(z_dir, LC_dir=LC_dir)
#LC_dir='FullSky'

##
#
#
# Data used in the model
#
#
##

t = Table.read(C_GAL.path_2_COSMOS)
good = (t['photoz']>0.01 )&( t['photoz']< 6. ) & ( t['K']>0 )& ( t['MK']<0 )&( t['MK']>-40 )&( t['mass_med']<14 )&( t['mass_med']>6 )
COSMOS = t[good]
area_COSMOS = 2.
completeness_COSMOS = 1.

##
KCORR_DATA = np.loadtxt( os.path.join( os.environ['GIT_STMOD_DATA'], 'data', 'models','model_GAL', 'VISTA_Ks_kcorrections.txt'), unpack = True)
kcorr_itp = interp1d(KCORR_DATA[0], KCORR_DATA[3])

##
t_GAMA = Table.read(C_GAL.path_2_GAMA)
kmag = 8.9 - 2.5*np.log10(t_GAMA['flux_Kt'])
zmag = 8.9 - 2.5*np.log10(t_GAMA['flux_Zt'])
rmag = 8.9 - 2.5*np.log10(t_GAMA['flux_rt'])
imag = 8.9 - 2.5*np.log10(t_GAMA['flux_it'])
keep_GAMA = (kmag>0) & (t_GAMA['Z']>0.01) & (t_GAMA['Z']<1.99) & (zmag>0) & (rmag>0)  & (imag>0)
GAMA = t_GAMA[keep_GAMA]
area_GAMA = 60.
completeness_GAMA = 0.98
kmag_abs = 8.9 - 2.5*np.log10(GAMA['flux_Kt']) - ( dm_itp(GAMA['Z']) + kcorr_itp(GAMA['Z']) )
rmag_abs = 8.9 - 2.5*np.log10(GAMA['flux_rt']) - dm_itp(GAMA['Z'])
GAMA['Kmag_abs'] = kmag_abs
GAMA['Rmag_abs'] = rmag_abs

##
SDSS = Table.read(C_GAL.path_2_SDSS)
area_SDSS = 8030.
completeness_SDSS = 0.95

##
t_KIDS = Table.read(C_GAL.path_2_KIDS)
kmag = 8.9 - 2.5*np.log10(t_KIDS['flux_Kt'])
zmag = 8.9 - 2.5*np.log10(t_KIDS['flux_Zt'])
rmag = 8.9 - 2.5*np.log10(t_KIDS['flux_rt'])
imag = 8.9 - 2.5*np.log10(t_KIDS['flux_it'])
keep_KIDS = (kmag>0) & (t_KIDS['z_peak']>0.01) & (t_KIDS['z_peak']<1.99) & (zmag>0) & (rmag>0)  & (imag>0)
KIDS = t_KIDS[keep_KIDS]
kmag_abs = 8.9 - 2.5*np.log10(KIDS['flux_Kt']) - ( dm_itp(KIDS['z_peak']) + kcorr_itp(KIDS['z_peak']) )
KIDS['Kmag_abs'] = kmag_abs
rmag_abs = 8.9 - 2.5*np.log10(KIDS['flux_rt']) - dm_itp(KIDS['z_peak'])
KIDS['Rmag_abs'] = rmag_abs
area_KIDS = 180.
completeness_KIDS = 1.0

##
#
#
# Literature measurements
#
#
##

# SDSS

# GAMA
Mk_m_5logh_z0, log10_phi_h3_Mpcm3_magm1 = np.loadtxt( os.path.join(validation_dir, 'GAMA', 'driver-2012-K-LF.txt'), unpack = True )
Msun_h70m2,  phi_All, phi_All_err = np.loadtxt( os.path.join(validation_dir, 'GAMA', 'driver-2022-SMF.txt'), unpack = True )
L15_ngal,   L15_M,     L15_phi,    L15_Err = np.loadtxt(os.path.join(validation_dir, 'GAMA', 'loveday_2015', 'smfs.txt'), unpack=True)

# COSMOS
smf_ilbert13 = lambda M, M_star, phi_1s, alpha_1s, phi_2s, alpha_2s : ( phi_1s * (M/M_star) ** alpha_1s + phi_2s * (M/M_star) ** alpha_2s ) * np.e ** (-M/M_star) * (M/ M_star)

path_ilbert13_SMF = os.path.join(validation_dir, 'COSMOS', "ilbert_2013_mass_function_params.txt")
zmin, zmax, N, M_comp, M_star, phi_1s, alpha_1s, phi_2s, alpha_2s, log_rho_s = np.loadtxt(path_ilbert13_SMF, unpack=True)

smf_ilbert_fun = np.array([
lambda mass : smf_ilbert13( mass , 10**M_star[0], phi_1s[0]*10**(-3), alpha_1s[0], phi_2s[0]*10**(-3), alpha_2s[0] )
, lambda mass : smf_ilbert13( mass , 10**M_star[1], phi_1s[1]*10**(-3), alpha_1s[1], phi_2s[1]*10**(-3), alpha_2s[1] )
, lambda mass : smf_ilbert13( mass , 10**M_star[2], phi_1s[2]*10**(-3), alpha_1s[2], phi_2s[2]*10**(-3), alpha_2s[2] )
, lambda mass : smf_ilbert13( mass , 10**M_star[3], phi_1s[3]*10**(-3), alpha_1s[3], phi_2s[3]*10**(-3), alpha_2s[3] )
, lambda mass : smf_ilbert13( mass , 10**M_star[4], phi_1s[4]*10**(-3), alpha_1s[4], phi_2s[4]*10**(-3), alpha_2s[4] )
, lambda mass : smf_ilbert13( mass , 10**M_star[5], phi_1s[5]*10**(-3), alpha_1s[5], phi_2s[5]*10**(-3), alpha_2s[5] )
, lambda mass : smf_ilbert13( mass , 10**M_star[6], phi_1s[6]*10**(-3), alpha_1s[6], phi_2s[6]*10**(-3), alpha_2s[6] )
, lambda mass : smf_ilbert13( mass , 10**M_star[7], phi_1s[7]*10**(-3), alpha_1s[7], phi_2s[7]*10**(-3), alpha_2s[7] )
])


smf_ilbert_zmin = np.array([
0.2
, 0.5
, 0.8
, 1.1
, 1.5
, 2.0
, 2.5
, 3.0 ])

smf_ilbert_zmax = np.array([
0.5
, 0.8
, 1.1
, 1.5
, 2.0
, 2.5
, 3.0
, 4.0 ])

smf_ilbert_name = np.array([ "Ilbert 2013 "+str(zmin)+"<z<"+str(zmax) for zmin, zmax in zip(smf_ilbert_zmin,smf_ilbert_zmax) ])


##
#
#
# tabulate smfs predicted by simulation
#
#
##
mdex = 0.1
mbins = np.arange(8, 12.1, mdex)
kdex = 0.1
kbins = np.arange(-30, -10, kdex)
rbins = np.arange(12, 26, kdex)

enough_area = (C_GAL.LC_MetaData['area_DC_max']>=0.5*np.max(C_GAL.LC_MetaData['area_DC_max'])) & (C_GAL.LC_MetaData['area_DC_max']>0)
C_GAL.LC_MetaData['mean_area'] = (C_GAL.LC_MetaData['area_DC_min'] + C_GAL.LC_MetaData['area_DC_max']) / 2.
small_difference_minmax_1 = ( C_GAL.LC_MetaData['area_DC_min'] / C_GAL.LC_MetaData['mean_area'] >= 0.8 ) & ( C_GAL.LC_MetaData['area_DC_min'] / C_GAL.LC_MetaData['mean_area'] <= 1.2 )
small_difference_minmax_2 = ( C_GAL.LC_MetaData['area_DC_max'] / C_GAL.LC_MetaData['mean_area'] >= 0.8 ) & ( C_GAL.LC_MetaData['area_DC_max'] / C_GAL.LC_MetaData['mean_area'] <= 1.2 )
for meta in C_GAL.LC_MetaData[(enough_area)&(small_difference_minmax_1)&(small_difference_minmax_2)]:
    #
    print(meta)
    # retrieve the resulting catalogues and meta data
    p_2_catalogue = os.path.join(os.environ['UCHUU'], LC_dir, z_dir, 'replication_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']), 'glist.fits')
    p_2_catal_MAG = os.path.join(os.environ['UCHUU'], LC_dir, z_dir, 'replication_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']), 'zMsSFRmatch_mags.fits')
    GAL = Table.read(p_2_catalogue)
    MAG = Table.read(p_2_catal_MAG)
    z_min, z_max = np.min(GAL['redshift_S']), np.max(GAL['redshift_S'])
    print('z_min, z_max=', z_min, z_max, 'area=', meta['mean_area'] )
    volume_mock = (( cosmo.comoving_volume(z_max) - cosmo.comoving_volume(z_min) ) * meta['mean_area'] * np.pi / 129600.).value
    z_mean = np.mean(GAL['redshift_S'])
    #
    z_COSMOS = (COSMOS['photoz']>z_min) & (COSMOS['photoz']<z_max)
    volume_COSMOS =  (( cosmo.comoving_volume(z_max) - cosmo.comoving_volume(z_min) ) * area_COSMOS * np.pi / 129600.).value
    SM_arr_COSMOS = COSMOS['mass_med'][z_COSMOS]
    MK_arr_COSMOS = COSMOS['MK'][z_COSMOS]
    Mr_arr_COSMOS = COSMOS['MR'][z_COSMOS]
    magr_arr_COSMOS = COSMOS['R'][z_COSMOS]
    #
    z_GAMA = (GAMA['Z']>z_min) & (GAMA['Z']<z_max)
    volume_GAMA =  (( cosmo.comoving_volume(z_max) - cosmo.comoving_volume(z_min) ) * area_GAMA * np.pi / 129600.).value
    SM_arr_GAMA = np.log10(GAMA['StellarMass_50'][z_GAMA])
    MK_arr_GAMA = GAMA['Kmag_abs'][z_GAMA]
    MR_arr_GAMA = GAMA['Rmag_abs'][z_GAMA]
    magr_arr_GAMA = 8.9 - 2.5*np.log10(GAMA['flux_rt'][z_GAMA])
    #
    z_SDSS = (SDSS['Z']>z_min) & (SDSS['Z']<z_max)
    volume_SDSS =  (( cosmo.comoving_volume(z_max) - cosmo.comoving_volume(z_min) ) * area_SDSS * np.pi / 129600.).value
    SM_arr_SDSS = SDSS['log10M_star'][z_SDSS]
    MR_arr_SDSS = SDSS['Mag_r'][z_SDSS]
    #
    z_KIDS = (KIDS['z_peak']>z_min) & (KIDS['z_peak']<z_max)
    volume_KIDS =  (( cosmo.comoving_volume(z_max) - cosmo.comoving_volume(z_min) ) * area_KIDS * np.pi / 129600.).value
    MK_arr_KIDS = KIDS['Kmag_abs'][z_KIDS]
    MR_arr_KIDS = KIDS['Rmag_abs'][z_KIDS]
    magr_arr_KIDS = 8.9 - 2.5*np.log10(KIDS['flux_rt'][z_KIDS])

    #
    # RLF
    #
    os.system('mkdir -p ' + os.path.join( validation_dir_Rband, z_dir ))
    fig_out = os.path.join(validation_dir_Rband, z_dir, LC_dir+'_RLF_galaxies_replication_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']) +'.png' )
    plt.figure(1, (5., 5.))
    #plt.plot(Mk_m_5logh_z0 - 5 * np.log10(h), 10**log10_phi_h3_Mpcm3_magm1 * 0.7**3, lw=3 , label='Driver 2012, z=0') # +5*np.log10(h)
    #
    # literature DATA
    # COSMOS
    if z_mean<4:
        plt.hist(Mr_arr_COSMOS, lw=1, weights=np.ones_like(Mr_arr_COSMOS) / ( mdex * volume_COSMOS * completeness_COSMOS ) ,
                 bins = kbins, histtype = 'step', rasterized = True, label = 'COSMOS data')
    # GAMA G09
    if z_mean<0.4:
        plt.hist(MR_arr_GAMA, lw=1, weights=np.ones_like(MR_arr_GAMA) / ( mdex * volume_GAMA * completeness_GAMA ) ,
                 bins = kbins, histtype = 'step', rasterized = True, label = 'GAMA data')
    # KIDS
    if z_mean<0.6 and z_mean>0.1:
        plt.hist(MR_arr_KIDS, lw=1, weights=np.ones_like(MR_arr_KIDS) / ( mdex * volume_KIDS * completeness_KIDS ) ,
                 bins = kbins, histtype = 'step', rasterized = True, label = 'KIDS data')
    # SDSS
    if z_mean<0.25:
        plt.hist(MR_arr_SDSS, lw=1, weights=np.ones_like(MR_arr_SDSS) / ( mdex * volume_SDSS * completeness_SDSS ) ,
                 bins = kbins, histtype = 'step', rasterized = True, label = 'SDSS data')
    # Mock catalog
    plt.hist(MAG['r_GP'] - dm_itp(GAL['redshift_S']), lw=2, weights=np.ones_like(MAG['r_GP']) / ( mdex * volume_mock ) ,
            bins = kbins, histtype = 'step', rasterized = True, label = 'Mock', color='grey')

    #
    plt.title(r'$\bar{z}$=' + str(np.round(z_mean, 3)))
    plt.ylabel(r'$\Phi/\; mag^{-1}\; Mpc^{-3}$')
    plt.xlabel('Absolute magnitude R')
    plt.yscale('log')
    plt.xlim(( -19, -26 ))
    plt.ylim((1e-6, 0.005))
    plt.legend(loc=0, fontsize=12)#, frameon=False)
    plt.tight_layout()
    plt.savefig(fig_out)
    plt.clf()
    print(fig_out, 'written')


    #
    # log N log R
    #
    os.system('mkdir -p ' + os.path.join( validation_dir_Rband, z_dir ))
    fig_out = os.path.join(validation_dir_Rband, z_dir, LC_dir+'_logNlogR_galaxies_replication_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']) +'.png' )
    plt.figure(1, (5., 5.))
    #plt.plot(Mk_m_5logh_z0 - 5 * np.log10(h), 10**log10_phi_h3_Mpcm3_magm1 * 0.7**3, lw=3 , label='Driver 2012, z=0') # +5*np.log10(h)
    #
    # literature DATA
    # COSMOS
    if z_mean<4 and z_mean>0.5:
        plt.hist(magr_arr_COSMOS, lw=1, weights=np.ones_like(magr_arr_COSMOS) / (  area_COSMOS * completeness_COSMOS ) ,
                 bins = rbins, histtype = 'step', rasterized = True, label = 'COSMOS data', cumulative = True)
    # GAMA G09
    if z_mean<0.4:
        plt.hist(magr_arr_GAMA, lw=1, weights=np.ones_like(magr_arr_GAMA) / (  area_GAMA * completeness_GAMA ) ,
                 bins = rbins, histtype = 'step', rasterized = True, label = 'GAMA data', cumulative = True)
    # KIDS
    if z_mean<0.6 and z_mean>0.1:
        plt.hist(magr_arr_KIDS, lw=1, weights=np.ones_like(magr_arr_KIDS) / (  area_KIDS * completeness_KIDS ) ,
                 bins = rbins, histtype = 'step', rasterized = True, label = 'KIDS data', cumulative = True)
    # SDSS
    if z_mean<0.25:
        plt.hist(magr_arr_SDSS, lw=1, weights=np.ones_like(magr_arr_SDSS) / ( area_SDSS * completeness_SDSS ) ,
                 bins = mbins, histtype = 'step', rasterized = True, label = 'SDSS data', cumulative = True)
    # Mock catalog
    plt.hist(MAG['r_GP'], lw=2, weights=np.ones_like(MAG['r_GP']) / (  meta['mean_area'] ) ,
            bins = rbins, histtype = 'step', rasterized = True, label = 'Mock', color='grey', cumulative = True)

    #
    plt.title(r'$\bar{z}$=' + str(np.round(z_mean, 3)))
    plt.ylabel(r'$N(>m)\; deg^{-2}$')
    plt.xlabel('Observed magnitude R')
    plt.yscale('log')
    plt.xlim(( 12, 26 ))
    #plt.ylim((1e-6, 0.005))
    plt.legend(loc=0, fontsize=12)#, frameon=False)
    plt.tight_layout()
    plt.savefig(fig_out)
    plt.clf()
    print(fig_out, 'written')

