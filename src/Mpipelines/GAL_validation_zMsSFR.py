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

validation_dir_Uband = os.path.join(validation_dir, 'GP_UbandLuminosityFunction')
validation_dir_Gband = os.path.join(validation_dir, 'GP_GbandLuminosityFunction')
validation_dir_Rband = os.path.join(validation_dir, 'GP_RbandLuminosityFunction')
validation_dir_Iband = os.path.join(validation_dir, 'GP_IbandLuminosityFunction')
validation_dir_Zband = os.path.join(validation_dir, 'GP_ZbandLuminosityFunction')
validation_dir_Yband = os.path.join(validation_dir, 'GP_YbandLuminosityFunction')
validation_dir_Jband = os.path.join(validation_dir, 'GP_JbandLuminosityFunction')
validation_dir_Hband = os.path.join(validation_dir, 'GP_HbandLuminosityFunction')
validation_dir_Kband = os.path.join(validation_dir, 'GP_KbandLuminosityFunction')

os.system('mkdir -p ' + validation_dir       )

os.system('mkdir -p ' + validation_dir_Uband )
os.system('mkdir -p ' + validation_dir_Gband )
os.system('mkdir -p ' + validation_dir_Rband )
os.system('mkdir -p ' + validation_dir_Iband )
os.system('mkdir -p ' + validation_dir_Zband )
os.system('mkdir -p ' + validation_dir_Yband )
os.system('mkdir -p ' + validation_dir_Jband )
os.system('mkdir -p ' + validation_dir_Hband )
os.system('mkdir -p ' + validation_dir_Kband )

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
#kmag_abs = 8.9 - 2.5*np.log10(GAMA['flux_Kt']) - ( dm_itp(GAMA['Z']) + kcorr_itp(GAMA['Z']) )
#rmag_abs = 8.9 - 2.5*np.log10(GAMA['flux_rt']) - dm_itp(GAMA['Z'])
#GAMA['Kmag_abs'] = kmag_abs
#GAMA['Rmag_abs'] = rmag_abs
GAMA['Umag_abs'] = 8.9 - 2.5*np.log10(GAMA['flux_ut']) - dm_itp(GAMA['Z'])
GAMA['Gmag_abs'] = 8.9 - 2.5*np.log10(GAMA['flux_gt']) - dm_itp(GAMA['Z'])
GAMA['Rmag_abs'] = 8.9 - 2.5*np.log10(GAMA['flux_rt']) - dm_itp(GAMA['Z'])
GAMA['Imag_abs'] = 8.9 - 2.5*np.log10(GAMA['flux_it']) - dm_itp(GAMA['Z'])
GAMA['Zmag_abs'] = 8.9 - 2.5*np.log10(GAMA['flux_Zt']) - dm_itp(GAMA['Z'])
GAMA['Ymag_abs'] = 8.9 - 2.5*np.log10(GAMA['flux_Yt']) - dm_itp(GAMA['Z'])
GAMA['Jmag_abs'] = 8.9 - 2.5*np.log10(GAMA['flux_Jt']) - dm_itp(GAMA['Z'])
GAMA['Hmag_abs'] = 8.9 - 2.5*np.log10(GAMA['flux_Ht']) - dm_itp(GAMA['Z'])
GAMA['Kmag_abs'] = 8.9 - 2.5*np.log10(GAMA['flux_Kt']) - ( dm_itp(GAMA['Z']) + kcorr_itp(GAMA['Z']) )

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

area_KIDS = 180.
completeness_KIDS = 1.0

KIDS['Umag_abs'] = 8.9 - 2.5*np.log10(KIDS['flux_ut']) - dm_itp(KIDS['z_peak'])
KIDS['Gmag_abs'] = 8.9 - 2.5*np.log10(KIDS['flux_gt']) - dm_itp(KIDS['z_peak'])
KIDS['Rmag_abs'] = 8.9 - 2.5*np.log10(KIDS['flux_rt']) - dm_itp(KIDS['z_peak'])
KIDS['Imag_abs'] = 8.9 - 2.5*np.log10(KIDS['flux_it']) - dm_itp(KIDS['z_peak'])
KIDS['Zmag_abs'] = 8.9 - 2.5*np.log10(KIDS['flux_Zt']) - dm_itp(KIDS['z_peak'])
KIDS['Ymag_abs'] = 8.9 - 2.5*np.log10(KIDS['flux_Yt']) - dm_itp(KIDS['z_peak'])
KIDS['Jmag_abs'] = 8.9 - 2.5*np.log10(KIDS['flux_Jt']) - dm_itp(KIDS['z_peak'])
KIDS['Hmag_abs'] = 8.9 - 2.5*np.log10(KIDS['flux_Ht']) - dm_itp(KIDS['z_peak'])
KIDS['Kmag_abs'] = 8.9 - 2.5*np.log10(KIDS['flux_Kt']) - ( dm_itp(KIDS['z_peak']) + kcorr_itp(KIDS['z_peak']) )


##
#
#
# tabulate smfs predicted by simulation
#
#
##
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
    mock_area = meta['mean_area']
    if LC_dir =='LC1800':
        mock_area = 1800.
    if LC_dir =='LC0060':
        mock_area = 60.
    if LC_dir =='LC0002':
        mock_area = 2.
    print('z_min, z_max=', z_min, z_max, 'area=', mock_area )
    volume_mock = (( cosmo.comoving_volume(z_max) - cosmo.comoving_volume(z_min) ) * mock_area * np.pi / 129600.).value
    z_mean = np.mean(GAL['redshift_S'])
    #
    z_COSMOS = (COSMOS['photoz']>z_min) & (COSMOS['photoz']<z_max)
    volume_COSMOS =  (( cosmo.comoving_volume(z_max) - cosmo.comoving_volume(z_min) ) * area_COSMOS * np.pi / 129600.).value
    magr_arr_COSMOS = COSMOS['R'][z_COSMOS]
    #
    z_GAMA = (GAMA['Z']>z_min) & (GAMA['Z']<z_max)
    volume_GAMA =  (( cosmo.comoving_volume(z_max) - cosmo.comoving_volume(z_min) ) * area_GAMA * np.pi / 129600.).value
    magr_arr_GAMA = 8.9 - 2.5*np.log10(GAMA['flux_rt'][z_GAMA])
    #
    z_SDSS = (SDSS['Z']>z_min) & (SDSS['Z']<z_max)
    volume_SDSS =  (( cosmo.comoving_volume(z_max) - cosmo.comoving_volume(z_min) ) * area_SDSS * np.pi / 129600.).value
    MR_arr_SDSS = SDSS['Mag_r'][z_SDSS]
    #
    z_KIDS = (KIDS['z_peak']>z_min) & (KIDS['z_peak']<z_max)
    volume_KIDS =  (( cosmo.comoving_volume(z_max) - cosmo.comoving_volume(z_min) ) * area_KIDS * np.pi / 129600.).value
    magr_arr_KIDS = 8.9 - 2.5*np.log10(KIDS['flux_rt'][z_KIDS])

    #
    # abs_mag_LF
    #
    def plot_LF_gr(
            mag_name = 'R'
            , mock_mag = MAG['r_GP']-MAG['dist_mod']
            , Mr_arr_COSMOS = COSMOS['MR'][z_COSMOS]
            , MR_arr_GAMA = GAMA['Rmag_abs'][z_GAMA]
            , MR_arr_KIDS = KIDS['Rmag_abs'][z_KIDS]
            , MR_arr_SDSS = SDSS['Mag_r'][z_SDSS]
            , topdir = validation_dir_Rband
                    ):
        os.system('mkdir -p ' + os.path.join( topdir, z_dir ))
        fig_out = os.path.join(topdir, z_dir, LC_dir+'_LF_galaxies_replication_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']) +'.png' )
        plt.figure(1, (5., 5.))
        #plt.plot(Mk_m_5logh_z0 - 5 * np.log10(h), 10**log10_phi_h3_Mpcm3_magm1 * 0.7**3, lw=3 , label='Driver 2012, z=0') # +5*np.log10(h)
        #
        # literature DATA
        # COSMOS
        if z_mean<4:
            plt.hist(Mr_arr_COSMOS, lw=1, weights=np.ones_like(Mr_arr_COSMOS) / ( kdex * volume_COSMOS * completeness_COSMOS ) ,
                    bins = kbins, histtype = 'step', rasterized = True, label = 'COSMOS data')
        # GAMA G09
        if z_mean<0.4:
            plt.hist(MR_arr_GAMA, lw=1, weights=np.ones_like(MR_arr_GAMA) / ( kdex * volume_GAMA * completeness_GAMA ) ,
                    bins = kbins, histtype = 'step', rasterized = True, label = 'GAMA data')
        # KIDS
        if z_mean<0.6 and z_mean>0.1:
            plt.hist(MR_arr_KIDS, lw=1, weights=np.ones_like(MR_arr_KIDS) / ( kdex * volume_KIDS * completeness_KIDS ) ,
                    bins = kbins, histtype = 'step', rasterized = True, label = 'KIDS data')
        # SDSS
        if z_mean<0.25:
            plt.hist(MR_arr_SDSS, lw=1, weights=np.ones_like(MR_arr_SDSS) / ( kdex * volume_SDSS * completeness_SDSS ) ,
                    bins = kbins, histtype = 'step', rasterized = True, label = 'SDSS data')
        # Mock catalog
        plt.hist(mock_mag , lw=2, weights=np.ones_like(mock_mag) / ( kdex * volume_mock ) ,
                bins = kbins, histtype = 'step', rasterized = True, label = 'Mock', color='grey')

        #
        plt.title(r'$\bar{z}$=' + str(np.round(z_mean, 3)))
        plt.ylabel(r'$\Phi/\; mag^{-1}\; Mpc^{-3}$')
        plt.xlabel('Absolute magnitude '+mag_name)
        plt.yscale('log')
        plt.xlim(( -19, -26 ))
        plt.ylim((1e-6, 0.005))
        plt.legend(loc=0, fontsize=12)#, frameon=False)
        plt.tight_layout()
        plt.savefig(fig_out)
        plt.clf()
        print(fig_out, 'written')

    plot_LF_gr(
        mag_name = 'R'
        , mock_mag = MAG['r_GP']-MAG['dist_mod']
        , Mr_arr_COSMOS = COSMOS['MR'][z_COSMOS]
        , MR_arr_GAMA = GAMA['Rmag_abs'][z_GAMA]
        , MR_arr_KIDS = KIDS['Rmag_abs'][z_KIDS]
        , MR_arr_SDSS = SDSS['Mag_r'][z_SDSS]
        , topdir = validation_dir_Rband)
    plot_LF_gr(
        mag_name = 'G'
        , mock_mag = MAG['g_GP']-MAG['dist_mod']
        , Mr_arr_COSMOS = COSMOS['MG'][z_COSMOS]
        , MR_arr_GAMA = GAMA['Gmag_abs'][z_GAMA]
        , MR_arr_KIDS = KIDS['Gmag_abs'][z_KIDS]
        , MR_arr_SDSS = SDSS['Mag_g'][z_SDSS]
        , topdir = validation_dir_Gband)


    def plot_LF_uizyjhk(
            mag_name = 'R'
            , mock_mag = MAG['r_GP']-MAG['dist_mod']
            , Mr_arr_COSMOS = COSMOS['MR'][z_COSMOS]
            , MR_arr_GAMA = GAMA['Rmag_abs'][z_GAMA]
            , MR_arr_KIDS = KIDS['Rmag_abs'][z_KIDS]
            , topdir = validation_dir_Rband
                    ):
        os.system('mkdir -p ' + os.path.join( topdir, z_dir ))
        fig_out = os.path.join(topdir, z_dir, LC_dir+'_LF_galaxies_replication_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']) +'.png' )
        plt.figure(1, (5., 5.))
        #plt.plot(Mk_m_5logh_z0 - 5 * np.log10(h), 10**log10_phi_h3_Mpcm3_magm1 * 0.7**3, lw=3 , label='Driver 2012, z=0') # +5*np.log10(h)
        #
        # literature DATA
        # COSMOS
        if z_mean<4:
            plt.hist(Mr_arr_COSMOS, lw=1, weights=np.ones_like(Mr_arr_COSMOS) / ( kdex * volume_COSMOS * completeness_COSMOS ) ,
                    bins = kbins, histtype = 'step', rasterized = True, label = 'COSMOS data')
        # GAMA G09
        if z_mean<0.4:
            plt.hist(MR_arr_GAMA, lw=1, weights=np.ones_like(MR_arr_GAMA) / ( kdex * volume_GAMA * completeness_GAMA ) ,
                    bins = kbins, histtype = 'step', rasterized = True, label = 'GAMA data')
        # KIDS
        if z_mean<0.6 and z_mean>0.1:
            plt.hist(MR_arr_KIDS, lw=1, weights=np.ones_like(MR_arr_KIDS) / ( kdex * volume_KIDS * completeness_KIDS ) ,
                    bins = kbins, histtype = 'step', rasterized = True, label = 'KIDS data')

        # Mock catalog
        plt.hist(mock_mag , lw=2, weights=np.ones_like(mock_mag) / ( kdex * volume_mock ) ,
                bins = kbins, histtype = 'step', rasterized = True, label = 'Mock', color='grey')

        #
        plt.title(r'$\bar{z}$=' + str(np.round(z_mean, 3)))
        plt.ylabel(r'$\Phi/\; mag^{-1}\; Mpc^{-3}$')
        plt.xlabel('Absolute magnitude '+mag_name)
        plt.yscale('log')
        plt.xlim(( -19, -26 ))
        plt.ylim((1e-6, 0.005))
        plt.legend(loc=0, fontsize=12)#, frameon=False)
        plt.tight_layout()
        plt.savefig(fig_out)
        plt.clf()
        print(fig_out, 'written')

    plot_LF_uizyjhk(
            mag_name = 'U'
            , mock_mag = MAG['u_GP']-MAG['dist_mod']
            , Mr_arr_COSMOS = COSMOS['MU'][z_COSMOS]
            , MR_arr_GAMA = GAMA['Umag_abs'][z_GAMA]
            , MR_arr_KIDS = KIDS['Umag_abs'][z_KIDS]
            , topdir = validation_dir_Uband
                    )

    plot_LF_uizyjhk(
            mag_name = 'I'
            , mock_mag = MAG['i_GP']-MAG['dist_mod']
            , Mr_arr_COSMOS = COSMOS['MI'][z_COSMOS]
            , MR_arr_GAMA = GAMA['Imag_abs'][z_GAMA]
            , MR_arr_KIDS = KIDS['Imag_abs'][z_KIDS]
            , topdir = validation_dir_Iband
                    )

    plot_LF_uizyjhk(
            mag_name = 'Z'
            , mock_mag = MAG['z_GP']-MAG['dist_mod']
            , Mr_arr_COSMOS = COSMOS['MZ'][z_COSMOS]
            , MR_arr_GAMA = GAMA['Zmag_abs'][z_GAMA]
            , MR_arr_KIDS = KIDS['Zmag_abs'][z_KIDS]
            , topdir = validation_dir_Zband
                    )

    #plot_LF_uizyjhk(
            #mag_name = 'Y'
            #, mock_mag = MAG['y_GP']-MAG['dist_mod']
            #, Mr_arr_COSMOS = COSMOS['MY'][z_COSMOS]
            #, MR_arr_GAMA = GAMA['Ymag_abs'][z_GAMA]
            #, MR_arr_KIDS = KIDS['Ymag_abs'][z_KIDS]
            #, topdir = validation_dir_Yband
                    #)

    plot_LF_uizyjhk(
            mag_name = 'J'
            , mock_mag = MAG['j_GP']-MAG['dist_mod']
            , Mr_arr_COSMOS = COSMOS['MJ'][z_COSMOS]
            , MR_arr_GAMA = GAMA['Jmag_abs'][z_GAMA]
            , MR_arr_KIDS = KIDS['Jmag_abs'][z_KIDS]
            , topdir = validation_dir_Jband
                    )

    #plot_LF_uizyjhk(
            #mag_name = 'H'
            #, mock_mag = MAG['h_GP']-MAG['dist_mod']
            #, Mr_arr_COSMOS = COSMOS['MH'][z_COSMOS]
            #, MR_arr_GAMA = GAMA['Hmag_abs'][z_GAMA]
            #, MR_arr_KIDS = KIDS['Hmag_abs'][z_KIDS]
            #, topdir = validation_dir_Hband
                    #)

    plot_LF_uizyjhk(
            mag_name = 'K'
            , mock_mag = MAG['k_GP']-MAG['dist_mod']
            , Mr_arr_COSMOS = COSMOS['MK'][z_COSMOS]
            , MR_arr_GAMA = GAMA['Kmag_abs'][z_GAMA]
            , MR_arr_KIDS = KIDS['Kmag_abs'][z_KIDS]
            , topdir = validation_dir_Kband
                    )
    ##
    ## log N log R
    ##
    #os.system('mkdir -p ' + os.path.join( validation_dir_Rband, z_dir ))
    #fig_out = os.path.join(validation_dir_Rband, z_dir, LC_dir+'_logNlogR_galaxies_replication_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']) +'.png' )
    #plt.figure(1, (5., 5.))
    ##plt.plot(Mk_m_5logh_z0 - 5 * np.log10(h), 10**log10_phi_h3_Mpcm3_magm1 * 0.7**3, lw=3 , label='Driver 2012, z=0') # +5*np.log10(h)
    ##
    ## literature DATA
    ## COSMOS
    #if z_mean<4 and z_mean>0.5:
        #plt.hist(magr_arr_COSMOS, lw=1, weights=np.ones_like(magr_arr_COSMOS) / (  area_COSMOS * completeness_COSMOS ) ,
                 #bins = rbins, histtype = 'step', rasterized = True, label = 'COSMOS data', cumulative = True)
    ## GAMA G09
    #if z_mean<0.4:
        #plt.hist(magr_arr_GAMA, lw=1, weights=np.ones_like(magr_arr_GAMA) / (  area_GAMA * completeness_GAMA ) ,
                 #bins = rbins, histtype = 'step', rasterized = True, label = 'GAMA data', cumulative = True)
    ## KIDS
    #if z_mean<0.6 and z_mean>0.1:
        #plt.hist(magr_arr_KIDS, lw=1, weights=np.ones_like(magr_arr_KIDS) / (  area_KIDS * completeness_KIDS ) ,
                 #bins = rbins, histtype = 'step', rasterized = True, label = 'KIDS data', cumulative = True)
    ## SDSS
    ##if z_mean<0.25:
        ##plt.hist(magr_arr_SDSS, lw=1, weights=np.ones_like(magr_arr_SDSS) / ( area_SDSS * completeness_SDSS ) ,
                 ##bins = mbins, histtype = 'step', rasterized = True, label = 'SDSS data', cumulative = True)
    ## Mock catalog
    #plt.hist(MAG['r_GP'], lw=2, weights=np.ones_like(MAG['r_GP']) / (  meta['mean_area'] ) ,
            #bins = rbins, histtype = 'step', rasterized = True, label = 'Mock', color='grey', cumulative = True)

    ##
    #plt.title(r'$\bar{z}$=' + str(np.round(z_mean, 3)))
    #plt.ylabel(r'$N(>m)\; deg^{-2}$')
    #plt.xlabel('Observed magnitude R')
    #plt.yscale('log')
    #plt.xlim(( 12, 26 ))
    ##plt.ylim((1e-6, 0.005))
    #plt.legend(loc=0, fontsize=12)#, frameon=False)
    #plt.tight_layout()
    #plt.savefig(fig_out)
    #plt.clf()
    #print(fig_out, 'written')

