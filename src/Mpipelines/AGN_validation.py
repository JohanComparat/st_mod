"""
Validation of the AGN model
Figures :
 * stellar mass function
 * luminosity functions
 * logN-logS
"""
import time
t0 = time.time()
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as plt
import sys, os, glob
from scipy.integrate import quad

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

#
# AGN MODEL
#
sys.path.append( os.path.join(os.environ['GIT_STMOD'], 'src') )
from models import AGN as GG

z_dir = sys.argv[1]
LC_dir = sys.argv[2] # 'FullSky'
C_AGN = GG.AGN(z_dir, LC_dir=LC_dir)
#LC_dir='FullSky'

str_scatter_0 = '0.8'
str_fsat = '8.0'
agn_data_dir = os.path.join( os.environ['GIT_STMOD_DATA'], 'data', 'validation/validation_AGN/literature_data' )


##

validation_dir           = os.path.join(os.environ['GIT_STMOD_DATA'], 'data', 'validation','validation_AGN')
validation_dir_Stell     = os.path.join(validation_dir, 'StellarMassFunction', z_dir)
validation_dir_XLF_hard  = os.path.join(validation_dir, 'XrayhardLuminosityFunction', z_dir)
validation_dir_XLF_soft  = os.path.join(validation_dir, 'XraySoftLuminosityFunction', z_dir)
validation_dir_lNlS_hard = os.path.join(validation_dir, 'XrayhardLogNlogS', z_dir)
validation_dir_lNlS_soft = os.path.join(validation_dir, 'XraySoftLogNlogS', z_dir)
validation_dir_Rband     = os.path.join(validation_dir, 'RbandLuminosityFunction', z_dir)
validation_dir_LSAR     = os.path.join(validation_dir, 'XrayLSAR', z_dir)
validation_dir_DC     = os.path.join(validation_dir, 'DutyCycle', z_dir)
#validation_dir_Kband     = os.path.join(validation_dir, 'KbandLuminosityFunction')

os.system('mkdir -p ' + validation_dir            )
os.system('mkdir -p ' + validation_dir_Stell      )
os.system('mkdir -p ' + validation_dir_XLF_hard   )
os.system('mkdir -p ' + validation_dir_XLF_soft   )
os.system('mkdir -p ' + validation_dir_lNlS_hard  )
os.system('mkdir -p ' + validation_dir_lNlS_soft  )
os.system('mkdir -p ' + validation_dir_Rband      )
os.system('mkdir -p ' + validation_dir_LSAR      )
os.system('mkdir -p ' + validation_dir_DC      )

##

f_duty = C_AGN.f_duty

# redshift binning
#DZ = z_max - z_min
#all_zs = np.arange(z_min, z_max, DZ*.99)

# LF binning
dlog_lx = 0.05*6
lx_bins = np.arange(36, 48, dlog_lx)
x_lx = lx_bins[:-1] + dlog_lx / 2.

# LSAR histogram
dlog_LSAR = 0.1
LSAR_bins = np.arange(30, 38, dlog_LSAR)
x_LSAR = LSAR_bins[:-1] + dlog_LSAR / 2.

# SMF and duty cycle
dlogM = 0.1
bins_SMF = np.arange(8, 13, dlogM)
x_SMF = (bins_SMF[1:] + bins_SMF[:-1]) * 0.5

# rLF
dlogR = 0.05
bins_rLF = np.arange(-30, -8, dlogM)
x_rLF = (bins_rLF[1:] + bins_rLF[:-1]) * 0.5



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

enough_area =  (C_AGN.LC_MetaData['area_DC_max']>0)  & (C_AGN.LC_MetaData['area_DC_min']>0)
#(C_AGN.LC_MetaData['area_DC_max']>=0.5*np.max(C_AGN.LC_MetaData['area_DC_max']))
C_AGN.LC_MetaData['mean_area'] = (C_AGN.LC_MetaData['area_DC_min'] + C_AGN.LC_MetaData['area_DC_max']) / 2.
small_difference_minmax_1 = ( C_AGN.LC_MetaData['area_DC_min'] / C_AGN.LC_MetaData['mean_area'] >= 0.8 ) & ( C_AGN.LC_MetaData['area_DC_min'] / C_AGN.LC_MetaData['mean_area'] <= 1.2 )
small_difference_minmax_2 = ( C_AGN.LC_MetaData['area_DC_max'] / C_AGN.LC_MetaData['mean_area'] >= 0.8 ) & ( C_AGN.LC_MetaData['area_DC_max'] / C_AGN.LC_MetaData['mean_area'] <= 1.2 )

for meta in C_AGN.LC_MetaData[(enough_area)]: # &(small_difference_minmax_1)&(small_difference_minmax_2)]:
    #
    print(meta)
    # retrieve the resulting catalogues and meta data
    p_2_catalogue = os.path.join(os.environ['UCHUU'], LC_dir, z_dir, 'replication_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']), 'glist.fits')
    p_2_catal_MAG = os.path.join(os.environ['UCHUU'], LC_dir, z_dir, 'replication_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']), 'Kmatch_mags.fits')
    p_2_catal_AGNs = np.array( glob.glob( os.path.join( os.path.dirname(p_2_catalogue), 'AGN_list_sigma_'+str_scatter_0+'_fsat_'+str_fsat+'.fits' ) ) )
    print(not os.path.isfile(p_2_catalogue))
    print(not os.path.isfile(p_2_catal_MAG))
    print(len(p_2_catal_AGNs)==0)
    if not os.path.isfile(p_2_catalogue) or not os.path.isfile(p_2_catal_MAG) or len(p_2_catal_AGNs)==0 :
        continue
    #
    sky_frac_Dmax = meta['area_DC_max'] / C_AGN.LC_MetaData['area_DC_max'].sum()
    C_AGN.sky_frac = sky_frac_Dmax * meta['N_frac_'+LC_dir] # ( sky_frac_Dmin + sky_frac_Dmax ) / 2.
    #
    # create a dictionnary of catalogs
    AGN_cat_names = {}
    for p_2_catal_AGN in p_2_catal_AGNs:
        bn = os.path.basename(p_2_catal_AGN)[9:-5]
        AGN_cat_names[p_2_catal_AGN] = bn
    # loop over
    #p_2_catal_AGN = os.path.join( os.path.dirname(p_2_catalogue), 'AGN_list_sigma_'+str_scatter_0+'_fsat_'+str_fsat+'.fits' )

    GAL = Table.read(p_2_catalogue)
    MAG = Table.read(p_2_catal_MAG)

    z_mean = np.mean(GAL['redshift_S'])
    z_min, z_max = np.round(z_mean, 2), np.round(z_mean, 2)+0.01
    #z_min, z_max = np.min(GAL['redshift_S']), np.max(GAL['redshift_S'])
    print('z_min, z_max=', z_min, z_max)
    volume_mock = (( cosmo.comoving_volume(z_max) - cosmo.comoving_volume(z_min) ) * C_AGN.sky_frac).value
    # galaxy data
    GAL = GAL[(GAL['redshift_S']>=z_min)&(GAL['redshift_S']<z_max)]
    MAG = MAG[(MAG['redshift_S']>=z_min)&(MAG['redshift_S']<z_max)]
    logm_gal = np.log10(GAL['sm'])
    z_gal = GAL['redshift_S']
    rmag_gal = MAG['rmag']
    kmagABS_gal = MAG['K_mag_abs']

    AGNs = {}
    for p_2_catal_AGN in p_2_catal_AGNs:
        t_i = Table.read(p_2_catal_AGN)
        AGNs[AGN_cat_names[p_2_catal_AGN]] = t_i[(t_i['redshift_S']>=z_min)&(t_i['redshift_S']<z_max)]


    # hard X-ray luminosity function Aird 2015)
    def kz_h(z): return 10**(-4.03 - 0.19 * (1 + z))
    def Ls_h(z): return 10**(44.84 - np.log10(((1 + 2.0) / (1 + z))** 3.87 + ((1 + 2.0) / (1 + z))**(-2.12)))
    def phi_h(L, z): return kz_h(z) / ((L / Ls_h(z))**0.48 + (L / Ls_h(z))**2.27)

    fig_out = os.path.join(validation_dir_XLF_hard, LC_dir+'_XLF_hard_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']) +'.png' )
    plt.figure(1, (6, 6))
    plt.axes([0.2, 0.18, 0.75, 0.75])
    plt.plot(x_lx, phi_h(10**x_lx, z_mean), c='cyan', ls='dashed', lw=2, label='Ai15')  # Aird 2-10 keV LADE')

    for p_2_catal_AGN in p_2_catal_AGNs:
        AGN = AGNs[AGN_cat_names[p_2_catal_AGN]]
        z = AGN['redshift_S']
        AGN = AGN[(z>=z_min)&(z<z_max)]
        z = AGN['redshift_S']
        lx = AGN['LX_hard']
        logNH = AGN['logNH']
        #fx = AGN['FX_soft']
        #lx_0520 = AGN['LX_soft']
        #logm = np.log10(AGN['obs_sm'])
        #lsar = lx - logm
        #mag_r = AGN['SDSS_r_AB']

        # Selections for the histogram
        NH20 = (logNH < 22)
        NH22 = (logNH >= 22) & (logNH < 24)
        NH24 = (logNH >= 24)

        N_nh20 = np.histogram(lx[NH20], lx_bins)[0] / volume_mock / dlog_lx
        N_nh22 = np.histogram(lx[NH22], lx_bins)[0] / volume_mock / dlog_lx
        N_nh24 = np.histogram(lx[NH24], lx_bins)[0] / volume_mock / dlog_lx

        nhar = np.histogram(lx, lx_bins)[0] / volume_mock / dlog_lx
        nharN = np.histogram(lx, lx_bins)[0]  # /vol/dlog_lx
        plt.fill_between(x_lx,
                        y1=(nhar) * (1 - (nharN)**(-0.5)),
                        y2=(nhar) * (1 + (nharN)**(-0.5)),
                        #color='g',
                        alpha=0.7,
                        label='this work', #AGN_cat_names[p_2_catal_AGN],
                        lw=2)  # 2-10  keV')

    # Aird 2015
    plt.plot(x_lx, N_nh20, label='nH<22', ls='dashed', lw=3)
    plt.plot(x_lx, N_nh22, label='22<nH<24', ls='dashed', lw=3)
    plt.plot(x_lx, N_nh24, label='24<nH', ls='dashed', lw=3)
    plt.xlabel(r'$\log_{10}(L^{2-10\, keV}_X/[erg/s])$')
    plt.ylabel(r'$\Phi$=dN/dlogL/dV [1/Mpc$^3$/dex]')
    plt.legend(frameon=False, loc=0, fontsize=10)
    plt.yscale('log')
    plt.xlim((37., 46.5))
    plt.ylim((1 / (2 * volume_mock), 1e-2))
    plt.title(str(np.round(z_min, 2)) + "<z<" + str(np.round(z_max, 2)))
    #plt.grid()
    plt.savefig(fig_out)
    plt.clf()
    print(fig_out, 'written')

    #DATA_XLF = np.transpose([x_lx, nhar, (nhar)*(1-(nharN)**(-0.5)), (nhar)*(1+(nharN)**(-0.5)), N_nh20, N_nh22, N_nh24, phi_h(10**x_lx,z_mean)])
    #np.savetxt(os.path.join(validation_dir_XLF_hard, fig_base +'XLF_soft_'+z_str+'.ascii'), DATA_XLF  )

    fig_out = os.path.join(validation_dir_XLF_hard, LC_dir+'_XLF_ratio_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']) +'.png' )
    # XLF_ratio_
    plt.figure(1, (6, 6))
    plt.axes([0.18, 0.18, 0.75, 0.75])

    for p_2_catal_AGN in p_2_catal_AGNs:
        AGN = AGNs[AGN_cat_names[p_2_catal_AGN]]
        z = AGN['redshift_S']
        AGN = AGN[(z>=z_min)&(z<z_max)]
        z = AGN['redshift_S']
        lx = AGN['LX_hard']
        logNH = AGN['logNH']
        #fx = AGN['FX_soft']
        #lx_0520 = AGN['LX_soft']
        #logm = np.log10(AGN['obs_sm'])
        #lsar = lx - logm
        #mag_r = AGN['SDSS_r_AB']
        nhar = np.histogram(lx, lx_bins)[0] / volume_mock / dlog_lx
        nharN = np.histogram(lx, lx_bins)[0]  # /vol/dlog_lx

        plt.plot(x_lx, nhar / phi_h(10**x_lx, z_mean), label=AGN_cat_names[p_2_catal_AGN])
        plt.fill_between(x_lx,
                        y1=1 - (nharN)**(-0.5),
                        y2=1 + (nharN)**(-0.5),
                        #color='green',
                        alpha=0.3,
                        label=AGN_cat_names[p_2_catal_AGN]+' ERR')
    plt.xlabel(r'$\log_{10}(L^{2-10\, keV}_X/[erg/s])$')
    plt.ylabel(r'mock/model')
    plt.legend(frameon=False, loc=0, fontsize=10)
    plt.xlim((37., 46.5))
    plt.ylim((0.7, 1.3))
    plt.title(str(np.round(z_min, 2)) + "<z<" + str(np.round(z_max, 2)))
    #plt.grid()
    plt.savefig(fig_out)
    plt.clf()
    print(fig_out, 'written')

    if z_mean<0.25:
        fig_out = os.path.join(validation_dir_XLF_soft, LC_dir+'_XLF_soft_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']) +'.png' )
        plt.figure(1, (6, 6))
        plt.axes([0.2, 0.18, 0.75, 0.75])
        #plt.plot(x_lx, phi_h(10**x_lx, z_mean), c='cyan', ls='dashed', lw=2, label='Ai15')  # Aird 2-10 keV LADE')
        h70 = 67.74/70
        H05 = np.loadtxt(os.path.join( validation_dir, 'literature_data', 'hasinger05_z015-02.ascii'), unpack = True)
        x_h05 = ( H05[2] + H05[3] ) * 0.5 - 2*np.log10(h70)
        x_h05_lo = H05[2]
        x_h05_up = H05[3]
        y_H05 = H05[6] * H05[9] * h70**3.
        y_H05_up = (H05[6]+H05[7]) * H05[9] * h70**3.
        y_H05_lo = (H05[6]-H05[8]) * H05[9] * h70**3.
        plt.errorbar(x_h05, y_H05, yerr = [y_H05 - y_H05_lo, y_H05_up-y_H05], xerr=[x_h05-x_h05_lo,x_h05_up-x_h05], label='Ha05', color='k')
        # Model
        ff0 = lambda lxR, g1, g2, A0 : A0 /( (lxR)**g1 + (lxR)**g2 )
        #phiMODHa05z01 = ff0(10**(lx_bins-43.45), 0.35, 2.1 ,3.64*1e-6)
        #plt.plot(lx_bins- 2*np.log10(h70),phiMODHa05z01* h70**3., 'k--' )
        phiMODHa05z01 = ff0(10**(lx_bins-43.45), 0.35, 2.1 ,6*1e-6)
        plt.plot(lx_bins- 2*np.log10(h70),phiMODHa05z01* h70**3., 'k--' )
        #phiMODHa05z01 = ff0(10**(lx_bins-43.45), 0.35, 2.1 ,1e-5)
        #plt.plot(lx_bins- 2*np.log10(h70),phiMODHa05z01* h70**3., 'k--' )
        for p_2_catal_AGN in p_2_catal_AGNs:
            AGN = AGNs[AGN_cat_names[p_2_catal_AGN]]
            z = AGN['redshift_S']
            AGN = AGN[(z>=z_min)&(z<z_max)]
            z = AGN['redshift_S']
            lx = AGN['LX_soft']
            logNH = AGN['logNH']
            #fx = AGN['FX_soft']
            #lx_0520 = AGN['LX_soft']
            #logm = np.log10(AGN['obs_sm'])
            #lsar = lx - logm
            #mag_r = AGN['SDSS_r_AB']

            # Selections for the histogram
            NH20 = (logNH < 22)
            NH22 = (logNH >= 22) & (logNH < 24)
            NH24 = (logNH >= 24)

            N_nh20 = np.histogram(lx[NH20], lx_bins)[0] / volume_mock / dlog_lx
            N_nh22 = np.histogram(lx[NH22], lx_bins)[0] / volume_mock / dlog_lx
            N_nh24 = np.histogram(lx[NH24], lx_bins)[0] / volume_mock / dlog_lx

            nhar = np.histogram(lx, lx_bins)[0] / volume_mock / dlog_lx
            nharN = np.histogram(lx, lx_bins)[0]  # /vol/dlog_lx
            plt.fill_between(x_lx,
                            y1=(nhar) * (1 - (nharN)**(-0.5)),
                            y2=(nhar) * (1 + (nharN)**(-0.5)),
                            #color='g',
                            alpha=0.7,
                            label='AGN, this work',#AGN_cat_names[p_2_catal_AGN],
                            lw=2)  # 2-10  keV')

        plt.plot(x_lx, N_nh20, label='nH<22', ls='dashed', lw=3)
        plt.plot(x_lx, N_nh22, label='22<nH<24', ls='dashed', lw=3)
        plt.plot(x_lx, N_nh24, label='24<nH', ls='dashed', lw=3)
        plt.xlabel(r'$\log_{10}(L^{0.5-2\, keV}_X/[erg/s])$')
        plt.ylabel(r'$\Phi$=dN/dlogL/dV [1/Mpc$^3$/dex]')
        plt.legend(frameon=False, loc=0, fontsize=12)
        plt.yscale('log')
        plt.xlim((37., 46.5))
        plt.ylim((1 / (2 * volume_mock), 1e-2))
        plt.title(str(np.round(z_min, 2)) + "<z<" + str(np.round(z_max, 2)))
        #plt.grid()
        plt.savefig(fig_out)
        plt.clf()
        print(fig_out, 'written')

    if z_mean<0.5:
        path_2_GAMA = os.path.join(os.environ['GIT_STMOD_DATA'], 'data', 'benchmark', 'GAMA', 'forJohan.fits')
        t_gama = Table.read(path_2_GAMA)
        kmag = 8.9 - 2.5*np.log10(t_gama['flux_Kt'])
        zmag = 8.9 - 2.5*np.log10(t_gama['flux_Zt'])
        rmag = 8.9 - 2.5*np.log10(t_gama['flux_rt'])
        imag = 8.9 - 2.5*np.log10(t_gama['flux_it'])
        keep_GAMA = (kmag>0) & (t_gama['Z']>z_min) & (t_gama['Z']<z_max) & (rmag>0)
        t_gama = Table(t_gama[keep_GAMA])
        kmag = 8.9 - 2.5*np.log10(t_gama['flux_Kt'])
        #zmag = 8.9 - 2.5*np.log10(t_gama['flux_Zt'])
        rmag = 8.9 - 2.5*np.log10(t_gama['flux_rt'])
        #imag = 8.9 - 2.5*np.log10(t_gama['flux_it'])

        KCORR_DATA = np.loadtxt( os.path.join( os.environ['GIT_STMOD_DATA'], 'data/models/model_GAL/VISTA_Ks_kcorrections.txt'), unpack = True)
        # redshift distmod bandpass kcorr_median kcorr_16pc kcorr_84pc
        #bandpass_itp = interp1d(KCORR_DATA[0], KCORR_DATA[2])
        kcorr_itp = interp1d(
            np.hstack(( 0., KCORR_DATA[0])),
            np.hstack((KCORR_DATA[3][0], KCORR_DATA[3])) )
        #kmag = absmag + distmod + bandpass + kcorrection
        # absmag = appmag - distmod - kcorr
        # is properly zero-centred; or at least zero-ish; the median +/- NMAD is 0.05 +/- 0.12.

        dm_values = dm_itp(t_gama['Z'].data)
        #kmag_abs = kmag - ( dm_values + bandpass_itp(t_gama['z_peak']) + kcorr_itp(t_gama['z_peak']) )
        kmag_abs = kmag - ( dm_values + kcorr_itp(t_gama['Z'].data) )
        rmag_abs = rmag - ( dm_values + kcorr_itp(t_gama['Z'].data) )

        t_gama['Kmag_abs'] = np.array(kmag_abs)
        t_gama['Rmag_abs'] = np.array(rmag_abs)
        t_gama['dist_mod'] = dm_values
        t_gama['ms'] = np.array(np.log10( t_gama['StellarMass_50'] ))
        volume_GAMA = (cosmo.comoving_volume(z_max).value - cosmo.comoving_volume(z_min).value) * np.pi * 60. / 129600.

    if z_mean>0.25:
        path_2_COSMOS = os.path.join(os.environ['GIT_STMOD_DATA'], 'data', 'benchmark', 'COSMOS','photoz_vers2.0_010312.fits')
        t = Table.read(path_2_COSMOS)
        good = (t['photoz']>z_min )&( t['photoz']< z_max ) & ( t['MK']<0 )&( t['MK']>-40 )&( t['mass_med']<14 )&( t['mass_med']>4 )
        t_cosmos = t[good]
        volume_cosmos = (cosmo.comoving_volume(z_max).value - cosmo.comoving_volume(z_min).value) * np.pi * 2. / 129600.


    # r band LF

    #Schechter_M_z = interp1d(M, phi * 0.7**3.)
    #mrange = np.arange(-24.6, -12.2, 0.01)
    #total_number = np.cumsum(Schechter_M_z(mrange)) * volume * fraction

    plt.figure(1, (6, 6))
    plt.axes([0.2, 0.18, 0.75, 0.75])
    for p_2_catal_AGN in p_2_catal_AGNs:
        AGN = AGNs[AGN_cat_names[p_2_catal_AGN]]
        z = AGN['redshift_S']
        AGN = AGN[(z>=z_min)&(z<z_max)]
        z = AGN['redshift_S']
        #lx = AGN['LX_hard']
        #logNH = AGN['logNH']
        #fx = AGN['FX_soft']
        #lx_0520 = AGN['LX_soft']
        #logm = np.log10(AGN['obs_sm'])
        #lsar = lx - logm
        mag_r = AGN['SDSS_r_AB']

        DM = cosmo.distmod(z).value
        magABS_r = mag_r - DM
        Ragn  = np.histogram(magABS_r, bins_rLF)[0] / volume_mock / dlogR
        RagnN  = np.histogram(magABS_r, bins_rLF)[0]
        plt.fill_between(x_rLF,
                        y1=(Ragn) * (1 - (RagnN)**(-0.5)),
                        y2=(Ragn) * (1 + (RagnN)**(-0.5)),
                        #color='b',
                        alpha=0.7,
                        label='AGN, this work', #AGN_cat_names[p_2_catal_AGN],
                        lw=2)
    # galaxy data
    DM_gal = cosmo.distmod(z_gal).value
    magABS_r_gal = rmag_gal - DM_gal
    Rgal = np.histogram(magABS_r_gal, bins_rLF)[0] / volume_mock / dlogR
    RgalN = np.histogram(magABS_r_gal, bins_rLF)[0]
    plt.fill_between(x_rLF,
                    y1=(Rgal) * (1 - (RgalN)**(-0.5)),
                    y2=(Rgal) * (1 + (RgalN)**(-0.5)),
                    color='k',
                    alpha=0.7,
                    label='UM galaxies',
                    lw=2)
    #if z_mean>0.25:
        #rmagABS_galCOS = t_cosmos['MR']
        #Rgal_cos = np.histogram(rmagABS_galCOS, bins_rLF)[0] / volume_cosmos / dlogR
        #RgalN_cos = np.histogram(rmagABS_galCOS, bins_rLF)[0]
        #plt.fill_between(x_rLF,
                        #y1=(Rgal_cos) * (1 - (RgalN_cos)**(-0.5)),
                        #y2=(Rgal_cos) * (1 + (RgalN_cos)**(-0.5)),
                        #color='g',
                        #alpha=0.7,
                        #label='COSMOS',
                        #lw=2)
    #if z_mean<0.5:
        #ngal, M, phi, Err = np.loadtxt(os.path.join(os.environ['GIT_STMOD_DATA'], 'data', 'validation/validation_GAL/GAMA/loveday_2015', 'LF_loveday_2015.txt'), unpack=True)
        #plt.plot(M, phi, label='loveday 2015 z=0.1')
        #rmagABS_galGAMA = t_gama['Rmag_abs']
        #Rgal_GAMA = np.histogram(rmagABS_galGAMA, bins_rLF)[0] / volume_GAMA / dlogR
        #RgalN_GAMA = np.histogram(rmagABS_galGAMA, bins_rLF)[0]
        #plt.fill_between(x_rLF,
                        #y1=(Rgal_GAMA) * (1 - (RgalN_GAMA)**(-0.5)),
                        #y2=(Rgal_GAMA) * (1 + (RgalN_GAMA)**(-0.5)),
                        #color='orange',
                        #alpha=0.7,
                        #label='GAMA',
                        #lw=2)
    plt.xlabel(r'$M_R$')
    plt.ylabel(r'$\Phi$=dN/dlogL/dV [1/Mpc$^3$/dex]')
    plt.legend(loc=2, fontsize=8)
    plt.yscale('log')
    plt.xlim((-28, -8))
    plt.ylim((1 / (2 * volume_mock), 1))
    plt.title(str(np.round(z_min, 2)) + "<z<" + str(np.round(z_max, 2)))
    #plt.grid()
    fig_out = os.path.join(validation_dir_Rband, LC_dir+'_RLF_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']) +'.png' )
    plt.savefig(fig_out)
    plt.clf()
    print(fig_out, 'written')


    plt.figure(1, (6, 6))
    plt.axes([0.16, 0.15, 0.8, 0.8])

    if z_mean<0.5:
        # Ge17
        x, y_min, y_max = np.loadtxt( os.path.join( agn_data_dir,  'lsar_hist_G17_z025.ascii'), unpack=True)
        fun = interp1d(x, 10**(0.5*(y_min+y_max)) )
        nrm = quad(fun,x.min(), x.max())[0]
        plt.fill_between(x+34, y1=10**(y_min)/nrm,  y2=10**(y_max)/nrm, label='Ge 17, M>8', alpha=0.5,color='grey')
        # A18
        x, y_min, y_max = np.loadtxt( os.path.join( agn_data_dir,  'lsar_hist_A17_010z050_095M100.ascii'), unpack=True)
        fun = interp1d(x, 0.5*(y_min+y_max) )
        nrm = quad(fun,x.min(), x.max())[0]
        plt.fill_between(x+34, y1=(y_min)/nrm,  y2=(y_max)/nrm, label='Ai 18, 9.5<M<10', alpha=0.5,color='green')
        x, y_min, y_max = np.loadtxt( os.path.join( agn_data_dir,  'lsar_hist_A17_010z050_100M105.ascii'), unpack=True)
        fun = interp1d(x, 0.5*(y_min+y_max) )
        nrm = quad(fun,x.min(), x.max())[0]
        plt.fill_between(x+34, y1=(y_min)/nrm,  y2=(y_max)/nrm, label='Ai 18, 10<M<10.5', alpha=0.5,color='red')
        # We17
        #edd_r = 10**np.arange(-6,0,0.1)
        #plt.plot(np.arange(-6,0,0.1)+34, xi_LSAR(edd_r)/dlogf, ls='dashed', lw=2, label='We17 z=0.1' )

    if z_mean>=0.5 and z_mean<1. :
        #G17
        x, y_min, y_max = np.loadtxt( os.path.join( agn_data_dir,  'lsar_hist_G17_z075.ascii'), unpack=True)
        fun = interp1d(x, 10**(0.5*(y_min+y_max)) )
        nrm = quad(fun,x.min(), x.max())[0]
        plt.fill_between(x+34, y1=10**(y_min)/nrm,  y2=10**(y_max)/nrm, alpha=0.5,color='grey')#, label='G17, M>8'
        #A18
        x, y_min, y_max = np.loadtxt( os.path.join( agn_data_dir,  'lsar_hist_A17_050z100_095M100.ascii'), unpack=True)
        fun = interp1d(x, 0.5*(y_min+y_max) )
        nrm = quad(fun,x.min(), x.max())[0]
        plt.fill_between(x+34, y1=(y_min)/nrm,  y2=(y_max)/nrm,alpha=0.5,color='green')# label='A18, 9.5<M<10',
        x, y_min, y_max = np.loadtxt( os.path.join( agn_data_dir,  'lsar_hist_A17_050z100_100M105.ascii'), unpack=True)
        fun = interp1d(x, 0.5*(y_min+y_max) )
        nrm = quad(fun,x.min(), x.max())[0]
        plt.fill_between(x+34, y1=(y_min)/nrm,  y2=(y_max)/nrm,alpha=0.5,color='red')# label='A18, 10<M<10.5',

    if z_mean>=1. and z_mean<1.5 :
        #G17
        x, y_min, y_max = np.loadtxt( os.path.join( agn_data_dir,  'lsar_hist_G17_z125.ascii'), unpack=True)
        fun = interp1d(x, 10**(0.5*(y_min+y_max)) )
        nrm = quad(fun,x.min(), x.max())[0]
        plt.fill_between(x+34, y1=10**(y_min)/nrm,  y2=10**(y_max)/nrm, alpha=0.5,color='grey')# label='G17, M>8',
        x, y_min, y_max = np.loadtxt( os.path.join( agn_data_dir,  'lsar_hist_A17_100z150_095M100.ascii'), unpack=True)
        fun = interp1d(x, 0.5*(y_min+y_max) )
        nrm = quad(fun,x.min(), x.max())[0]
        plt.fill_between(x+34, y1=(y_min)/nrm,  y2=(y_max)/nrm,alpha=0.5,color='green')# label='A18, 9.5<M<10',
        x, y_min, y_max = np.loadtxt( os.path.join( agn_data_dir,  'lsar_hist_A17_100z150_100M105.ascii'), unpack=True)
        fun = interp1d(x, 0.5*(y_min+y_max) )
        nrm = quad(fun,x.min(), x.max())[0]
        plt.fill_between(x+34, y1=(y_min)/nrm,  y2=(y_max)/nrm, alpha=0.5,color='red')#label='A18, 10<M<10.5',

    if z_mean>=1.5 and z_mean<2. :
        x, y_min, y_max = np.loadtxt( os.path.join( agn_data_dir,  'lsar_hist_G17_z175.ascii'), unpack=True)
        fun = interp1d(x, 10**(0.5*(y_min+y_max)) )
        nrm = quad(fun,x.min(), x.max())[0]
        plt.fill_between(x+34, y1=10**(y_min)/nrm,  y2=10**(y_max)/nrm, alpha=0.5,color='grey')#, label='G17, M>8'
        x, y_min, y_max = np.loadtxt( os.path.join( agn_data_dir,  'lsar_hist_A17_150z200_095M100.ascii'), unpack=True)
        fun = interp1d(x, 0.5*(y_min+y_max) )
        nrm = quad(fun,x.min(), x.max())[0]
        plt.fill_between(x+34, y1=(y_min)/nrm,  y2=(y_max)/nrm, alpha=0.5,color='green')#, label='A18, 9.5<M<10'
        x, y_min, y_max = np.loadtxt( os.path.join( agn_data_dir,  'lsar_hist_A17_150z200_100M105.ascii'), unpack=True)
        fun = interp1d(x, 0.5*(y_min+y_max) )
        nrm = quad(fun,x.min(), x.max())[0]
        plt.fill_between(x+34, y1=(y_min)/nrm,  y2=(y_max)/nrm, alpha=0.5,color='red')#, label='A18, 10<M<10.5'

    if z_mean>=2. and z_mean<2.5 :
        x, y_min, y_max = np.loadtxt( os.path.join( agn_data_dir,  'lsar_hist_G17_z225.ascii'), unpack=True)
        fun = interp1d(x, 10**(0.5*(y_min+y_max)) )
        nrm = quad(fun,x.min(), x.max())[0]
        plt.fill_between(x+34, y1=10**(y_min)/nrm,  y2=10**(y_max)/nrm, alpha=0.5,color='grey')#, label='G17, M>8'
        x, y_min, y_max = np.loadtxt( os.path.join( agn_data_dir,  'lsar_hist_A17_200z250_100M105.ascii'), unpack=True)
        fun = interp1d(x, 0.5*(y_min+y_max) )
        nrm = quad(fun,x.min(), x.max())[0]
        plt.fill_between(x+34, y1=(y_min)/nrm,  y2=(y_max)/nrm, alpha=0.5,color='red')#, label='A18, 10<M<10.5'

    if z_mean>=2.5 and z_mean<3. :
        x, y_min, y_max = np.loadtxt( os.path.join( agn_data_dir,  'lsar_hist_G17_z275.ascii'), unpack=True)
        fun = interp1d(x, 10**(0.5*(y_min+y_max)) )
        nrm = quad(fun,x.min(), x.max())[0]
        plt.fill_between(x+34, y1=10**(y_min)/nrm,  y2=10**(y_max)/nrm, alpha=0.5,color='grey')#, label='G17, M>8'
        x, y_min, y_max = np.loadtxt( os.path.join( agn_data_dir,  'lsar_hist_A17_250z300_100M105.ascii'), unpack=True)
        fun = interp1d(x, 0.5*(y_min+y_max) )
        nrm = quad(fun,x.min(), x.max())[0]
        plt.fill_between(x+34, y1=(y_min)/nrm,  y2=(y_max)/nrm, alpha=0.5,color='red')#, label='A18, 10<M<10.5'

    if z_mean>=3 :
        x, y_min, y_max = np.loadtxt( os.path.join( agn_data_dir,  'lsar_hist_G17_z350.ascii'), unpack=True)
        fun = interp1d(x, 10**(0.5*(y_min+y_max)) )
        nrm = quad(fun,x.min(), x.max())[0]
        plt.fill_between(x+34, y1=10**(y_min)/nrm,  y2=10**(y_max)/nrm, alpha=0.5,color='grey')#, label='G17, M>8'
        x, y_min, y_max = np.loadtxt( os.path.join( agn_data_dir,  'lsar_hist_A17_300z400_100M105.ascii'), unpack=True)
        fun = interp1d(x, 0.5*(y_min+y_max) )
        nrm = quad(fun,x.min(), x.max())[0]
        plt.fill_between(x+34, y1=(y_min)/nrm,  y2=(y_max)/nrm,  alpha=0.5,color='red')#label='A18, 10<M<10.5',


    # LSAR histogram
    for p_2_catal_AGN in p_2_catal_AGNs:
        AGN = AGNs[AGN_cat_names[p_2_catal_AGN]]
        z = AGN['redshift_S']
        AGN = AGN[(z>=z_min)&(z<z_max)]
        z = AGN['redshift_S']
        lx = AGN['LX_hard']
        #logNH = AGN['logNH']
        #fx = AGN['FX_soft']
        #lx_0520 = AGN['LX_soft']
        logm = np.log10(AGN['obs_sm'])
        lsar = lx - logm
        #mag_r = AGN['SDSS_r_AB']

        nall = np.histogram(lsar, LSAR_bins)[0] / volume_mock / dlog_LSAR
        nallN = np.histogram(lsar, LSAR_bins)[0]

        zsel = (logm >= 12)
        nall_12 = np.histogram(lsar[zsel], LSAR_bins)[0] / volume_mock / dlog_LSAR
        nallN_12 = np.histogram(lsar[zsel], LSAR_bins)[0]

        zsel = (logm >= 11) & (logm < 12)
        nall_11 = np.histogram(lsar[zsel], LSAR_bins)[0] / volume_mock / dlog_LSAR
        nallN_11 = np.histogram(lsar[zsel], LSAR_bins)[0]

        zsel = (logm >= 10) & (logm < 11)
        nall_10 = np.histogram(lsar[zsel], LSAR_bins)[0] / volume_mock / dlog_LSAR
        nallN_10 = np.histogram(lsar[zsel], LSAR_bins)[0]

        zsel = (logm >= 9) & (logm < 10)
        nall_9 = np.histogram(lsar[zsel], LSAR_bins)[0] / volume_mock / dlog_LSAR
        nallN_9 = np.histogram(lsar[zsel], LSAR_bins)[0]

        zsel = (logm >= 8) & (logm < 9)
        nall_8 = np.histogram(lsar[zsel], LSAR_bins)[0] / volume_mock / dlog_LSAR
        nallN_8 = np.histogram(lsar[zsel], LSAR_bins)[0]

        fun = interp1d(x_LSAR, nall)
        nrm = quad(fun, x_LSAR.min(), x_LSAR.max())[0]
        plt.plot(x_LSAR, nall / nrm, 'k', lw=3 , label=AGN_cat_names[p_2_catal_AGN])

        fun = interp1d(x_LSAR, nall_8)
        nrm = quad(fun, x_LSAR.min(), x_LSAR.max())[0]
        plt.plot(x_LSAR, nall_8 / nrm, 'magenta', lw=2 , label='8<M*<9')

        fun = interp1d(x_LSAR, nall_9)
        nrm = quad(fun, x_LSAR.min(), x_LSAR.max())[0]
        plt.plot(x_LSAR, nall_9 / nrm, 'g', lw=2 , label='9<M*<10')

        fun = interp1d(x_LSAR, nall_10)
        nrm = quad(fun, x_LSAR.min(), x_LSAR.max())[0]
        plt.plot(x_LSAR, nall_10 / nrm, 'r', lw=2 , label='10<M*<11')

        fun = interp1d(x_LSAR, nall_11)
        nrm = quad(fun, x_LSAR.min(), x_LSAR.max())[0]
        plt.plot(x_LSAR, nall_11 / nrm, 'y', lw=2 , label='11<M*<12')

    plt.xlabel(r'$\log_{10}(\lambda_{SAR})$')
    plt.ylabel(r'probability distribution function')
    plt.legend(frameon=False, loc=3)
    plt.yscale('log')
    plt.xlim((30., 35.5))
    plt.ylim((1e-4, 4))
    plt.title('Specific accretion rate, ' + str(np.round(z_min, 2)) + r"<z<" + str(np.round(z_max, 2)))
    #plt.grid()
    fig_out = os.path.join(validation_dir_LSAR, LC_dir+'_LSAR_hist_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']) +'.png' )
    plt.savefig(fig_out)
    plt.clf()
    print(fig_out, 'written')

    # DUTY CYCLE
    # print("duty_cycle_AGN")
    N_gal = np.histogram(logm_gal, bins=bins_SMF)[0]

    plt.figure(2, (6, 6))
    plt.axes([0.16, 0.15, 0.8, 0.8])

    dx = np.log10(0.6777**2)

    if z_mean<=0.35:
        x_41, y_41, y_41_up, y_41_low = np.loadtxt( os.path.join( agn_data_dir, 'duty_cycle_G11_z01_LXhardgt41.ascii'), unpack=True)
        plt.fill_between(x_41+dx, y1=10**(y_41_low),  y2=10**(y_41_up), label=r'G11 z=0.1 $L_X>10^{41}$erg s$^{-1}$', alpha=0.5, color='green')

    if z_mean<=0.35:
        x_42, y_42, y_42_up, y_42_low = np.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_S08_z025_LXhardgt42.ascii'), unpack=True)
        plt.fill_between(x_42+dx, y1=10**(y_42_low),  y2=10**(y_42_up), label=r'S08 z=0.25 $L_X>10^{42}$erg s$^{-1}$', alpha=0.5, color='brown')

    if z_mean<=0.5:
        x_41, y_41, y_41_up, y_41_low = np.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z025_LXhardgt41.ascii'), unpack=True)
        x_42, y_42, y_42_up, y_42_low = np.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z025_LXhardgt42.ascii'), unpack=True)
        x_43, y_43, y_43_up, y_43_low = np.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z025_LXhardgt43.ascii'), unpack=True)
        lab_bib = 'G17 z=0.25'
        #plt.fill_between(x_44, y1=10**(y_44_low),  y2=10**(y_44_up), label=lab_bib+r' $L_X>10^{44}$', alpha=0.5, color='magenta')
        plt.fill_between(x_43+dx, y1=10**(y_43_low),  y2=10**(y_43_up), label=lab_bib+r' $L_X>10^{43}$erg s$^{-1}$', alpha=0.5,color='blue')
        plt.fill_between(x_42+dx, y1=10**(y_42_low),  y2=10**(y_42_up), label=r'$L_X>10^{42}$erg s$^{-1}$', alpha=0.5,color='red')
        plt.fill_between(x_41+dx, y1=10**(y_41_low),  y2=10**(y_41_up), label=r'$L_X>10^{41}$erg s$^{-1}$', alpha=0.5,color='black')

    if z_mean>0.5 and z_mean<1. :
        x_41, y_41, y_41_up, y_41_low = np.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z075_LXhardgt41.ascii'), unpack=True)
        x_42, y_42, y_42_up, y_42_low = np.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z075_LXhardgt42.ascii'), unpack=True)
        x_43, y_43, y_43_up, y_43_low = np.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z075_LXhardgt43.ascii'), unpack=True)
        x_44, y_44, y_44_up, y_44_low = np.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z075_LXhardgt44.ascii'), unpack=True)
        lab_bib = 'G17 z=0.75'
        plt.fill_between(x_44+dx, y1=10**(y_44_low),  y2=10**(y_44_up), label=lab_bib+r' $L_X>10^{44}$erg s$^{-1}$', alpha=0.5, color='magenta')
        plt.fill_between(x_43+dx, y1=10**(y_43_low),  y2=10**(y_43_up), label=r'$L_X>10^{43}$erg s$^{-1}$', alpha=0.5,color='blue')
        plt.fill_between(x_42+dx, y1=10**(y_42_low),  y2=10**(y_42_up), label=r'$L_X>10^{42}$erg s$^{-1}$', alpha=0.5,color='red')
        plt.fill_between(x_41+dx, y1=10**(y_41_low),  y2=10**(y_41_up), label=r'$L_X>10^{41}$erg s$^{-1}$', alpha=0.5,color='black')

    if z_mean>1. and z_mean<1.5 :
        x_42, y_42, y_42_up, y_42_low = np.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z125_LXhardgt42.ascii'), unpack=True)
        x_43, y_43, y_43_up, y_43_low = np.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z125_LXhardgt43.ascii'), unpack=True)
        x_44, y_44, y_44_up, y_44_low = np.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z125_LXhardgt44.ascii'), unpack=True)
        lab_bib = 'G17 z=1.25'
        plt.fill_between(x_44+dx, y1=10**(y_44_low),  y2=10**(y_44_up), label=lab_bib+r' $L_X>10^{44}$erg s$^{-1}$', alpha=0.5, color='magenta')
        plt.fill_between(x_43+dx, y1=10**(y_43_low),  y2=10**(y_43_up), label=r'$L_X>10^{43}$erg s$^{-1}$', alpha=0.5,color='blue')
        plt.fill_between(x_42+dx, y1=10**(y_42_low),  y2=10**(y_42_up), label=r'$L_X>10^{42}$erg s$^{-1}$', alpha=0.5,color='red')
        #plt.fill_between(x_41, y1=10**(y_41_low),  y2=10**(y_41_up),label=r'$L_X>10^{41}$', alpha=0.5,color='black')

    if z_mean>1.5 and z_mean<2. :
        x_42, y_42, y_42_up, y_42_low = np.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z175_LXhardgt42.ascii'), unpack=True)
        x_43, y_43, y_43_up, y_43_low = np.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z175_LXhardgt43.ascii'), unpack=True)
        x_44, y_44, y_44_up, y_44_low = np.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z175_LXhardgt44.ascii'), unpack=True)
        lab_bib = 'G17 z=1.75'
        plt.fill_between(x_44+dx, y1=10**(y_44_low),  y2=10**(y_44_up), label=lab_bib+r' $L_X>10^{44}$erg s$^{-1}$', alpha=0.5, color='magenta')
        plt.fill_between(x_43+dx, y1=10**(y_43_low),  y2=10**(y_43_up), label=r'$L_X>10^{43}$erg s$^{-1}$', alpha=0.5,color='blue')
        plt.fill_between(x_42+dx, y1=10**(y_42_low),  y2=10**(y_42_up), label=r'$L_X>10^{42}$erg s$^{-1}$', alpha=0.5,color='red')
        #plt.fill_between(x_41, y1=10**(y_41_low),  y2=10**(y_41_up),label=r'$L_X>10^{41}$', alpha=0.5,color='black')

    if z_mean>2. and z_mean<2.5 :
        x_42, y_42, y_42_up, y_42_low = np.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z225_LXhardgt42.ascii'), unpack=True)
        x_43, y_43, y_43_up, y_43_low = np.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z225_LXhardgt43.ascii'), unpack=True)
        x_44, y_44, y_44_up, y_44_low = np.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z225_LXhardgt44.ascii'), unpack=True)
        lab_bib = 'G17 z=2.25'
        plt.fill_between(x_44+dx, y1=10**(y_44_low),  y2=10**(y_44_up), label=lab_bib+r' $L_X>10^{44}$erg s$^{-1}$', alpha=0.5, color='magenta')
        plt.fill_between(x_43+dx, y1=10**(y_43_low),  y2=10**(y_43_up), label=r'$L_X>10^{43}$erg s$^{-1}$', alpha=0.5,color='blue')
        plt.fill_between(x_42+dx, y1=10**(y_42_low),  y2=10**(y_42_up), label=r'$L_X>10^{42}$erg s$^{-1}$', alpha=0.5,color='red')
        #plt.fill_between(x_41, y1=10**(y_41_low),  y2=10**(y_41_up),label=r'$L_X>10^{41}$', alpha=0.5,color='black')

    if z_mean>2.5 and z_mean<3. :
        x_42, y_42, y_42_up, y_42_low = np.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z275_LXhardgt42.ascii'), unpack=True)
        x_43, y_43, y_43_up, y_43_low = np.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z275_LXhardgt43.ascii'), unpack=True)
        x_44, y_44, y_44_up, y_44_low = np.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z275_LXhardgt44.ascii'), unpack=True)
        lab_bib = 'G17 z=2.75'
        plt.fill_between(x_44+dx, y1=10**(y_44_low),  y2=10**(y_44_up), label=lab_bib+r' $L_X>10^{44}$erg s$^{-1}$', alpha=0.5, color='magenta')
        plt.fill_between(x_43+dx, y1=10**(y_43_low),  y2=10**(y_43_up), label=r'$L_X>10^{43}$erg s$^{-1}$', alpha=0.5,color='blue')
        plt.fill_between(x_42+dx, y1=10**(y_42_low),  y2=10**(y_42_up), label=r'$L_X>10^{42}$erg s$^{-1}$', alpha=0.5,color='red')
        #plt.fill_between(x_41, y1=10**(y_41_low),  y2=10**(y_41_up),label=r'$L_X>10^{41}$', alpha=0.5,color='black')

    if z_mean>3 :
        lab_bib = 'G17 z=3.5'
        x_43, y_43, y_43_up, y_43_low = np.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z350_LXhardgt43.ascii'), unpack=True)
        x_44, y_44, y_44_up, y_44_low = np.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z350_LXhardgt44.ascii'), unpack=True)
        #x_42, y_42, y_42_up, y_42_low = np.loadtxt( os.path.join( agn_data_dir, 'duty_cycle_H10_z02_LXhardgt42.ascii'), unpack=True)
        plt.fill_between(x_44+dx, y1=10**(y_44_low),  y2=10**(y_44_up), label=lab_bib+r' $L_X>10^{44}$erg s$^{-1}$', alpha=0.5, color='magenta')
        plt.fill_between(x_43+dx, y1=10**(y_43_low),  y2=10**(y_43_up), label=r'$L_X>10^{43}$erg s$^{-1}$', alpha=0.5,color='blue')
        #plt.fill_between(x_42, y1=10**(y_42_low),  y2=10**(y_42_up),label=r'$L_X>10^{42}$', alpha=0.5,color='red')
        #plt.fill_between(x_41, y1=10**(y_41_low),  y2=10**(y_41_up),label=r'$L_X>10^{41}$', alpha=0.5,color='black')

    plt.plot(x_SMF, f_duty(x_SMF), ls='dashed', lw=2, label='input='+str(np.round(f_duty(x_SMF)[0],2)) )

    for p_2_catal_AGN in p_2_catal_AGNs:
        AGN = AGNs[AGN_cat_names[p_2_catal_AGN]]
        z = AGN['redshift_S']
        AGN = AGN[(z>=z_min)&(z<z_max)]
        z = AGN['redshift_S']
        lx = AGN['LX_hard']
        #logNH = AGN['logNH']
        #fx = AGN['FX_soft']
        #lx_0520 = AGN['LX_soft']
        logm = np.log10(AGN['obs_sm'])
        lsar = lx - logm
        #mag_r = AGN['SDSS_r_AB']

        N_agn_a = np.histogram(logm, bins=bins_SMF)[0]
        y = N_agn_a * 1. / N_gal  # *dc_val
        yerr = y * N_agn_a**(-0.5)  # *dc_val
        plt.errorbar(x_SMF, y, yerr=yerr, color='grey')#, label=AGN_cat_names[p_2_catal_AGN])

        tsel = (lx > 41)
        N_agn_41 = np.histogram(logm[tsel], bins=bins_SMF)[0]
        y = N_agn_41 * 1. / N_gal  # *dc_val
        yerr = y * N_agn_41**(-0.5)  # *dc_val
        plt.errorbar(x_SMF, y, yerr=yerr, color='black')#, label=r'L$_X>10^{41}$')

        tsel = (lx > 42)
        N_agn_42 = np.histogram(logm[tsel], bins=bins_SMF)[0]
        y = N_agn_42 * 1. / N_gal  # *dc_val
        yerr = y * N_agn_42**(-0.5)  # *dc_val
        plt.errorbar(x_SMF, y, yerr=yerr, color='red')#, label=r'L$_X>10^{42}$')

        tsel = (lx > 43)
        N_agn_43 = np.histogram(logm[tsel], bins=bins_SMF)[0]
        y = N_agn_43 * 1. / N_gal  # *dc_val
        yerr = y * N_agn_43**(-0.5)  # *dc_val
        plt.errorbar(x_SMF, y, yerr=yerr, color='blue')#, label=r'L$_X>10^{43}$')

        tsel = (lx > 44)
        N_agn_44 = np.histogram(logm[tsel], bins=bins_SMF)[0]
        y = N_agn_44 * 1. / N_gal  # *dc_val
        yerr = y * N_agn_44**(-0.5)  # *dc_val
        plt.errorbar(x_SMF, y, yerr=yerr, color='magenta')#, label=r'L$_X>10^{44}$')

    plt.errorbar(x_SMF-100000, y-100000, yerr=yerr, color='grey', label='this work')#AGN_cat_names[p_2_catal_AGN])
    plt.errorbar(x_SMF-100000, y-100000, yerr=yerr, color='black', label=r'L$_X>10^{41}$')
    plt.errorbar(x_SMF-100000, y-100000, yerr=yerr, color='red', label=r'L$_X>10^{42}$')
    plt.errorbar(x_SMF-100000, y-100000, yerr=yerr, color='blue', label=r'L$_X>10^{43}$')
    plt.errorbar(x_SMF-100000, y-100000, yerr=yerr, color='magenta', label=r'L$_X>10^{44}$')


    plt.xlabel(r'$\log_{10}(M^*/M_\odot)$')
    plt.ylabel(r'$f_{AGN}(M^*, ' + str(np.round(z_min, 2)) + r"<z<" + str(np.round(z_max, 2)) + r')$')
    plt.yscale('log')
    plt.ylim((5e-5, 0.4))
    plt.xlim((8.0, 12.))
    #plt.grid()
    # , '+str(np.round(z_min,2))+"<z<"+str(np.round(z_max,2)))
    plt.title('Duty cycle')
    #plt.legend( loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=2, fancybox=True, title=str(np.round(z_min,2))+"<z<"+str(np.round(z_max,2)) )
    plt.legend(frameon=False, loc=2, ncol=1, fontsize=10)
    #plt.yscale('linear')
    #plt.ylim((0, 0.5))
    fig_out = os.path.join(validation_dir_DC, LC_dir+'_linear_duty_cycle_AGN_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']) +'.png' )
    plt.savefig(fig_out)
    plt.clf()
    print(fig_out, 'written')



    plt.figure(1, (6, 6))
    plt.axes([0.17, 0.15, 0.73, 0.73])
    #for fun, name in zip(smf_ilbert_fun, smf_ilbert_name):
        #plt.plot(mbins, fun(10**mbins)/(0.7**3), label=name, ls='dashed', lw=0.5)
    if z_mean>0.15 and z_mean<0.8:
        d1 = np.loadtxt( os.path.join(os.environ['GIT_STMOD_DATA'], 'data/validation/validation_AGN/literature_data', 'AGN_HGSMF_BO16_03_z_08_LX_430.ascii'), unpack = True , delimiter = ',')
        plt.plot(d1[0][d1[0]<12], 10**d1[1][d1[0]<12], ls='dashed', color='darkgreen', label='Bo16 43, 43.5, 44, 44.5')
        d1 = np.loadtxt( os.path.join(os.environ['GIT_STMOD_DATA'], 'data/validation/validation_AGN/literature_data', 'AGN_HGSMF_BO16_03_z_08_LX_435.ascii'), unpack = True , delimiter = ',')
        plt.plot(d1[0][d1[0]<12], 10**d1[1][d1[0]<12], ls='dashed', color='darkgreen')
        d1 = np.loadtxt( os.path.join(os.environ['GIT_STMOD_DATA'], 'data/validation/validation_AGN/literature_data', 'AGN_HGSMF_BO16_03_z_08_LX_440.ascii'), unpack = True , delimiter = ',')
        plt.plot(d1[0][d1[0]<12], 10**d1[1][d1[0]<12], ls='dashed', color='darkgreen', lw=3)
        d1 = np.loadtxt( os.path.join(os.environ['GIT_STMOD_DATA'], 'data/validation/validation_AGN/literature_data', 'AGN_HGSMF_BO16_03_z_08_LX_445.ascii'), unpack = True , delimiter = ',')
        plt.plot(d1[0][d1[0]<12], 10**d1[1][d1[0]<12], ls='dashed', color='darkgreen', lw=3)

    if z_mean>0.8 and z_mean<1.5:
        d1 = np.loadtxt( os.path.join(os.environ['GIT_STMOD_DATA'], 'data/validation/validation_AGN/literature_data', 'AGN_HGSMF_BO16_08_z_15_LX_430.ascii'), unpack = True , delimiter = ',')
        plt.plot(d1[0][d1[0]<12], 10**d1[1][d1[0]<12], ls='dashed', color='darkgreen', label='Bo16 43, 43.5, 44, 44.5')
        d1 = np.loadtxt( os.path.join(os.environ['GIT_STMOD_DATA'], 'data/validation/validation_AGN/literature_data', 'AGN_HGSMF_BO16_08_z_15_LX_435.ascii'), unpack = True , delimiter = ',')
        plt.plot(d1[0][d1[0]<12], 10**d1[1][d1[0]<12], ls='dashed', color='darkgreen')
        d1 = np.loadtxt( os.path.join(os.environ['GIT_STMOD_DATA'], 'data/validation/validation_AGN/literature_data', 'AGN_HGSMF_BO16_08_z_15_LX_440.ascii'), unpack = True , delimiter = ',')
        plt.plot(d1[0][d1[0]<12], 10**d1[1][d1[0]<12], ls='dashed', color='darkgreen', lw=3)
        d1 = np.loadtxt( os.path.join(os.environ['GIT_STMOD_DATA'], 'data/validation/validation_AGN/literature_data', 'AGN_HGSMF_BO16_08_z_15_LX_445.ascii'), unpack = True , delimiter = ',')
        plt.plot(d1[0][d1[0]<12], 10**d1[1][d1[0]<12], ls='dashed', color='darkgreen', lw=3)

    if z_mean>1.5 and z_mean<2.5:
        d1 = np.loadtxt( os.path.join(os.environ['GIT_STMOD_DATA'], 'data/validation/validation_AGN/literature_data', 'AGN_HGSMF_BO16_15_z_25_LX_430.ascii'), unpack = True , delimiter = ',')
        plt.plot(d1[0][d1[0]<12], 10**d1[1][d1[0]<12], ls='dashed', color='darkgreen', label='Bo16 43, 43.5, 44, 44.5')
        d1 = np.loadtxt( os.path.join(os.environ['GIT_STMOD_DATA'], 'data/validation/validation_AGN/literature_data', 'AGN_HGSMF_BO16_15_z_25_LX_435.ascii'), unpack = True , delimiter = ',')
        plt.plot(d1[0][d1[0]<12], 10**d1[1][d1[0]<12], ls='dashed', color='darkgreen')
        d1 = np.loadtxt( os.path.join(os.environ['GIT_STMOD_DATA'], 'data/validation/validation_AGN/literature_data', 'AGN_HGSMF_BO16_15_z_25_LX_440.ascii'), unpack = True , delimiter = ',')
        plt.plot(d1[0][d1[0]<12], 10**d1[1][d1[0]<12], ls='dashed', color='darkgreen', lw=3)
        d1 = np.loadtxt( os.path.join(os.environ['GIT_STMOD_DATA'], 'data/validation/validation_AGN/literature_data', 'AGN_HGSMF_BO16_15_z_25_LX_445.ascii'), unpack = True , delimiter = ',')
        plt.plot(d1[0][d1[0]<12], 10**d1[1][d1[0]<12], ls='dashed', color='darkgreen', lw=3)

    #if z_mean>0.25:
        #Mgal_COS = np.histogram(t_cosmos['mass_med'], bins=bins_SMF)[0] / volume_cosmos / dlogM
        #Mgal_COSN = np.histogram(t_cosmos['mass_med'], bins=bins_SMF)[0]
        #plt.fill_between(x_SMF,
                        #y1=(Mgal_COS) * (1 - (Mgal_COSN)**(-0.5)),
                        #y2=(Mgal_COS) * (1 + (Mgal_COSN)**(-0.5)),
                        #color='g',
                        #alpha=0.7,
                        #label='COSMOS',
                        #lw=2)
    #if z_mean<0.5:
        #Mgal_GAMA = np.histogram(t_gama['ms'], bins=bins_SMF)[0] / volume_GAMA / dlogM
        #Mgal_GAMAN = np.histogram(t_gama['ms'], bins=bins_SMF)[0]
        #plt.fill_between(x_SMF,
                        #y1=(Mgal_GAMA) * (1 - (Mgal_GAMAN)**(-0.5)),
                        #y2=(Mgal_GAMA) * (1 + (Mgal_GAMAN)**(-0.5)),
                        #color='orange',
                        #alpha=0.7,
                        #label='GAMA',
                        #lw=2)

    plt.plot(x_SMF, N_gal / (volume_mock * dlogM), color='green', label='UM galaxies')

    for p_2_catal_AGN in p_2_catal_AGNs:
        AGN = AGNs[AGN_cat_names[p_2_catal_AGN]]
        z = AGN['redshift_S']
        AGN = AGN[(z>=z_min)&(z<z_max)]
        z = AGN['redshift_S']
        lx = AGN['LX_hard']
        #logNH = AGN['logNH']
        #fx = AGN['FX_soft']
        #lx_0520 = AGN['LX_soft']
        logm = np.log10(AGN['obs_sm'])
        lsar = lx - logm
        #mag_r = AGN['SDSS_r_AB']

        N_agn_410 = np.histogram(logm[(lx > 41.0)], bins=bins_SMF)[0]
        N_agn_420 = np.histogram(logm[(lx > 42.0)], bins=bins_SMF)[0]
        N_agn_430 = np.histogram(logm[(lx > 43.0)], bins=bins_SMF)[0]
        N_agn_435 = np.histogram(logm[(lx > 43.5)], bins=bins_SMF)[0]
        N_agn_440 = np.histogram(logm[(lx > 44.0)], bins=bins_SMF)[0]
        N_agn_445 = np.histogram(logm[(lx > 44.5)], bins=bins_SMF)[0]

        plt.plot(x_SMF, N_agn_a / (volume_mock * dlogM), color='grey', label='AGN, this work') #AGN_cat_names[p_2_catal_AGN])
        plt.plot(x_SMF, N_agn_410 / (volume_mock * dlogM), ls='dashed', label=r'$L_X>41$')
        plt.plot(x_SMF, N_agn_420 / (volume_mock * dlogM), ls='dashed', label=r'$L_X>42$')
        plt.plot(x_SMF, N_agn_430 / (volume_mock * dlogM), ls='dashed', label=r'$L_X>43$')
        #plt.plot(x_SMF, N_agn_435 / (volume_mock * dlogM), ls='dashed', label=r'$L_X>43.5$')
        #plt.plot(x_SMF, N_agn_440 / (volume_mock * dlogM), ls='dashed', label=r'$L_X>44$')
        #plt.plot(x_SMF, N_agn_445 / (volume_mock * dlogM), ls='dashed', label=r'$L_X>44.5$')

    plt.xlabel(r'$\log_{10}(M^*/[M_\odot])$')
    plt.ylabel(r'$\Phi$=dN/dlogL/dV [1/Mpc$^3$/dex]')
    plt.legend(frameon=False, loc=0, fontsize=12, ncol=2)
    plt.yscale('log')
    plt.xlim((8, 12.5))
    plt.ylim((1e-8, 1e-1))
    plt.title('Stellar mass function, ' + str(np.round(z_min, 2)) + "<z<" + str(np.round(z_max, 2)))
    #plt.grid()
    fig_out = os.path.join(validation_dir_Stell, LC_dir+'_SMF_AGN_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']) +'.png' )
    plt.savefig(fig_out)
    plt.clf()
    print(fig_out, 'written')
