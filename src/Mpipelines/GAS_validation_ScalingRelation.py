"""
Validation of the gas model
Figures :
 * scaling relation M500 - LX soft, stellar mass -- LX soft
 * luminosity function (LX soft)
 * logNlogS

"""
import time
t0 = time.time()

import sys, os, glob

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as p

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

from scipy.interpolate import interp1d
from scipy.interpolate import RegularGridInterpolator

from scipy.stats import scoreatpercentile

from astropy.table import Table
import astropy.io.fits as fits

import numpy as np

cs = {}
cs ['ALL_ANY'] = 'black'
cs ['CEN_ANY'] = 'darkgreen'
cs ['PS_CEN_ANY'] = 'red'
cs ['CEN_RS']  = 'firebrick'
cs ['CEN_BC']  = 'navy'
cs ['CEN_RS_hotGAS']  = 'goldenrod'
cs ['CEN_BC_hotGAS']  = 'aqua'
cs ['SAT_ANY'] = 'aquamarine'
cs ['SAT_BC']  = 'cornflowerblue'
cs ['SAT_RS']  = 'lightcoral'

cs ['ALL_Mstar']  = 'deeppink'
cs ['CEN_Mstar']  = 'fuchsia'
cs ['SAT_Mstar']  = 'plum'
cs ['ALL_Mhalo']  = 'forestgreen'
cs ['CEN_Mhalo']  = 'chartreuse'
cs ['SAT_Mhalo']  = 'turquoise'

cs['AGN'] = 'darkmagenta'
cs['XRB'] = 'hotpink'
cs['hotGAS'] = 'darkorange'
cs ['CEN_ANY_MS'] = 'lime'
cs['hotGAS_MS'] = 'moccasin'


cosmoUNIT = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
h = 0.6774
L_box = 2000.0 / h
cosmo = cosmoUNIT

nl = lambda sel : len(sel.nonzero()[0])
zs = np.arange(0.0000001, 7.1, 0.001)
dm_itp = interp1d(zs, cosmo.distmod(zs).value)

#
# GAS MODEL
#
sys.path.append( os.path.join(os.environ['GIT_STMOD'], 'src') )
from models import GAS as GG

z_dir = sys.argv[1]
LC_dir=sys.argv[2] #'FullSky'
C_GAS = GG.GAS(z_dir, b_HS=0.8, logM500c_min=11., logFX_min=-18, LC_dir=LC_dir)

validation_dir       = os.path.join(os.environ['GIT_STMOD'], 'data', 'validation','validation_GAS')
validation_dir_SR = os.path.join(validation_dir, 'ScalingRelation', z_dir)
validation_dir_XLF = os.path.join(validation_dir, 'XsoftLuminosityFunction', z_dir)
#validation_dir_lNlS = os.path.join(validation_dir, 'logNlogS', z_dir)

os.system('mkdir -p ' + validation_dir       )
os.system('mkdir -p ' + validation_dir_SR )
os.system('mkdir -p ' + validation_dir_XLF )
#os.system('mkdir -p ' + validation_dir_lNlS )


# # #
#
#
# Literature
#
#
# # #

SR_dir = os.path.join(os.environ['GIT_STMOD'], 'data/validation/validation_GAS', 'scaling_relations')
itp_z, itp_kt, itp_frac_obs = np.loadtxt( os.path.join( SR_dir, "fraction_05_20_01_24_no_nH.txt"), unpack=True )

nh_vals = 10**np.arange(-2,4+0.01,0.5)
z_vals = np.hstack(( np.arange(0.,0.7,0.05), np.arange(0.8, 4.5, 0.1)))
kT_vals = np.hstack(([0.1, 0.2], np.arange(0.5,8,0.5), [10, 20, 30, 40, 50] ))

XX_nh, YY_z, ZZ_kt = np.meshgrid(nh_vals, z_vals, kT_vals)

shape_i = XX_nh.shape

matrix_z_nh_kt = itp_frac_obs.reshape(shape_i)

attenuation_3d = RegularGridInterpolator((z_vals, np.log10(nh_vals*1e22), kT_vals), matrix_z_nh_kt)

kT_kev, kT_kev_err, Rspec, R500, R500_err, M500, M500_err, Mgas500, Mgas500_err, R2500,  R2500_err, M2500, M2500_err, Mgas2500, Mgas2500_err, tcool, tcool_err, LXxmm, LXxmm_err = np.loadtxt(os.path.join(SR_dir, 'lovisari_2015_table2.ascii'), unpack=True)

redshift_s18, nh, LX_s18, Rktmax, M500NFWFreeze, M500kTextrp, M500NFWHudson, M500NFWAll, M200NFWFreeze, M200kTextrp, M500PlanckSZ = np.loadtxt(os.path.join(SR_dir, 'schllenberger_2018_tableB2B3.ascii'), unpack=True)

s_mi20, l_mi20, b_mi20, T_mi20, LX_mi20, sig_LX_mi20, f_mi20, NH_mi20, Z_mi20 = np.loadtxt(os.path.join(SR_dir, 'migkas_2020_tableC1.ascii'), unpack=True)

B18_id,    B18_z, B18_R500, B18_LXcin, B18_LXcinbol, B18_TXcin, B18_ZXcin, B18_LXcexbol, B18_LXcex, B18_TXcex, B18_ZXcex, B18_MICM, B18_YXcin, B18_M500 = np.loadtxt(os.path.join(SR_dir, 'bulbul_2018_table1_2.ascii'), unpack=True)

Lo20_planckName, Lo20_z, Lo20_M500, Lo20_Mg500, Lo20_kT, Lo20_kTexc, Lo20_LX, Lo20_LXexc, Lo20_Lbol, Lo20_Lbolexc, Lo20_NT, Lo20_fT, Lo20_Nsb, Lo20_fsb = np.loadtxt(os.path.join(SR_dir, 'lovisari_2020_tableA1.ascii'), unpack=True)

XXL_i = fits.open(os.path.join(SR_dir, 'xxl365gc.fits'))[1].data
XXL = XXL_i[(XXL_i['Class']==1)]

WtG = fits.open(os.path.join(SR_dir,'Mantz16_Table2.fits'))[1].data

Zh24ms_MS0, Zh24ms_MS1, Zh24ms_z0, Zh24ms_z1, Zh24ms_LX_T_NM, Zh24ms_LX_T_NM_err, Zh24ms_LX_T_NM_exp, Zh24ms_LX_T_M, Zh24ms_LX_T_M_err, Zh24ms_LX_T_M_exp, Zh24ms_LX_CGM, Zh24ms_LX_CGM_err, Zh24ms_LX_CGM_exp, Zh24ms_LX_CONT, Zh24ms_LX_CONT_err, Zh24ms_LX_CONT_exp = np.loadtxt(os.path.join(SR_dir, 'zhang_2024_StellarMass.ascii'), unpack=True)

Zh24mh_M200m0, Zh24mh_M200m1, Zh24mh_M500c0, Zh24mh_M500c1, Zh24mh_z0, Zh24mh_z1, Zh24mh_LX_T_NM, Zh24mh_LX_T_NM_err, Zh24mh_LX_T_NM_exp, Zh24mh_LX_T_M, Zh24mh_LX_T_M_err, Zh24mh_LX_T_M_exp, Zh24mh_LX_CGM, Zh24mh_LX_CGM_err, Zh24mh_LX_CGM_exp, Zh24mh_LX_CONT, Zh24mh_LX_CONT_err, Zh24mh_LX_CONT_exp = np.loadtxt(os.path.join(SR_dir, 'zhang_2024_HaloMass.ascii'), unpack=True)

efeds = Table.read(os.path.join(SR_dir,'efeds-trimmed.fits'))
sz = ( efeds['z_final']>0.01 ) & ( efeds['z_final']<0.35 )&(efeds['M500']*cosmo.efunc(efeds['z_final'])>4e13)
efeds = efeds[sz]

MERGE_Mh = Table.read(os.path.join(SR_dir, 'FULL_SUMMARY_Mhalobin_centrals_Z3P.fits'))
t_mh = MERGE_Mh[np.unique(MERGE_Mh['file_ID'], return_index=True)[1]]
t_mh = t_mh [(t_mh['N_gal'] > 1000)&(t_mh['M_min'] >= 12.0)&(t_mh['M_max'] <=14.0)]
t = t_mh


def get_mean_scaling_relation(mass_array_log10, mass_proxy_array):#, percents = percents ):
    """
    extracts the scaling relation from simulated data
    """
    #percents = np.arange(0,101,1)
    m_array = np.arange(np.min(mass_array_log10), np.max(mass_array_log10)+1, 0.1)
    mass_mins = m_array[:-1]
    mass_maxs = m_array[1:]
    # define empty columns
    mean_mass_proxy = np.zeros_like(mass_mins)
    std_mass_proxy  = np.zeros_like(mass_mins)
    mean_mass       = np.zeros_like(mass_mins)
    #percentile_values = np.zeros((len(mass_mins), len(percents)))
    # in each bin compute the mean and std of the mass proxy, mass and redshift
    for jj, (m_min, m_max) in enumerate(zip(mass_mins, mass_maxs)):
        selection = ( mass_array_log10 >= m_min ) & ( mass_array_log10 < m_max )
        mproxy_values = mass_proxy_array[ selection ]
        mm_values = mass_array_log10[ selection ]
        mean_mass_proxy[jj] = np.mean( mproxy_values )
        std_mass_proxy [jj] = np.std( mproxy_values )
        #percentile_values[jj] = scoreatpercentile(mproxy_values, percents)
        mean_mass      [jj] = np.mean(mm_values )
    return mass_mins, mean_mass, mass_maxs, mean_mass_proxy, std_mass_proxy

#
#
#


MERGE_Mh = Table.read(os.path.join(SR_dir, 'FULL_SUMMARY_Mhalobin_centrals_Z3P.fits'))
t_mh = MERGE_Mh[np.unique(MERGE_Mh['file_ID'], return_index=True)[1]]
t_mh = t_mh [(t_mh['N_gal'] > 1000)&(t_mh['M_min'] >= 12.0)&(t_mh['M_max'] <=14.0)]

MERGE_Ms = Table.read(os.path.join(SR_dir, 'FULL_SUMMARY_Mstarbin_centrals_Z3P.fits'))
t_ms = MERGE_Ms[np.unique(MERGE_Ms['file_ID'], return_index=True)[1]]
t_ms = t_ms [(t_ms['N_gal'] > 1000)&(t_ms['M_min'] >= 10.75)&(t_ms['M_max'] <=11.75)]

MS_A15, LX_max_A15, LX_min_A15 = np.loadtxt(  os.path.join( os.environ['GIT_EVTXCORR'], 'data', 'anderson2015', 'anderson2015_fig5_total.ascii'), unpack = True)
MS_A15 = np.arange(10.1,11.95, 0.1)+0.05

##
#
# XLF Clusters, literature
#
#
##


dlog10lx = 0.25
lx_bins = np.arange(41.,50.,dlog10lx)
lx_val = ( lx_bins[1:] + lx_bins[:-1] )/2.

# CODEX Finoguenov 2019
f19_LX, f19_n_up, f19_n_low = np.loadtxt(os.path.join(os.environ['GIT_STMOD'], 'data/validation/validation_GAS/luminosity_functions', 'finoguenov_2019.txt'), unpack=True)
L_1e42, dndl, Err_percent = np.loadtxt(os.path.join(os.environ['GIT_STMOD'], 'data/validation/validation_GAS/luminosity_functions', 'pacaud_2015.txt'), unpack=True)
# REFLEX Boehringer 2002
schechter_fun = lambda L, Lstar, alpha, n0 : n0 * np.e**(-L/Lstar) * (L/Lstar)**(-alpha) #/ Lstar
schechter_mod = lambda L, Lstar, alpha, n0 : n0 * np.e**(-L/Lstar) * (L/Lstar)**(-alpha) * (1-(1+L/0.25)**(-1.7))#/ Lstar

boehringer2002 = lambda L : schechter_fun(L, 8e44, 1.3, 8e-8)
boehringer2014 = lambda L : schechter_mod(L, (5*0.7**(-2))*10**44, 2.13, 1.7e-7 * 0.7**3)
# REFLEX SOUTH
boehringer2002s = lambda L : schechter_fun(L, (5*0.5**(-2))*10**44, 1.55, 2.6e-7 * 0.5**3)/np.log(10)
# REFLEX norTH
boehringer2002n = lambda L : schechter_fun(L, (9.4*0.5**(-2))*10**44, 1.79, 0.9e-7 * 0.5**3)/np.log(10)


nBH02S = boehringer2002s(10**lx_val)
nBH02N = boehringer2002n(10**lx_val)
nBH14 = boehringer2014(10**lx_val)
#boehringer2014 = lambda L : schechter_fun(L, 1.08e44*0.7**-2, 1.7, 1.8e-6)
quickFit = boehringer2002(10**lx_val)




#####
# Figures
#

enough_area = (C_GAS.LC_MetaData['area_DC_max']>=0.5*np.max(C_GAS.LC_MetaData['area_DC_max'])) & (C_GAS.LC_MetaData['area_DC_max']>0)
C_GAS.LC_MetaData['mean_area'] = (C_GAS.LC_MetaData['area_DC_min'] + C_GAS.LC_MetaData['area_DC_max']) / 2.
small_difference_minmax_1 = ( C_GAS.LC_MetaData['area_DC_min'] / C_GAS.LC_MetaData['mean_area'] >= 0.8 ) & ( C_GAS.LC_MetaData['area_DC_min'] / C_GAS.LC_MetaData['mean_area'] <= 1.2 )
small_difference_minmax_2 = ( C_GAS.LC_MetaData['area_DC_max'] / C_GAS.LC_MetaData['mean_area'] >= 0.8 ) & ( C_GAS.LC_MetaData['area_DC_max'] / C_GAS.LC_MetaData['mean_area'] <= 1.2 )

for meta in C_GAS.LC_MetaData:#[(enough_area)&(small_difference_minmax_1)&(small_difference_minmax_2)]:
    #
    print(meta)
    # retrieve the resulting catalogues and meta data
    p_2_catalogue = os.path.join(os.environ['UCHUU'], LC_dir, z_dir, 'replication_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']), 'glist.fits')
    p_2_catal_GAS_b06 = os.path.join(os.environ['UCHUU'], LC_dir, z_dir, 'replication_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']), 'Xgas_bHS0.6.fits')
    p_2_catal_GAS_b08 = os.path.join(os.environ['UCHUU'], LC_dir, z_dir, 'replication_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']), 'Xgas_bHS0.8.fits')
    p_2_catal_GAS = p_2_catal_GAS_b08
    p_2_catal_GAS_b10 = os.path.join(os.environ['UCHUU'], LC_dir, z_dir, 'replication_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']), 'Xgas_bHS1.0.fits')
    if os.path.isfile(p_2_catal_GAS)==False :
        continue
    XGA = Table.read(p_2_catal_GAS)
    #XGA_b06 = Table.read(p_2_catal_GAS_b06)
    XGA_b08 = XGA
    #XGA_b10 = Table.read(p_2_catal_GAS_b10)
    z_min, z_max = np.min(XGA['redshift_S']), np.max(XGA['redshift_S'])
    print('z_min, z_max=', z_min, z_max)
    volume_mock = (( cosmo.comoving_volume(z_max) - cosmo.comoving_volume(z_min) ) * meta['mean_area'] * np.pi / 129600.).value
    z_mean = np.mean(XGA['redshift_S'])
    #EZ_mean = cosmo.efunc(z_mean)
    EZ_mean = cosmo.efunc(XGA['redshift_S'])
    xx = np.log10(XGA['M500c'])+np.log10(EZ_mean)
    yy = XGA['CLUSTER_LX_soft_RF_R500c']-np.log10(EZ_mean)
    mass_mins, mean_mass, mass_maxs, mean_mass_proxy, std_mass_proxy = get_mean_scaling_relation(xx, yy)

    ###
    #
    # M500 - LX
    #
    ##

    fig_out = os.path.join(validation_dir_SR, LC_dir+'_M500c-LX-z_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']) +'.png' )

    p.figure(0, (5.5, 5.))

    # popesso 23
    z_po23 = 0.1
    po23_x_m500  = np.array([7.5e12, 2.2e13])*cosmo.efunc(z_po23)
    po23_y_LX500 = np.array([2.8e41, 1.0e42])/cosmo.efunc(z_po23)
    po23_yLO_LX500 = np.array([1e38, 2.2e41])/cosmo.efunc(z_po23)
    po23_yUP_LX500 = np.array([5.5e41, 1.8e42])/cosmo.efunc(z_po23)
    p.errorbar( po23_x_m500 ,
        po23_y_LX500 ,
        yerr= [po23_y_LX500-po23_yLO_LX500 , po23_yUP_LX500-po23_y_LX500],
        marker='v', ls='', mfc='none', label='Po23', markersize=10)

    # Bulbul 18
    p.plot(B18_M500*1e14*cosmo.efunc(B18_z),
        B18_LXcin * 1e44 / cosmo.efunc(B18_z),
        marker='*', ls='', mfc='none', label='Bu19')

    # Mantz 16
    k_correction_3d_Ma16 = attenuation_3d( np.transpose([WtG['Ma16_z'], np.ones_like(WtG['Ma16_z'])*20.1, WtG['Ma16_kT_keV']]))
    p.plot(  WtG['Ma16_Mlen_1e15'] * 1e15 ,
        WtG['Ma16_LX_1e44'] * k_correction_3d_Ma16 * 1e44 / cosmo.efunc(WtG['Ma16_z'] )  ,
        marker='s', ls='', mfc='none', label='Ma16')

    # Liu 22
    p.plot(efeds['M500']*cosmo.efunc(efeds['z_final']), efeds['L500'].T[1]/cosmo.efunc(efeds['z_final']), marker='^', ls='', mfc='none', label='Li22', rasterized = True)

    # Lovsari 2020
    k_correction_3d_Lo20 = attenuation_3d( np.transpose([Lo20_z, np.ones_like(Lo20_z)*20.1, Lo20_kT]))
    p.plot( Lo20_M500*1e14*cosmo.efunc(Lo20_z) ,
        Lo20_LX * k_correction_3d_Lo20 * 1e44 / cosmo.efunc(Lo20_z)   ,
        marker='o', ls='', mfc='none', label='Lo20')

    # Adami 18, XXL
    ok = (XXL['Mgas500kpc']*1e11*10**(1.1) * cosmo.efunc(XXL['z'])>4e13)
    p.plot( XXL['Mgas500kpc'][ok]*1e11*10**(1.1) * cosmo.efunc(XXL['z'][ok]),
        XXL['LXXL500MT'][ok] * 1e42/cosmo.efunc(XXL['z'][ok]) ,
        marker='o', ls='', label='Ad18', mfc='none')

    # Lovisari 2015, groups
    k_correction_3d_Lo15 = attenuation_3d( np.transpose([np.ones_like(kT_kev)*0.02, np.ones_like(kT_kev)*20.1, kT_kev]))
    ok = (M500 * 1e13 / 0.7>4e13)
    p.plot( M500[ok] * 1e13 / 0.7 ,
        LXxmm[ok] * k_correction_3d_Lo15[ok] * 1e43,
        marker='s', ls='', mfc='none', label='Lo15')

    # Schellenberger 18
    k_correction_3d_Sc17 = attenuation_3d( np.transpose([redshift_s18, np.ones_like(redshift_s18)*20.1, np.ones_like(redshift_s18)*2]))
    ok = (M500NFWFreeze*1e14*cosmo.efunc(redshift_s18) >4e13)
    p.plot(M500NFWFreeze[ok]*1e14*cosmo.efunc(redshift_s18[ok]),
        LX_s18[ok] * k_correction_3d_Sc17[ok] * 1e43 / cosmo.efunc(redshift_s18[ok]),
        marker='o', ls='', label='Sc17', mfc='none')

    xx = np.log10(XGA['M500c'])+np.log10(EZ_mean)
    yy = XGA['CLUSTER_LX_soft_RF_R500c']-np.log10(EZ_mean)
    mass_mins, mean_mass, mass_maxs, mean_mass_proxy, std_mass_proxy = get_mean_scaling_relation(xx, yy)
    p.plot(10**mean_mass, 10**(mean_mass_proxy), color='r')
    p.fill_between(10**mean_mass, 10**(mean_mass_proxy-std_mass_proxy), 10**(mean_mass_proxy+std_mass_proxy), alpha=0.2, color='r', label='eqLX')

    #xx = np.log10(XGA['M500c'])+np.log10(EZ_mean)
    #yy = XGA['CLUSTER_LX_soft_RF_R500c_kTcorr']-np.log10(EZ_mean)
    #mass_mins, mean_mass, mass_maxs, mean_mass_proxy, std_mass_proxy = get_mean_scaling_relation(xx, yy)
    #p.plot(10**mean_mass, 10**(mean_mass_proxy), color='g')
    #p.fill_between(10**mean_mass, 10**(mean_mass_proxy-std_mass_proxy), 10**(mean_mass_proxy+std_mass_proxy), alpha=0.2, color='g', label='eqkT')

    xx = np.log10(XGA['M500c'])+np.log10(EZ_mean)
    yy = XGA['CBP_CLUSTER_LX_soft_RF_R500c']-np.log10(EZ_mean)
    mass_mins, mean_mass, mass_maxs, mean_mass_proxy, std_mass_proxy = get_mean_scaling_relation(xx, yy)
    p.plot(10**mean_mass, 10**(mean_mass_proxy), color='b')
    p.fill_between(10**mean_mass, 10**(mean_mass_proxy-std_mass_proxy), 10**(mean_mass_proxy+std_mass_proxy), alpha=0.2, color='b', label='Mock Co20')



    EF = cosmo.efunc((Zh24mh_z0+ Zh24mh_z1)/2.)
    p.errorbar(EF*10**((Zh24mh_M500c0 + Zh24mh_M500c1)/2.), Zh24mh_LX_CGM*10**Zh24mh_LX_CGM_exp / EF,
               xerr = [EF*10**((Zh24mh_M500c0 + Zh24mh_M500c1)/2.)-EF*10**Zh24mh_M500c0, EF*10**Zh24mh_M500c1-EF*10**((Zh24mh_M500c0 + Zh24mh_M500c1)/2.)],
               yerr = Zh24mh_LX_CGM_err*10**Zh24mh_LX_CGM_exp/ EF, color='purple', label='CGM Zhang 24'
    )

    p.errorbar(EF*10**((Zh24mh_M500c0 + Zh24mh_M500c1)/2.), Zh24mh_LX_T_NM*10**Zh24mh_LX_T_NM_exp/ EF,
               xerr = [EF*10**((Zh24mh_M500c0 + Zh24mh_M500c1)/2.)-EF*10**Zh24mh_M500c0, EF*10**Zh24mh_M500c1-EF*10**((Zh24mh_M500c0 + Zh24mh_M500c1)/2.)],
               yerr = Zh24mh_LX_T_NM_err*10**Zh24mh_LX_T_NM_exp/ EF, color='grey', label='CGM+AGN+\n XRB Zhang 24'
    )
    x_val = np.arange(11, 15.5, 0.1)
    #p.plot(10**x_val, 10**(1.7*x_val+19.3), 'k--')#, label='1.7x+19.3\n 1.5x+22')
    p.plot(10**x_val, 10**(1.5*x_val+22), 'k--', label='1.5x+22')

    p.grid()
    p.xlabel(r'$E(z) M_{500c}\; [M_\odot]$')
    p.ylabel(r'$L^{<R_{500c}}_{X,\; 0.5-2 keV}/E(z)$ [erg s$^{-1}$]')
    p.legend(frameon=True, loc=2, ncol=2, fontsize=10, title=z_dir)
    p.xscale('log')
    p.yscale('log')
    p.xlim((1e11, 2e15))
    p.ylim((1e38, 1e47))
    p.tight_layout()
    p.savefig(fig_out)
    p.clf()
    print(fig_out)
    print('-'*30)


    ###
    #
    # stellar mass - LX
    #
    ##

    fig_out = os.path.join(validation_dir_SR, LC_dir+'_StellarMass-LX_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']) +'.png' )

    p.figure(0, (5.5, 5.))
    xx = np.log10(XGA['obs_sm'])
    yy = XGA['CLUSTER_LX_soft_RF_R500c']
    mass_mins, mean_mass, mass_maxs, mean_mass_proxy, std_mass_proxy = get_mean_scaling_relation(xx, yy)
    p.plot(mean_mass, 10**(mean_mass_proxy), color='r')
    p.fill_between(mean_mass, 10**(mean_mass_proxy-std_mass_proxy), 10**(mean_mass_proxy+std_mass_proxy), alpha=0.2, color='r', label='EqLX')

    #xx = np.log10(XGA['obs_sm'])
    #yy = XGA['CLUSTER_LX_soft_RF_R500c_kTcorr']
    #mass_mins, mean_mass, mass_maxs, mean_mass_proxy, std_mass_proxy = get_mean_scaling_relation(xx, yy)
    #p.plot(mean_mass, 10**(mean_mass_proxy), color='g')
    #p.fill_between(mean_mass, 10**(mean_mass_proxy-std_mass_proxy), 10**(mean_mass_proxy+std_mass_proxy), alpha=0.2, color='g', label='EqkT')

    xx = np.log10(XGA['obs_sm'])
    yy = XGA['CBP_CLUSTER_LX_soft_RF_R500c']
    mass_mins, mean_mass, mass_maxs, mean_mass_proxy, std_mass_proxy = get_mean_scaling_relation(xx, yy)
    p.plot(mean_mass, 10**(mean_mass_proxy), color='b')
    p.fill_between(mean_mass, 10**(mean_mass_proxy-std_mass_proxy), 10**(mean_mass_proxy+std_mass_proxy), alpha=0.2, color='b', label='Mock Co20')

    y1 = 10** ( (LX_min_A15+LX_max_A15)/2. )
    yerr_up = 10** LX_max_A15 - y1
    yerr_lo = - 10** LX_min_A15 + y1
    p.errorbar( MS_A15-np.log10(0.7), y1 , xerr=0.05, yerr=  [yerr_lo, yerr_up],  label=r'Anderson 2015', ls='', fmt='s',markersize=8,markerfacecolor='none', color='grey')

    p.errorbar(((Zh24ms_MS0+ Zh24ms_MS1)/2.), Zh24ms_LX_CGM*10**Zh24ms_LX_CGM_exp,
               xerr = [((Zh24ms_MS0+ Zh24ms_MS1)/2.)-Zh24ms_MS0, Zh24ms_MS1-((Zh24ms_MS0+ Zh24ms_MS1)/2.)],
               yerr = Zh24ms_LX_CGM_err*10**Zh24ms_LX_CGM_exp, color='purple', label='CGM Zhang 24'
    )

    p.errorbar(((Zh24ms_MS0+ Zh24ms_MS1)/2.), Zh24ms_LX_T_NM*10**Zh24ms_LX_T_NM_exp,
               xerr = [((Zh24ms_MS0+ Zh24ms_MS1)/2.)-Zh24ms_MS0, Zh24ms_MS1-((Zh24ms_MS0+ Zh24ms_MS1)/2.)],
               yerr = Zh24ms_LX_T_NM_err*10**Zh24ms_LX_T_NM_exp, color='grey', label='CGM+AGN+XRB Zhang 24'
    )

    p.xlabel(r'$\log_{10}(M_*)$ [M$_\odot$]')
    p.ylabel(r'L$^{<R_{500c}}_{X}$ [erg/s]')
    p.legend(fontsize=10, loc='upper left', ncol=2, title=z_dir)#, title='Centrals')
    p.yscale('log')
    p.xlim((10, 12.2))
    p.ylim((1e38, 5e46))
    p.tight_layout()
    p.savefig( fig_out )
    p.clf()
    print(fig_out)
    print('-'*30)

    fig_out = os.path.join(validation_dir_XLF, LC_dir+'_XLF_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']) +'.png' )
    p.figure(0, (5.5, 5.))
    #p.fill_between(np.log10(f19_LX), f19_n_low, f19_n_up, alpha=0.4, label='F19 '+r'$0.1<z<0.2$')
    #p.fill_between(np.log10(L_1e42*1e42*0.697**(-2)), dndl*0.697**(5)*(1-Err_percent/100.), dndl*0.697**(5)*(1+Err_percent/100.), alpha=0.4, label=r'P15 $z<1.2$')
    p.plot(lx_val+np.log10(0.67), nBH02S, 'k--', label='BH02 S')
    p.plot(lx_val+np.log10(0.67), nBH02N, 'g--', label='BH02 N')
    p.plot(lx_val+np.log10(0.67), nBH14, 'b--', label='BH14')
    p.plot(lx_val, quickFit, ls='dotted', label=r'$\alpha=-1.3$')
    #
    #N_GAS = np.histogram( XGA_b06['CLUSTER_LX_soft_RF_R500c'], lx_bins)[0]
    #p.plot(lx_val, N_GAS/dlog10lx/volume_mock/np.log(10), lw=3, label='Mock b=0.6')
    #
    N_GAS = np.histogram( XGA_b08['CLUSTER_LX_soft_RF_R500c'], lx_bins)[0]
    p.plot(lx_val, N_GAS/dlog10lx/volume_mock/np.log(10), lw=3, label='Mock b=0.8')
    #
    #N_GAS = np.histogram( XGA_b10['CLUSTER_LX_soft_RF_R500c'], lx_bins)[0]
    #p.plot(lx_val, N_GAS/dlog10lx/volume_mock/np.log(10), lw=3, label='Mock b=1.0')
    p.xlabel(r'$\log_{10}(L_X/[erg/s])$')
    p.ylabel(r'$dn/d\log_{10}(L_X)/dV$ [Mpc$^{-3}$ dex$^{-1}$]')
    p.legend(frameon=True, loc=3, fontsize=10, title=z_dir)
    #p.xscale('log')
    p.yscale('log')
    p.xlim((38, 45.25))
    p.ylim((2e-9, 2e-2))
    p.tight_layout()
    p.savefig(fig_out)
    p.clf()
    print(fig_out)
    print('-'*30)


    fig_out = os.path.join(validation_dir_SR, LC_dir+'_M500-kT_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']) +'.png' )
    p.figure(0, (5.5, 5.))
    p.plot(np.log10( Lo20_M500*1e14), np.log10(Lo20_kT / 0.77 / cosmo.efunc(Lo20_z)**(2./3.)), marker='o', ls='', mfc='none', label='Lo20')
    p.plot( np.log10(XXL['Mgas500kpc']*1e11*10**(1.1) * cosmo.efunc(XXL['z'])), np.log10(XXL['T300kpc'] /cosmo.efunc(XXL['z'])**(2./3.) ), marker='o', ls='', label='Ad18', mfc='none')
    p.plot(np.log10(M500*1e13/0.7), np.log10(kT_kev), marker='s', ls='', mfc='none', label='Lo15')
    p.plot(np.log10(B18_M500*1e14), np.log10(B18_TXcin/0.77 / cosmo.efunc(B18_z)**(2./3.)), marker='*', ls='', mfc='none', label='Bu18', alpha=0.7)
    p.plot(np.log10(WtG['Ma16_Mlen_1e15']*1e15), np.log10(WtG['Ma16_kT_keV']/0.77/cosmo.efunc(WtG['Ma16_z'])**(2./3.)), marker='s', ls='', mfc='none', label='Ma16', alpha=0.5)
    #
    xx = np.log10(XGA_b08['M500c'])
    yy = np.log10(XGA_b08['CBP_kT'])- 2./3.*np.log10(EZ_mean)
    mass_mins, mean_mass, mass_maxs, mean_mass_proxy, std_mass_proxy = get_mean_scaling_relation(xx, yy)
    p.plot(mean_mass, (mean_mass_proxy), color='b')
    p.fill_between(mean_mass, (mean_mass_proxy-std_mass_proxy), (mean_mass_proxy+std_mass_proxy), alpha=0.2, color='b', label='Mock Co20')

    xx = np.log10(XGA_b08['M500c'])
    yy = XGA_b08['CLUSTER_kT']- 2./3.*np.log10(EZ_mean)
    mass_mins, mean_mass, mass_maxs, mean_mass_proxy, std_mass_proxy = get_mean_scaling_relation(xx, yy)
    p.plot(mean_mass, (mean_mass_proxy), color='r')
    p.fill_between(mean_mass, (mean_mass_proxy-std_mass_proxy), (mean_mass_proxy+std_mass_proxy), alpha=0.2, color='r', label='Eq LX')

    #xx = np.log10(XGA_b08['M500c'])
    #yy = np.log10(XGA_b08['CLUSTER_kT_kTcorr'])- 2./3.*np.log10(EZ_mean)
    #mass_mins, mean_mass, mass_maxs, mean_mass_proxy, std_mass_proxy = get_mean_scaling_relation(xx, yy)
    #p.plot(mean_mass, (mean_mass_proxy), color='g')
    #p.fill_between(mean_mass, (mean_mass_proxy-std_mass_proxy), (mean_mass_proxy+std_mass_proxy), alpha=0.2, color='g', label='Eq kT')

    x_val = np.arange(11, 15, 0.01)
    p.plot(x_val, (0.6*x_val-8.), 'k--', label='0.6x-8')
    #corr=10**(XGA_b08['CLUSTER_LX_soft_RF_R500c']-XGA_b08['CBP_CLUSTER_LX_soft_RF_R500c'])
    #xx = np.log10(XGA_b08['M500c'])
    #yy = np.log10(XGA_b08['CBP_kT']*corr)- 2./3.*np.log10(EZ_mean)
    #mass_mins, mean_mass, mass_maxs, mean_mass_proxy, std_mass_proxy = get_mean_scaling_relation(xx, yy)
    #p.plot(mean_mass, (mean_mass_proxy), color='r')
    #p.fill_between(mean_mass, (mean_mass_proxy-std_mass_proxy), (mean_mass_proxy+std_mass_proxy), alpha=0.2, color='r', label='Mock Se22')
    p.xlabel(r'$\log_{10}(M_{500c}$ [M$_\odot$])')
    p.ylabel(r'$\log_{10}($kT/E$(z)^{2/3}$ [keV])')
    p.legend(frameon=True, loc=4, fontsize=10, title=z_dir)
    #p.xscale('log')
    p.grid()
    #p.yscale('log')
    p.xlim((11.5, 15.5))
    p.ylim((-2, 1.5))
    p.tight_layout()
    p.savefig(fig_out)
    p.clf()
    print(fig_out)

    print('-'*30)
    fig_out = os.path.join(validation_dir_SR, LC_dir+'_LX-kT_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']) +'.png' )
    p.figure(0, (5.5, 5.))
    # Lovsari 2020
    k_correction_3d_Lo20 = attenuation_3d( np.transpose([Lo20_z, np.ones_like(Lo20_z)*20.1, Lo20_kT]))
    p.plot(Lo20_kT*cosmo.efunc(Lo20_z)**(-2./3.)/0.77,  Lo20_LX* k_correction_3d_Lo20 * 1e44 / 1.62 / cosmo.efunc(Lo20_z)**2, marker='o', ls='', mfc='none', label='Lo20')
    # Adami 18, XXL
    p.plot(  XXL['T300kpc']*cosmo.efunc(XXL['z'])**(-2./3.) , XXL['LXXL500MT'] * 1e42 /cosmo.efunc(XXL['z'])**2 , marker='o', ls='', label='Ad18', mfc='none')
    # Lovisari 2015, groups
    k_correction_3d_Lo15 = attenuation_3d( np.transpose([np.ones_like(kT_kev)*0.02, np.ones_like(kT_kev)*20.1, kT_kev]))
    p.plot( kT_kev, LXxmm * k_correction_3d_Lo15 * 1e43, marker='s', ls='', mfc='none', label='Lo15')
    # Bulbul 18
    p.plot(B18_TXcin*cosmo.efunc(B18_z)**(-2./3.)/0.77, B18_LXcin * 1e44/cosmo.efunc(B18_z)**2 , marker='*', ls='', mfc='none', label='Bu19', alpha=0.7)
    # Mantz 16
    k_correction_3d_Ma16 = attenuation_3d( np.transpose([WtG['Ma16_z'], np.ones_like(WtG['Ma16_z'])*20.1, WtG['Ma16_kT_keV']]))
    p.plot(WtG['Ma16_kT_keV']*cosmo.efunc(WtG['Ma16_z'])**(-2./3.)/0.77, WtG['Ma16_LX_1e44'] * k_correction_3d_Ma16 * 1e44/cosmo.efunc(WtG['Ma16_z'])**2 , marker='s', ls='', mfc='none', label='Ma16', alpha=0.5)
    #
    xx = np.log10(XGA_b08['CBP_kT']) -2./3.*np.log10(EZ_mean)
    yy = XGA_b08['CBP_CLUSTER_LX_soft_RF_R500c']-np.log10(EZ_mean)
    mass_mins, mean_mass, mass_maxs, mean_mass_proxy, std_mass_proxy = get_mean_scaling_relation(xx, yy)
    p.plot(10**mean_mass, 10**(mean_mass_proxy), color='b')
    p.fill_between(10**mean_mass, 10**(mean_mass_proxy-std_mass_proxy), 10**(mean_mass_proxy+std_mass_proxy), alpha=0.2, color='b', label='Mock Co20')
    #
    xx = XGA_b08['kT_Mean_oEzm23']
    yy = XGA_b08['LX_Mean_eZm1']
    mass_mins, mean_mass, mass_maxs, mean_mass_proxy, std_mass_proxy = get_mean_scaling_relation(xx, yy)
    p.plot(10**mean_mass, 10**(mean_mass_proxy), color='orange', lw=3, ls='dotted')
    p.fill_between(10**mean_mass, 10**(mean_mass_proxy-std_mass_proxy), 10**(mean_mass_proxy+std_mass_proxy), alpha=0.2, color='orange', label='Eq mean')

    EZ = cosmo.efunc(XGA['redshift_S'])
    print(np.mean(EZ), np.std(EZ), EZ)
    print(np.mean(-2./3.*np.log10(EZ)), np.std(-2./3.*np.log10(EZ)), -2./3.*np.log10(EZ))
    print(np.mean(-2.*np.log10(EZ)), np.std(-2.*np.log10(EZ)), -2.*np.log10(EZ))
    xx = XGA_b08['CLUSTER_kT'] #-2./3.*np.log10(EZ)
    yy = XGA_b08['CLUSTER_LX_soft_RF_R500c'] # -2.*np.log10(EZ)
    mass_mins, mean_mass, mass_maxs, mean_mass_proxy, std_mass_proxy = get_mean_scaling_relation(xx, yy)
    p.plot(10**mean_mass, 10**(mean_mass_proxy), color='r')
    p.fill_between(10**mean_mass, 10**(mean_mass_proxy-std_mass_proxy), 10**(mean_mass_proxy+std_mass_proxy), alpha=0.2, color='r', label='Eq LX')

    #EZ = cosmo.efunc(self.CAT['redshift_S'])
    #self.CAT['kT_Mean_oEzm23'] = ( fun_M_T(np.log10(self.CAT['M500c'])) ) #+ 2./3. * np.log10(EZ)
    #self.CAT['LX_Mean_eZm1'] = fun_M_L( np.log10(self.CAT['M500c']) ) #+ 2*np.log10(EZ)
    ## generate values with correlated scatter
    #corr_scat = np.random.multivariate_normal([0,0], cov_kT_LX, size=len(self.CAT))
    #corr_scat_kT = corr_scat.T[0]
    #corr_scat_LX = corr_scat.T[1]
    #self.CAT['CLUSTER_kT'] = self.CAT['kT_Mean_oEzm23'] + corr_scat_kT + 2./3. * np.log10(EZ)
    #self.CAT['CLUSTER_LX_soft_RF_R500c'] = self.CAT['LX_Mean_eZm1'] + corr_scat_LX + 2*np.log10(EZ)


    x_val = np.arange(-0.7, 1.5, 0.01)
    p.plot(10**x_val, 10**(2.5*x_val+42.), 'k--', label='2.5x+42')
    #
    p.xlabel(r'kT/E$^{2/3}(z)$ [keV]')
    p.ylabel(r'$L^{0.5-2\; keV}_X/E^2(z) \; [erg/s]$')
    p.legend(frameon=True, loc=4, fontsize=10, title=z_dir)
    p.xscale('log')
    p.yscale('log')
    p.xlim((0.05, 40))
    p.ylim((1e39, 5e45))
    p.tight_layout()
    p.savefig(fig_out)
    p.clf()
    print(fig_out)

    ###
    #
    # stellar mass - LX
    #
    ##

    fig_out = os.path.join(validation_dir_SR, LC_dir+'_StellarMass-kT_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']) +'.png' )

    p.figure(0, (5.5, 5.))
    xx = np.log10(XGA['obs_sm'])
    yy = XGA['CLUSTER_kT']
    mass_mins, mean_mass, mass_maxs, mean_mass_proxy, std_mass_proxy = get_mean_scaling_relation(xx, yy)
    p.plot(mean_mass, 10**(mean_mass_proxy), color='r')
    p.fill_between(mean_mass, 10**(mean_mass_proxy-std_mass_proxy), 10**(mean_mass_proxy+std_mass_proxy), alpha=0.2, color='r', label='Eq LX')

    #xx = np.log10(XGA['obs_sm'])
    #yy = np.log10(XGA['CLUSTER_kT_kTcorr'])
    #mass_mins, mean_mass, mass_maxs, mean_mass_proxy, std_mass_proxy = get_mean_scaling_relation(xx, yy)
    #p.plot(mean_mass, 10**(mean_mass_proxy), color='g')
    #p.fill_between(mean_mass, 10**(mean_mass_proxy-std_mass_proxy), 10**(mean_mass_proxy+std_mass_proxy), alpha=0.2, color='g', label='Eq kT')

    xx = np.log10(XGA['obs_sm'])
    yy = np.log10(XGA['CBP_kT'])
    mass_mins, mean_mass, mass_maxs, mean_mass_proxy, std_mass_proxy = get_mean_scaling_relation(xx, yy)
    p.plot(mean_mass, 10**(mean_mass_proxy), color='b')
    p.fill_between(mean_mass, 10**(mean_mass_proxy-std_mass_proxy), 10**(mean_mass_proxy+std_mass_proxy), alpha=0.2, color='b', label='Mock Co20')

    p.xlabel(r'$\log_{10}(M_*)$ [M$_\odot$]')
    p.ylabel('kT [keV]')
    p.legend(fontsize=10, loc='upper left', ncol=2, title=z_dir)#, title='Centrals')
    p.yscale('log')
    p.xlim((10, 12.2))
    p.ylim((0.02, 40))
    p.tight_layout()
    p.savefig( fig_out )
    p.clf()
    print(fig_out, 'written')


    fig_out = os.path.join(validation_dir_SR, LC_dir+'_haloMass500c-StellarMass_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']) +'.png' )

    p.figure(0, (5.5, 5.))
    xx = np.log10(XGA['M500c'])
    yy = np.log10(XGA['obs_sm'])
    mass_mins, mean_mass, mass_maxs, mean_mass_proxy, std_mass_proxy = get_mean_scaling_relation(xx, yy)
    p.plot(10**mean_mass, 10**(mean_mass_proxy), color='k')
    p.fill_between(10**mean_mass, 10**(mean_mass_proxy-std_mass_proxy), 10**(mean_mass_proxy+std_mass_proxy), alpha=0.2, color='k', label='UniverseMachine Uchuu')
    p.plot(10**mass_mins, 10**(2*mass_mins-12.9), 'r--', label='2x-12.9')
    p.plot(10**mass_mins, 10**(0.5*mass_mins+4.6), 'g--', label='0.5x+4.6')

    p.xlabel(r'$M_{500c}$ [M$_\odot$])')
    p.ylabel(r'$M_*$ [M$_\odot$]')
    p.legend(fontsize=10, loc='upper left', ncol=2, title=z_dir)#, title='Centrals')
    p.xscale('log')
    p.yscale('log')
    p.xlim((1e11, 2e15))
    p.ylim((1e9, 2e12))
    p.tight_layout()
    p.savefig( fig_out )
    p.clf()
    print(fig_out, 'written')


    fig_out = os.path.join(validation_dir_SR, LC_dir+'_StellarMass-haloMass500c_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']) +'.png' )

    p.figure(0, (5.5, 5.))
    xx = np.log10(XGA['obs_sm'])
    yy = np.log10(XGA['M500c'])
    mass_mins, mean_mass, mass_maxs, mean_mass_proxy, std_mass_proxy = get_mean_scaling_relation(xx, yy)
    p.plot(10**mean_mass, 10**(mean_mass_proxy), color='k')
    p.fill_between(10**mean_mass, 10**(mean_mass_proxy-std_mass_proxy), 10**(mean_mass_proxy+std_mass_proxy), alpha=0.2, color='k', label='UniverseMachine Uchuu')
    p.plot(10**mass_mins, 10**(2*mass_mins-9.9), 'r--', label='2x-9.9')
    p.plot(10**mass_mins, 10**(0.3*mass_mins+8.4), 'g--', label='0.3x+8.4')

    p.ylabel(r'$M_{500c}$ [M$_\odot$])')
    p.xlabel(r'$M_*$ [M$_\odot$]')
    p.legend(fontsize=10, loc='upper left', ncol=2, title=z_dir)#, title='Centrals')
    p.xscale('log')
    p.yscale('log')
    p.ylim((8e10, 2e15))
    p.xlim((1e9, 2e12))
    p.tight_layout()
    p.savefig( fig_out )
    p.clf()
    print(fig_out, 'written')


