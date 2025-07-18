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

validation_dir       = os.path.join(os.environ['GIT_STMOD_DATA'], 'data', 'validation','validation_GAS')
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

SR_dir = os.path.join(os.environ['GIT_STMOD_DATA'], 'data/validation/validation_GAS', 'scaling_relations')
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

erass1 = Table.read(os.path.join(SR_dir,'erass1cl_main_v2.0.fits'))
sz = ( erass1['BEST_Z']>0.01 ) & ( erass1['BEST_Z']<0.35 )&(erass1['M500']>=0.01)&(erass1['L500']>=0.01)&(erass1['CR500']>0.1)
erass1 = erass1[sz]

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

MS_A15, LX_max_A15, LX_min_A15 = np.loadtxt(  os.path.join( os.environ['GIT_STMOD_DATA'], 'data/benchmark', 'anderson2015', 'anderson2015_fig5_total.ascii'), unpack = True)
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
f19_LX, f19_n_up, f19_n_low = np.loadtxt(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/validation/validation_GAS/luminosity_functions', 'finoguenov_2019.txt'), unpack=True)
L_1e42, dndl, Err_percent = np.loadtxt(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/validation/validation_GAS/luminosity_functions', 'pacaud_2015.txt'), unpack=True)
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
    p_2_catal_GAS_b08 = os.path.join(os.environ['UCHUU'], LC_dir, z_dir, 'replication_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']), 'Xgas_bHS0.8.fits')
    p_2_catal_GAS = p_2_catal_GAS_b08
    if os.path.isfile(p_2_catal_GAS)==False :
        continue
    XGA = Table.read(p_2_catal_GAS)
    z_min, z_max = np.min(XGA['redshift_S']), np.max(XGA['redshift_S'])
    print('z_min, z_max=', z_min, z_max)
    volume_mock = (( cosmo.comoving_volume(z_max) - cosmo.comoving_volume(z_min) ) * meta['mean_area'] * np.pi / 129600.).value
    z_mean = np.mean(XGA['redshift_S'])
    #EZ_mean = cosmo.efunc(z_mean)
    EZ_mean = cosmo.efunc(XGA['redshift_S'])
    xx = np.log10(XGA['M500c'])+np.log10(EZ_mean)
    yy = XGA['CLUSTER_LX_soft_RF_R500c']-np.log10(EZ_mean)
    mass_mins, mean_mass, mass_maxs, mean_mass_proxy, std_mass_proxy = get_mean_scaling_relation(xx, yy)

    BP_allz = np.loadtxt( os.path.join( SR_dir, 'BP-mass-lx.txt'), unpack = True)
    BP_01z03 = np.loadtxt( os.path.join( SR_dir, 'BP-mass-lx-01z03.txt'), unpack = True)

    ###
    #
    # M500 - LX
    #
    ##

    fig_out = os.path.join(validation_dir_SR, LC_dir+'_M500c-LX-z_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']) +'.png' )

    p.figure(0, (5.5, 5.))

    ## popesso 23
    #z_po23 = 0.1
    #po23_x_m500  = np.array([7.5e12, 2.2e13])*cosmo.efunc(z_po23)
    #po23_y_LX500 = np.array([2.8e41, 1.0e42])/cosmo.efunc(z_po23)
    #po23_yLO_LX500 = np.array([1e38, 2.2e41])/cosmo.efunc(z_po23)
    #po23_yUP_LX500 = np.array([5.5e41, 1.8e42])/cosmo.efunc(z_po23)
    #p.errorbar( po23_x_m500 ,
        #po23_y_LX500 ,
        #yerr= [po23_y_LX500-po23_yLO_LX500 , po23_yUP_LX500-po23_y_LX500],
        #marker='+', ls='', mfc='none', color='grey', zorder=1, rasterized=True)

    # Bulbul 18
    p.plot(B18_M500*1e14*cosmo.efunc(B18_z),
        B18_LXcin * 1e44 / cosmo.efunc(B18_z),
        marker='+', ls='', mfc='none',color='grey', zorder=1, rasterized=True)

    # Mantz 16
    k_correction_3d_Ma16 = attenuation_3d( np.transpose([WtG['Ma16_z'], np.ones_like(WtG['Ma16_z'])*20.1, WtG['Ma16_kT_keV']]))
    p.plot(  WtG['Ma16_Mlen_1e15'] * 1e15 ,
        WtG['Ma16_LX_1e44'] * k_correction_3d_Ma16 * 1e44 / cosmo.efunc(WtG['Ma16_z'] )  ,
        marker='+', ls='', mfc='none', color='grey', zorder=1, rasterized=True)

    # Liu 22
    p.plot(efeds['M500']*cosmo.efunc(efeds['z_final']), efeds['L500'].T[1]/cosmo.efunc(efeds['z_final']), marker='+', ls='', mfc='none',color='grey', zorder=1)
    # Bulbul 24
    p.plot(erass1['M500']*cosmo.efunc(erass1['BEST_Z'])*1e13, erass1['L500']*1e42/cosmo.efunc(erass1['BEST_Z']), marker='+', ls='', mfc='none',color='grey', zorder=1)

    # Lovsari 2020
    k_correction_3d_Lo20 = attenuation_3d( np.transpose([Lo20_z, np.ones_like(Lo20_z)*20.1, Lo20_kT]))
    p.plot( Lo20_M500*1e14*cosmo.efunc(Lo20_z) ,
        Lo20_LX * k_correction_3d_Lo20 * 1e44 / cosmo.efunc(Lo20_z)   ,
        marker='+', ls='', mfc='none', color='grey', zorder=1, rasterized=True)

    # Adami 18, XXL
    ok = (XXL['Mgas500kpc']*1e11*10**(1.1) * cosmo.efunc(XXL['z'])>4e13)
    p.plot( XXL['Mgas500kpc'][ok]*1e11*10**(1.1) * cosmo.efunc(XXL['z'][ok]),
        XXL['LXXL500MT'][ok] * 1e42/cosmo.efunc(XXL['z'][ok]) ,
        marker='+', ls='', color='grey', mfc='none', zorder=1, rasterized=True)

    # Lovisari 2015, groups
    k_correction_3d_Lo15 = attenuation_3d( np.transpose([np.ones_like(kT_kev)*0.02, np.ones_like(kT_kev)*20.1, kT_kev]))
    ok = (M500 * 1e13 / 0.7>4e13)
    p.plot( M500[ok] * 1e13 / 0.7 ,
        LXxmm[ok] * k_correction_3d_Lo15[ok] * 1e43,
        marker='+', ls='', mfc='none', color='grey', zorder=1, rasterized=True)

    # Schellenberger 18
    k_correction_3d_Sc17 = attenuation_3d( np.transpose([redshift_s18, np.ones_like(redshift_s18)*20.1, np.ones_like(redshift_s18)*2]))
    ok = (M500NFWFreeze*1e14*cosmo.efunc(redshift_s18) >4e13)
    p.plot(M500NFWFreeze[ok]*1e14*cosmo.efunc(redshift_s18[ok]),
        LX_s18[ok] * k_correction_3d_Sc17[ok] * 1e43 / cosmo.efunc(redshift_s18[ok]),
        marker='+', ls='',color='grey', mfc='none', zorder=1, rasterized=True, label='literature')

    log10_M500c = np.arange(11, 15.5,0.01)
    SR_slope = 1.612
    log10_LX_05_20_mean = np.log10( 10** ( (SR_slope) * (log10_M500c - 15) + 44.7 ) )
    SR_slope_lo = 1.612 - 0.072
    SR_slope_hi = 1.612 + 0.068
    log10_LX_05_20_mean_lo = np.log10( 10** ( (SR_slope_lo) * (log10_M500c - 15) + 44.7 ) )
    log10_LX_05_20_mean_hi = np.log10( 10** ( (SR_slope_hi) * (log10_M500c - 15) + 44.7 ) )
    p.fill_between(10**(log10_M500c), 10**log10_LX_05_20_mean_lo, 10**log10_LX_05_20_mean_hi, alpha=0.2, color='darkblue', label='Comparat 2025, X-corr', zorder=10)#, label='This work, X-corr, M*>'+m0_str, zorder=40 )
    p.plot(10**(log10_M500c), 10**log10_LX_05_20_mean, color='darkblue', ls='--', lw=1, zorder=10)

    logm500c, LX_soft_high, LX_soft_low = np.loadtxt( os.path.join(os.environ['GIT_STMOD_DATA'], 'data/benchmark/popesso2024', 'SR_fig6_left.ascii' ), unpack = True )
    ok=(logm500c>12.6)
    logm500c, LX_soft_high, LX_soft_low = logm500c[ok], LX_soft_high[ok], LX_soft_low[ok]
    p.errorbar(10**logm500c, 10**(0.5*(LX_soft_low+LX_soft_high)),
                xerr=[10**(logm500c)-10**(logm500c-0.25), 10**(logm500c+0.25)-10**(logm500c)],
                yerr=[10**(0.5*(LX_soft_low+LX_soft_high))-10**(LX_soft_low), 10**LX_soft_high-10**(0.5*(LX_soft_low+LX_soft_high))],
                    color='green', ls='none', label='Popesso 2024, IGrM stack', zorder=3)

    Zh24mh_M200m0, Zh24mh_M200m1, Zh24mh_M500c0, Zh24mh_M500c1, Zh24mh_z0, Zh24mh_z1, Zh24mh_LX_T_NM, Zh24mh_LX_T_NM_err, Zh24mh_LX_T_NM_exp, Zh24mh_LX_T_M, Zh24mh_LX_T_M_err, Zh24mh_LX_T_M_exp, Zh24mh_LX_CGM, Zh24mh_LX_CGM_err, Zh24mh_LX_CGM_exp, Zh24mh_LX_CONT, Zh24mh_LX_CONT_err, Zh24mh_LX_CONT_exp = np.loadtxt(os.path.join(SR_dir, 'zhang_2024_HaloMass.ascii'), unpack=True, max_rows=5)

    EF = cosmo.efunc((Zh24mh_z0+ Zh24mh_z1)/2.)
    p.errorbar(EF*10**((Zh24mh_M500c0 + Zh24mh_M500c1)/2.), Zh24mh_LX_CGM*10**Zh24mh_LX_CGM_exp / EF,
               xerr = [EF*10**((Zh24mh_M500c0 + Zh24mh_M500c1)/2.)-EF*10**Zh24mh_M500c0, EF*10**Zh24mh_M500c1-EF*10**((Zh24mh_M500c0 + Zh24mh_M500c1)/2.)],
               yerr = Zh24mh_LX_CGM_err*10**Zh24mh_LX_CGM_exp/ EF, color='purple', label='Zhang 2024, CGM stack', ls='none', zorder=3)

    xx = np.log10(XGA['M500c'])+np.log10(EZ_mean)
    yy = XGA['CLUSTER_LX_soft_RF_R500c']-np.log10(EZ_mean)
    mass_mins, mean_mass, mass_maxs, mean_mass_proxy, std_mass_proxy = get_mean_scaling_relation(xx, yy)
    p.plot(10**mean_mass, 10**(mean_mass_proxy), color='darkred', ls='--', zorder=100)
    p.fill_between(10**mean_mass, 10**(mean_mass_proxy-std_mass_proxy), 10**(mean_mass_proxy+std_mass_proxy), alpha=0.2, color='darkred', label='This work', zorder=100)


    #mass_mins, mean_mass, mass_maxs, mean_mass_proxy, std_mass_proxy = get_mean_scaling_relation(xx, yy)
    #BP_allz = np.loadtxt( os.path.join( SR_dir, 'BP-mass-lx.txt'), unpack = True)
    #BP_01z03 = np.loadtxt( os.path.join( SR_dir, 'BP-mass-lx-01z03.txt'), unpack = True)
    #mass_mins, mean_mass, mass_maxs, mean_mass_proxy, std_mass_proxy = BP_01z03
    #p.plot(10**mean_mass, 10**(mean_mass_proxy), color='darkgreen', ls='--', zorder=100)
    #p.fill_between(10**mean_mass, 10**(mean_mass_proxy-std_mass_proxy), 10**(mean_mass_proxy+std_mass_proxy), alpha=0.2, color='darkgreen', label='BP model 0.1<z<0.3', zorder=100)

    #p.grid()
    p.xlabel(r'$E(z) M_{500c}\; [M_\odot]$')
    p.ylabel(r'$L^{<R_{500c}}_{X,\; 0.5-2 keV}/E(z)$ [erg s$^{-1}$]')
    p.legend(frameon=True, loc=2, ncol=1, fontsize=11)#, title=z_dir)
    p.xscale('log')
    p.yscale('log')
    p.xlim((1e11, 2e15))
    p.ylim((1e38, 1e46))
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
    p.plot(mean_mass, 10**(mean_mass_proxy), color='darkred', ls='--', zorder=100)
    p.fill_between(mean_mass, 10**(mean_mass_proxy-std_mass_proxy), 10**(mean_mass_proxy+std_mass_proxy), alpha=0.2, color='darkred', label='This work', zorder=100)

    y1 = 10** ( (LX_min_A15+LX_max_A15)/2. )
    yerr_up = 10** LX_max_A15 - y1
    yerr_lo = - 10** LX_min_A15 + y1
    p.errorbar( MS_A15-np.log10(0.7), y1 , xerr=0.05, yerr=  [yerr_lo, yerr_up],  label=r'Anderson 2015', ls='', fmt='s',markersize=8,markerfacecolor='none', color='orange', zorder=1)

    p.errorbar(((Zh24ms_MS0+ Zh24ms_MS1)/2.), Zh24ms_LX_CGM*10**Zh24ms_LX_CGM_exp,
               xerr = [((Zh24ms_MS0+ Zh24ms_MS1)/2.)-Zh24ms_MS0, Zh24ms_MS1-((Zh24ms_MS0+ Zh24ms_MS1)/2.)],
               yerr = Zh24ms_LX_CGM_err*10**Zh24ms_LX_CGM_exp, color='purple', label='Zhang 2024', ls='none', zorder=10
    )

    p.xlabel(r'$\log_{10}($M$_*$ [M$_\odot$]$)$')
    p.ylabel(r'$L^{<R_{500c}}_{X,\; 0.5-2 keV}$ [erg s$^{-1}$]')
    p.legend(fontsize=11, loc='upper left', ncol=1)#, title=z_dir)#, title='Centrals')
    p.yscale('log')
    p.xlim((10, 12.2))
    p.ylim((1e38, 1e45))
    p.tight_layout()
    p.savefig( fig_out )
    p.clf()
    print(fig_out)
    print('-'*30)

    fig_out = os.path.join(validation_dir_XLF, LC_dir+'_XLF_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']) +'.png' )
    p.figure(0, (5.5, 5.))
    p.fill_between(np.log10(f19_LX), f19_n_low, f19_n_up, alpha=0.4, label='F19 '+r'$0.1<z<0.2$')
    p.fill_between(np.log10(L_1e42*1e42*0.697**(-2)), dndl*0.697**(5)*(1-Err_percent/100.), dndl*0.697**(5)*(1+Err_percent/100.), alpha=0.4, label=r'P15 $z<1.2$')
    p.plot(lx_val+np.log10(0.67), nBH02S, 'k--', label='BH02 S')
    p.plot(lx_val+np.log10(0.67), nBH02N, 'g--', label='BH02 N')
    p.plot(lx_val+np.log10(0.67), nBH14, 'b--', label='BH14')
    p.plot(lx_val, quickFit, ls='dotted', label=r'$\alpha=-1.3$')
    #
    #N_GAS = np.histogram( XGA_b06['CLUSTER_LX_soft_RF_R500c'], lx_bins)[0]
    #p.plot(lx_val, N_GAS/dlog10lx/volume_mock/np.log(10), lw=3, label='Mock b=0.6')
    #
    N_GAS = np.histogram( XGA['CLUSTER_LX_soft_RF_R500c'], lx_bins)[0]
    p.plot(lx_val, N_GAS/dlog10lx/volume_mock/np.log(10), lw=3, label='Mock')
    #
    #N_GAS = np.histogram( XGA_b10['CLUSTER_LX_soft_RF_R500c'], lx_bins)[0]
    #p.plot(lx_val, N_GAS/dlog10lx/volume_mock/np.log(10), lw=3, label='Mock b=1.0')
    p.xlabel(r'$\log_{10}(L_X/[erg/s])$')
    p.ylabel(r'$dn/d\log_{10}(L_X)/dV$ [Mpc$^{-3}$ dex$^{-1}$]')
    p.legend(frameon=True, loc=3, fontsize=10)#, title=z_dir)
    #p.xscale('log')
    p.yscale('log')
    p.xlim((41.5, 45.5))
    p.ylim((1e-10, 2e-2))
    p.tight_layout()
    p.savefig(fig_out)
    p.clf()
    print(fig_out)
    print('-'*30)

    ez01 = cosmo.efunc(0.1)
    log10_M500c, log10_M500c_up, log10_M500c_lo, kt_lo, kT_BGmod, kt, kt_up = np.transpose([
        [12.97878, 13.08223, 12.88064, 0.22239, 0.43046, 0.62596, 0.79247],
        [13.17507, 13.27851, 13.08488, 0.53541, 0.71267, 0.83320, 0.93197],
        [13.37401, 13.28382, 13.47480, 0.46203, 0.72753, 0.88641, 0.99736],
        [13.57029, 13.47215, 13.67374, 0.67186, 0.94861, 1.64151, 2.26364],
        [13.82759, 13.66578, 14.06631, 1.47621, 2.10900, 2.31084, 3.97526],
        [14.21751, 14.06897, 14.36870, 2.87422, 3.43040, 3.19606, 3.95189],
        [14.65252, 14.36340, 14.86472, 4.07013, 7.08477, 4.71660, 9.97363],])

    fig_out = os.path.join(validation_dir_SR, LC_dir+'_M500-kT_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']) +'.png' )
    p.figure(0, (5.5, 5.))
    p.errorbar(log10_M500c, np.log10(kt/ez01**(2./3.) ), xerr=[log10_M500c_up-log10_M500c, log10_M500c-log10_M500c_lo], yerr=[np.log10(kt_up/ez01**(2./3.))-np.log10(kt/ez01**(2./3.)), np.log10(kt/ez01**(2./3.))-np.log10(kt_lo/ez01**(2./3.))], ls='none', label='Toptun 2025, stack', zorder=20)

    p.plot(np.log10( Lo20_M500*1e14), np.log10(Lo20_kT / 0.77 / cosmo.efunc(Lo20_z)**(2./3.)), marker='+', color='grey', ls='', mfc='none', zorder=1)
    p.plot( np.log10(XXL['Mgas500kpc']*1e11*10**(1.1) * cosmo.efunc(XXL['z'])), np.log10(XXL['T300kpc'] /cosmo.efunc(XXL['z'])**(2./3.) ), marker='+', color='grey', ls='', mfc='none',zorder=1)
    p.plot(np.log10(M500*1e13/0.7), np.log10(kT_kev), marker='+', ls='', color='grey', mfc='none', label='literature',zorder=1)
    p.plot(np.log10(B18_M500*1e14), np.log10(B18_TXcin/0.77 / cosmo.efunc(B18_z)**(2./3.)), marker='+', color='grey', ls='', mfc='none', alpha=0.7,zorder=1)
    p.plot(np.log10(WtG['Ma16_Mlen_1e15']*1e15), np.log10(WtG['Ma16_kT_keV']/0.77/cosmo.efunc(WtG['Ma16_z'])**(2./3.)), marker='+', color='grey', ls='', mfc='none', alpha=0.5,zorder=1)
    #kt_err = (erass1['KT_H']-erass1['KT_L'])/2.
    #ok = (erass1['KT']>0)&(erass1['KT_L']/kt_err>7)
    #p.plot(np.log10(erass1['M500'][ok]*1e13), np.log10(erass1['KT'][ok]/cosmo.efunc(erass1['BEST_Z'][ok])**(2./3.) ), marker='+', ls='', mfc='none',color='orange', zorder=1)

    xx = np.log10(XGA['M500c'])
    yy = np.log10(XGA['CLUSTER_kT'])- 2./3.*np.log10(EZ_mean)
    mass_mins, mean_mass, mass_maxs, mean_mass_proxy, std_mass_proxy = get_mean_scaling_relation(xx, yy)
    p.plot(mean_mass, (mean_mass_proxy), color='darkred', ls='dashed',zorder=10)
    p.fill_between(mean_mass, (mean_mass_proxy-std_mass_proxy), (mean_mass_proxy+std_mass_proxy), alpha=0.4, color='darkred', label='This work',zorder=10)

    #x_val = np.arange(11, 15, 0.01)
    #p.plot(x_val, (0.6*x_val-8.), 'k--', label='0.6x-8')
    M500c_MW = np.log10(9.7e11)
    M500c_MW_up = np.log10(9.7e11+1.5e11)
    M500c_MW_lo = np.log10(9.7e11-1.5e11)
    kt_MW = np.log10( ( 0.153 + 0.178 ) / 2. )
    kt_MW_up = np.log10(0.178)
    kt_MW_lo = np.log10(0.153)

    p.plot([M500c_MW, M500c_MW], [kt_MW_lo, kt_MW_up], color='orange', lw=2, zorder=100)
    p.plot([M500c_MW_lo, M500c_MW_up], [kt_MW, kt_MW], color='orange', lw=2, label='Ponti 2023, MW',zorder=100)

    x_val = np.arange(11, 15, 0.01)
    p.plot(x_val, (0.6*x_val-8.), 'k:', label='0.6x-8')

    p.xlabel(r'$\log_{10}(M_{500c}$ [M$_\odot$])')
    p.ylabel(r'$\log_{10}($kT/E$(z)^{2/3}$ [keV])')
    p.legend(frameon=True, loc=4, fontsize=11)#, title=z_dir)
    #p.xscale('log')
    #p.grid()
    #p.yscale('log')
    p.xlim((11.5, 15.5))
    p.ylim((-1.5, 1.5))
    p.tight_layout()
    p.savefig(fig_out)
    p.clf()
    print(fig_out)


    ez01 = cosmo.efunc(0.1)
    log10_M500c, log10_M500c_up, log10_M500c_lo, kt_lo, kT_BGmod, kt, kt_up = np.transpose([
        [12.97878, 13.08223, 12.88064, 0.22239, 0.43046, 0.62596, 0.79247],
        [13.17507, 13.27851, 13.08488, 0.53541, 0.71267, 0.83320, 0.93197],
        [13.37401, 13.28382, 13.47480, 0.46203, 0.72753, 0.88641, 0.99736],
        [13.57029, 13.47215, 13.67374, 0.67186, 0.94861, 1.64151, 2.26364],
        [13.82759, 13.66578, 14.06631, 1.47621, 2.10900, 2.31084, 3.97526],
        [14.21751, 14.06897, 14.36870, 2.87422, 3.43040, 3.19606, 3.95189],
        [14.65252, 14.36340, 14.86472, 4.07013, 7.08477, 4.71660, 9.97363],])

    fig_out = os.path.join(validation_dir_SR, LC_dir+'_M500-kT_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']) +'-noTT.png' )
    p.figure(0, (5.5, 5.))
    #p.errorbar(log10_M500c, np.log10(kt/ez01**(2./3.) ), xerr=[log10_M500c_up-log10_M500c, log10_M500c-log10_M500c_lo], yerr=[np.log10(kt_up/ez01**(2./3.))-np.log10(kt/ez01**(2./3.)), np.log10(kt/ez01**(2./3.))-np.log10(kt_lo/ez01**(2./3.))], ls='none', label='Toptun 2025, stack', zorder=20)

    p.plot(np.log10( Lo20_M500*1e14), np.log10(Lo20_kT / 0.77 / cosmo.efunc(Lo20_z)**(2./3.)), marker='+', color='grey', ls='', mfc='none', zorder=1)
    p.plot( np.log10(XXL['Mgas500kpc']*1e11*10**(1.1) * cosmo.efunc(XXL['z'])), np.log10(XXL['T300kpc'] /cosmo.efunc(XXL['z'])**(2./3.) ), marker='+', color='grey', ls='', mfc='none',zorder=1)
    p.plot(np.log10(M500*1e13/0.7), np.log10(kT_kev), marker='+', ls='', color='grey', mfc='none', label='literature',zorder=1)
    p.plot(np.log10(B18_M500*1e14), np.log10(B18_TXcin/0.77 / cosmo.efunc(B18_z)**(2./3.)), marker='+', color='grey', ls='', mfc='none', alpha=0.7,zorder=1)
    p.plot(np.log10(WtG['Ma16_Mlen_1e15']*1e15), np.log10(WtG['Ma16_kT_keV']/0.77/cosmo.efunc(WtG['Ma16_z'])**(2./3.)), marker='+', color='grey', ls='', mfc='none', alpha=0.5,zorder=1)
    #kt_err = (erass1['KT_H']-erass1['KT_L'])/2.
    #ok = (erass1['KT']>0)&(erass1['KT_L']/kt_err>7)
    #p.plot(np.log10(erass1['M500'][ok]*1e13), np.log10(erass1['KT'][ok]/cosmo.efunc(erass1['BEST_Z'][ok])**(2./3.) ), marker='+', ls='', mfc='none',color='orange', zorder=1)

    xx = np.log10(XGA['M500c'])
    yy = np.log10(XGA['CLUSTER_kT'])- 2./3.*np.log10(EZ_mean)
    mass_mins, mean_mass, mass_maxs, mean_mass_proxy, std_mass_proxy = get_mean_scaling_relation(xx, yy)
    p.plot(mean_mass, (mean_mass_proxy), color='darkred', ls='dashed',zorder=10)
    p.fill_between(mean_mass, (mean_mass_proxy-std_mass_proxy), (mean_mass_proxy+std_mass_proxy), alpha=0.4, color='darkred', label='This work',zorder=10)

    #x_val = np.arange(11, 15, 0.01)
    #p.plot(x_val, (0.6*x_val-8.), 'k--', label='0.6x-8')
    M500c_MW = np.log10(9.7e11)
    M500c_MW_up = np.log10(9.7e11+1.5e11)
    M500c_MW_lo = np.log10(9.7e11-1.5e11)
    kt_MW = np.log10( ( 0.153 + 0.178 ) / 2. )
    kt_MW_up = np.log10(0.178)
    kt_MW_lo = np.log10(0.153)

    p.plot([M500c_MW, M500c_MW], [kt_MW_lo, kt_MW_up], color='orange', lw=2, zorder=100)
    p.plot([M500c_MW_lo, M500c_MW_up], [kt_MW, kt_MW], color='orange', lw=2, label='Ponti 2023, MW',zorder=100)

    x_val = np.arange(11, 15, 0.01)
    p.plot(x_val, (0.6*x_val-8.), 'k:', label='0.6x-8')

    p.xlabel(r'$\log_{10}(M_{500c}$ [M$_\odot$])')
    p.ylabel(r'$\log_{10}($kT/E$(z)^{2/3}$ [keV])')
    p.legend(frameon=True, loc=4, fontsize=11)#, title=z_dir)
    #p.xscale('log')
    #p.grid()
    #p.yscale('log')
    p.xlim((11.5, 15.5))
    p.ylim((-1.5, 1.5))
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
    EZ = cosmo.efunc(XGA['redshift_S'])
    print(np.mean(EZ), np.std(EZ), EZ)
    print(np.mean(-2./3.*np.log10(EZ)), np.std(-2./3.*np.log10(EZ)), -2./3.*np.log10(EZ))
    print(np.mean(-2.*np.log10(EZ)), np.std(-2.*np.log10(EZ)), -2.*np.log10(EZ))
    xx = np.log10(XGA['CLUSTER_kT']) #-2./3.*np.log10(EZ)
    yy = XGA['CLUSTER_LX_soft_RF_R500c'] # -2.*np.log10(EZ)
    mass_mins, mean_mass, mass_maxs, mean_mass_proxy, std_mass_proxy = get_mean_scaling_relation(xx, yy)
    p.plot(10**mean_mass, 10**(mean_mass_proxy), color='r')
    p.fill_between(10**mean_mass, 10**(mean_mass_proxy-std_mass_proxy), 10**(mean_mass_proxy+std_mass_proxy), alpha=0.2, color='r', label='Mock')

    #x_val = np.arange(-0.7, 1.5, 0.01)
    #p.plot(10**x_val, 10**(2.5*x_val+42.), 'k--', label='2.5x+42')
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
    yy = np.log10(XGA['CLUSTER_kT'])
    mass_mins, mean_mass, mass_maxs, mean_mass_proxy, std_mass_proxy = get_mean_scaling_relation(xx, yy)
    p.plot(mean_mass, 10**(mean_mass_proxy), color='r')
    p.fill_between(mean_mass, 10**(mean_mass_proxy-std_mass_proxy), 10**(mean_mass_proxy+std_mass_proxy), alpha=0.2, color='r', label='Mock')

    #xx = np.log10(XGA['obs_sm'])
    #yy = np.log10(XGA['CLUSTER_kT_kTcorr'])
    #mass_mins, mean_mass, mass_maxs, mean_mass_proxy, std_mass_proxy = get_mean_scaling_relation(xx, yy)
    #p.plot(mean_mass, 10**(mean_mass_proxy), color='g')
    #p.fill_between(mean_mass, 10**(mean_mass_proxy-std_mass_proxy), 10**(mean_mass_proxy+std_mass_proxy), alpha=0.2, color='g', label='Eq kT')

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


