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
LC_dir='FullSky'
C_GAS = GG.GAS(z_dir, b_HS=0.8, logM500c_min=11., logFX_min=-18, LC_dir=LC_dir)

validation_dir       = os.path.join(os.environ['GIT_STMOD_DATA'], 'data', 'validation','validation_GAS')
validation_dir_SR = os.path.join(validation_dir, 'ScalingRelation', z_dir)
validation_dir_XLF = os.path.join(validation_dir, 'XsoftLuminosityFunction', z_dir)
validation_dir_lNlS = os.path.join(validation_dir, 'logNlogS', z_dir)

os.system('mkdir -p ' + validation_dir       )
os.system('mkdir -p ' + validation_dir_SR )
os.system('mkdir -p ' + validation_dir_XLF )
os.system('mkdir -p ' + validation_dir_lNlS )


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


#efeds = Table.read(os.path.join(SR_dir,'efeds-trimmed.fits'))
#sz = ( efeds['z_final']>0.01 ) & ( efeds['z_final']<0.35 )&(efeds['M500']*cosmo.efunc(efeds['z_final'])>4e13)
#efeds = efeds[sz]

MERGE_Mh = Table.read(os.path.join(SR_dir, 'FULL_SUMMARY_Mhalobin_centrals_Z3P.fits'))
t_mh = MERGE_Mh[np.unique(MERGE_Mh['file_ID'], return_index=True)[1]]
t_mh = t_mh [(t_mh['N_gal'] > 1000)&(t_mh['M_min'] >= 12.0)&(t_mh['M_max'] <=14.0)]
t = t_mh


def get_mean_scaling_relation(mass_array_log10, mass_proxy_array):#, percents = percents ):
    """
    extracts the scaling relation from simulated data
    """
    #percents = np.arange(0,101,1)
    m_array = np.arange(int(np.min(mass_array_log10)), int(np.max(mass_array_log10))+1, 0.1)
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



enough_area = (C_GAS.LC_MetaData['area_DC_max']>=0.5*np.max(C_GAS.LC_MetaData['area_DC_max'])) & (C_GAS.LC_MetaData['area_DC_max']>0)
C_GAS.LC_MetaData['mean_area'] = (C_GAS.LC_MetaData['area_DC_min'] + C_GAS.LC_MetaData['area_DC_max']) / 2.
small_difference_minmax_1 = ( C_GAS.LC_MetaData['area_DC_min'] / C_GAS.LC_MetaData['mean_area'] >= 0.8 ) & ( C_GAS.LC_MetaData['area_DC_min'] / C_GAS.LC_MetaData['mean_area'] <= 1.2 )
small_difference_minmax_2 = ( C_GAS.LC_MetaData['area_DC_max'] / C_GAS.LC_MetaData['mean_area'] >= 0.8 ) & ( C_GAS.LC_MetaData['area_DC_max'] / C_GAS.LC_MetaData['mean_area'] <= 1.2 )

for meta in C_GAS.LC_MetaData[(enough_area)&(small_difference_minmax_1)&(small_difference_minmax_2)]:
    #
    print(meta)
    # retrieve the resulting catalogues and meta data
    p_2_catalogue = os.path.join(os.environ['UCHUU'], LC_dir, z_dir, 'replication_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']), 'glist.fits')
    p_2_catal_GAS = os.path.join(os.environ['UCHUU'], LC_dir, z_dir, 'replication_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']), 'Xgas_bHS0.8.fits')
    XGA = Table.read(p_2_catal_GAS)
    z_min, z_max = np.min(XGA['redshift_S']), np.max(XGA['redshift_S'])
    print('z_min, z_max=', z_min, z_max)
    volume_mock = (( cosmo.comoving_volume(z_max) - cosmo.comoving_volume(z_min) ) * meta['mean_area'] * np.pi / 129600.).value
    z_mean = np.mean(XGA['redshift_S'])
    EZ_mean = cosmo.efunc(z_mean)
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
    # Co20
    #p.fill_between(10**mean_mass, mean_mass_proxy-std_mass_proxy*0.5, mean_mass_proxy+std_mass_proxy*0.5, color='g', alpha=0.6, label='Co 20', rasterized = True)
    p.plot(10**mean_mass, 10**(mean_mass_proxy-39), color='b')
    p.fill_between(10**mean_mass, 10**(mean_mass_proxy-std_mass_proxy-39), 10**(mean_mass_proxy+std_mass_proxy-39), alpha=0.2, color='b', label='Mock')

    # popesso 23
    z_po23 = 0.1
    po23_x_m500  = np.array([7.5e12, 2.2e13])*cosmo.efunc(z_po23)
    po23_y_LX500 = np.array([2.8e41, 1.0e42])/1e39/cosmo.efunc(z_po23)
    po23_yLO_LX500 = np.array([1e38, 2.2e41])/1e39/cosmo.efunc(z_po23)
    po23_yUP_LX500 = np.array([5.5e41, 1.8e42])/1e39/cosmo.efunc(z_po23)
    p.errorbar( po23_x_m500 ,
        po23_y_LX500 ,
        yerr= [po23_y_LX500-po23_yLO_LX500 , po23_yUP_LX500-po23_y_LX500],
        marker='v', ls='', mfc='none', label='Po23', markersize=10)

    # Bulbul 18
    p.plot(B18_M500*1e14*cosmo.efunc(B18_z),
        B18_LXcin * 1e44 / cosmo.efunc(B18_z)/1e39,
        marker='*', ls='', mfc='none', label='Bu19')

    # Mantz 16
    k_correction_3d_Ma16 = attenuation_3d( np.transpose([WtG['Ma16_z'], np.ones_like(WtG['Ma16_z'])*20.1, WtG['Ma16_kT_keV']]))
    p.plot(  WtG['Ma16_Mlen_1e15'] * 1e15 ,
        WtG['Ma16_LX_1e44'] * k_correction_3d_Ma16 * 1e44 / cosmo.efunc(WtG['Ma16_z'] ) /1e39 ,
        marker='s', ls='', mfc='none', label='Ma16')

    # Liu 22
    #p.plot(efeds['M500']*cosmo.efunc(efeds['z_final']), efeds['L500'].T[1]/1e39/cosmo.efunc(efeds['z_final']), marker='^', ls='', mfc='none', label='Li22', rasterized = True)
    #p.plot(efeds['M500']*cosmo.efunc(efeds['z_final']), efeds['L500_Bulbul2019'].T[1]/1e39/cosmo.efunc(efeds['z_final']), marker='^', ls='', mfc='none', label='Li22 B', rasterized = True)
    #p.plot(efeds['M500']*cosmo.efunc(efeds['z_final']), efeds['L500_Sereno2020'].T[1]/1e39/cosmo.efunc(efeds['z_final']), marker='^', ls='', mfc='none', label='Li22 S', rasterized = True)
    #p.plot(efeds['M500']*cosmo.efunc(efeds['z_final']), efeds['L500_Lovisari2015'].T[1]/1e39/cosmo.efunc(efeds['z_final']), marker='^', ls='', mfc='none', label='Li22 L', rasterized = True)

    # Lovsari 2020
    k_correction_3d_Lo20 = attenuation_3d( np.transpose([Lo20_z, np.ones_like(Lo20_z)*20.1, Lo20_kT]))
    p.plot( Lo20_M500*1e14*cosmo.efunc(Lo20_z) ,
        Lo20_LX * k_correction_3d_Lo20 * 1e44 / cosmo.efunc(Lo20_z) /1e39  ,
        marker='o', ls='', mfc='none', label='Lo20')

    # Adami 18, XXL
    ok = (XXL['Mgas500kpc']*1e11*10**(1.1) * cosmo.efunc(XXL['z'])>4e13)
    p.plot( XXL['Mgas500kpc'][ok]*1e11*10**(1.1) * cosmo.efunc(XXL['z'][ok]),
        XXL['LXXL500MT'][ok] * 1e42/cosmo.efunc(XXL['z'][ok]) /1e39,
        marker='o', ls='', label='Ad18', mfc='none')

    # Lovisari 2015, groups
    k_correction_3d_Lo15 = attenuation_3d( np.transpose([np.ones_like(kT_kev)*0.02, np.ones_like(kT_kev)*20.1, kT_kev]))
    ok = (M500 * 1e13 / 0.7>4e13)
    p.plot( M500[ok] * 1e13 / 0.7 ,
        LXxmm[ok] * k_correction_3d_Lo15[ok] * 1e43/1e39,
        marker='s', ls='', mfc='none', label='Lo15')

    # Schellenberger 18
    k_correction_3d_Sc17 = attenuation_3d( np.transpose([redshift_s18, np.ones_like(redshift_s18)*20.1, np.ones_like(redshift_s18)*2]))
    ok = (M500NFWFreeze*1e14*cosmo.efunc(redshift_s18) >4e13)
    p.plot(M500NFWFreeze[ok]*1e14*cosmo.efunc(redshift_s18[ok]),
        LX_s18[ok] * k_correction_3d_Sc17[ok] * 1e43 / cosmo.efunc(redshift_s18[ok])/1e39,
        marker='o', ls='', label='Sc17', mfc='none')


    p.errorbar( t['M500c_mean']*cosmo.efunc(t['redshift_mean']), t['total_LX_R500c']/1e39/cosmo.efunc(t['redshift_mean']),
                yerr = [t['total_LX_R500c' ]/cosmo.efunc(t['redshift_mean'])/1e39-t['total_LX_R500c_lo']/cosmo.efunc(t['redshift_mean'])/1e39, t['total_LX_R500c_up']/cosmo.efunc(t['redshift_mean'])/1e39-t['total_LX_R500c']/cosmo.efunc(t['redshift_mean'])/1e39],
                xerr = t['M500c_std']*cosmo.efunc(t['redshift_mean']),
                lw=2, ls='', color='darkgreen', label='Total emission Zh23' )
    # Hot gas
    #p.errorbar( t['M500c_mean']*cosmo.efunc(t['redshift_mean']), t['hotgas_LX_R500c']/cosmo.efunc(t['redshift_mean'])/1e39,
                #yerr = [t['hotgas_LX_R500c' ]/cosmo.efunc(t['redshift_mean'])/1e39-t['hotgas_LX_R500c_lo']/cosmo.efunc(t['redshift_mean'])/1e39, t['hotgas_LX_R500c_up']/cosmo.efunc(t['redshift_mean'])/1e39-t['hotgas_LX_R500c']/cosmo.efunc(t['redshift_mean'])/1e39],
                #xerr = [ t['M500c_mean']*cosmo.efunc(t['redshift_mean'])-t['M500c_Q05']*cosmo.efunc(t['redshift_mean']),t['M500c_Q95']*cosmo.efunc(t['redshift_mean'])-t['M500c_mean']*cosmo.efunc(t['redshift_mean'])],
                #lw=1, ls='', color='orange' )
    p.errorbar( t['M500c_mean']*cosmo.efunc(t['redshift_mean']), t['hotgas_LX_R500c']/cosmo.efunc(t['redshift_mean'])/1e39,
                yerr = [t['hotgas_LX_R500c' ]/cosmo.efunc(t['redshift_mean'])/1e39-t['hotgas_LX_R500c_lo']/cosmo.efunc(t['redshift_mean'])/1e39, t['hotgas_LX_R500c_up']/cosmo.efunc(t['redshift_mean'])/1e39-t['hotgas_LX_R500c']/cosmo.efunc(t['redshift_mean'])/1e39],
                xerr = t['M500c_std']*cosmo.efunc(t['redshift_mean']),
                lw=2, ls='', color='orange', label='Hot gas Zh23' )


    #fun_ch22 = lambda M500 : np.log(3.36) + np.log(1e43) + (1.44 - 0.07*np.log(1.35/1.35) )*np.log(M500/1.4e14) + 2*np.log(cosmo.efunc(0.35)/cosmo.efunc(0.35))-0.51**np.log(1.35/1.35)
    #fun_ch22 = lambda M500 : np.log(3.36) + np.log(1e43) + 1.44*np.log(M500/1.4e14)
    #M500s = 10**np.arange(13,15,0.01)*cosmo.efunc(0.35)
    #lx = np.e**(fun_ch22(M500s)) /1e39/cosmo.efunc(0.35)
    #p.plot(M500s, lx, ls='dashed',color='k')
    slope = (44-np.log10(3e41))/(np.log10(4e14)-13)
    oo = np.log10(3e41)-slope*13
    fun = lambda x : slope * x + oo
    x_1 = np.arange(13.2,15,0.01)
    p.plot(10**x_1, 10**(fun(x_1)-39), ls='dashed',color='k', label='slope 1.57')#r'$\log_{10}(y)=1.57\log_{10}(x)+21$')
    #p.plot([1e13, 4e14], [3e41/1e39,1e44/1e39], ls='dashed',color='k')

    p.plot([5.5e11, 5.5e11]    , [1, 10], 'k--', alpha=0.8)
    p.plot([1.5e12, 1.5e12], [1, 10], 'k--', alpha=0.8)
    p.plot([5e12, 5e12]    , [1, 10], 'k--', alpha=0.8)
    p.plot([4e13, 4e13]    , [1, 10], 'k--', alpha=0.8)
    p.text( 6e11       , 1.5, 'MW' , alpha=0.8)
    p.text( 2e12       , 1.5, 'M31', alpha=0.8)
    p.text( 8e12       , 1.5, 'Groups', alpha=0.8)
    p.text( 8e13       , 1.5, 'Clusters', alpha=0.8)

    p.xlabel(r'$E(z) M_{500c}\; [M_\odot]$')
    p.ylabel(r'$L^{<R_{500c}}_X/E(z)$ [$10^{39}$ erg s$^{-1}$]')
    p.legend(frameon=True, loc=2, ncol=2, fontsize=10)
    p.xscale('log')
    p.yscale('log')
    p.xlim((5e11, 2e15))
    p.ylim((1, 3e6))
    #p.title(title_str)
    p.tight_layout()
    #p.grid()
    p.savefig(fig_out)
    p.clf()
    print(fig_out, 'written')


    ###
    #
    # stellar mass - LX
    #
    ##
    xx = np.log10(XGA['obs_sm'])
    yy = 10**XGA['CLUSTER_LX_soft_RF_R500c']
    mass_mins, mean_mass, mass_maxs, mean_mass_proxy, std_mass_proxy = get_mean_scaling_relation(xx, yy)

    fig_out = os.path.join(validation_dir_SR, LC_dir+'_StellarMass-LX_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']) +'.png' )

    p.figure(0, (5.5, 5.))
    # Mock
    #p.plot(xx, 10**yy, 'k,')
    p.plot(mean_mass, 10**(mean_mass_proxy), color='b')
    p.fill_between(mean_mass, 10**(mean_mass_proxy-std_mass_proxy), 10**(mean_mass_proxy+std_mass_proxy), alpha=0.2, color='b', label='Mock')

    y1 = 10** ( (LX_min_A15+LX_max_A15)/2. )
    yerr_up = 10** LX_max_A15 - y1
    yerr_lo = - 10** LX_min_A15 + y1
    p.errorbar( MS_A15-np.log10(0.7), y1 , xerr=0.05, yerr=  [yerr_lo, yerr_up],  label=r'Anderson 2015', ls='', fmt='s',markersize=8,markerfacecolor='none', color='grey')

    p.errorbar( t_ms['GAL_stellarMass_mean'],
                t_ms['total_LX_R500c'],
                yerr = [ t_ms['total_LX_R500c'] - t_ms['total_LX_R500c_lo'], t_ms['total_LX_R500c_up'] - t_ms['total_LX_R500c'] ],
                xerr = t_ms['GAL_stellarMass_std'],
                lw=2, ls='', color=cs ['CEN_ANY_MS'], label='M*, total' )

    p.errorbar( t_ms['GAL_stellarMass_mean'],
                t_ms['hotgas_LX_R500c'],
                yerr = [ t_ms['hotgas_LX_R500c'] - t_ms['hotgas_LX_R500c_lo'], t_ms['hotgas_LX_R500c_up'] - t_ms['hotgas_LX_R500c'] ],
                xerr = t_ms['GAL_stellarMass_std'],
                lw=2, ls='', color=cs['hotGAS_MS'], label='M*, hot gas' )

    p.errorbar( t_mh['GAL_stellarMass_mean'],
                t_mh['total_LX_R500c'],
                yerr = [ t_mh['total_LX_R500c'] - t_mh['total_LX_R500c_lo'], t_mh['total_LX_R500c_up'] - t_mh['total_LX_R500c'] ],
                xerr = t_mh['GAL_stellarMass_std'],
                lw=2, ls='', color=cs ['CEN_ANY'], label='Mh, total' )

    p.errorbar( t_mh['GAL_stellarMass_mean'],
                t_mh['hotgas_LX_R500c'],
                yerr = [ t_mh['hotgas_LX_R500c'] - t_mh['hotgas_LX_R500c_lo'], t_mh['hotgas_LX_R500c_up'] - t_mh['hotgas_LX_R500c'] ],
                xerr = t_mh['GAL_stellarMass_std'],
                lw=2, ls='', color=cs['hotGAS'], label='Mh, hot gas' )

    #p.plot([10.75, 10.75], [0.9e39, 3e39], 'k--', alpha=0.8)
    #p.plot([11., 11]     , [0.9e39, 3e39], 'k--', alpha=0.8)
    #p.plot([11.25, 11.25], [0.9e39, 3e39], 'k--', alpha=0.8)
    #p.plot([11.75, 11.75], [0.9e39, 3e39], 'k--', alpha=0.8)
    #p.text(10.8,  1e39, 'MW' , alpha=0.8)
    #p.text(11.02, 1e39, 'M31', alpha=0.8)
    #p.text(11.35, 1e39, 'Groups', alpha=0.8)
    p.xlabel(r'$\log_{10}(M_*)$ [M$_\odot$]')
    p.ylabel(r'L$^{<R_{500c}}_{X}$ [erg/s]')
    p.legend(fontsize=12, loc='upper left', ncol=1)#, title='Centrals')
    p.yscale('log')
    p.xlim((10, 12.2))
    p.ylim((1e38, 3e45))
    p.tight_layout()
    p.savefig( fig_out )
    p.clf()
    print(fig_out, 'written')



