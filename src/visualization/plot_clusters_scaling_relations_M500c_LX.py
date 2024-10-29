"""
python plot_clusters_scaling_relations_M500c_LX.py MD10
python plot_clusters_scaling_relations_M500c_LX.py MD04
python plot_clusters_scaling_relations_M500c_LX.py MD40
python plot_clusters_scaling_relations_M500c_LX.py UNIT_fA1i_DIR
python plot_clusters_scaling_relations_M500c_LX.py UNIT_fA1_DIR
python plot_clusters_scaling_relations_M500c_LX.py UNIT_fA2_DIR

"""
sys.path.append( os.path.join(os.environ['GIT_STMOD'], 'src') )
from io.Profile import *

from scipy.stats import norm
from astropy.table import Table, Column, hstack, vstack

import sys
import os
from scipy.interpolate import interp1d
import astropy.io.fits as fits
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
from scipy.stats import norm
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
from scipy.stats import scoreatpercentile
import matplotlib.pyplot as p

import numpy as n
#import h5py
import time
print('Plots clusters catalog, scaling relATIONS')
t0 = time.time()

#hdu_0 = get_cat(env='MD04')
#hdu_1 = get_cat(env='MD40')
#env = 'MD10'

SR_dir = os.path.join(os.environ['GIT_AGN_MOCK'], 'data', 'scaling_relations')

itp_z, itp_kt, itp_frac_obs = n.loadtxt( os.path.join( os.environ['GIT_AGN_MOCK'], "data", "xray_k_correction", "fraction_05_20_01_24_no_nH.txt"), unpack=True )

nh_vals = 10**n.arange(-2,4+0.01,0.5)#0.05)
z_vals = n.hstack(( n.arange(0.,0.7,0.05), n.arange(0.8, 4.5, 0.1)))#[0.8, 0.9, 1, 1.1, 1.2, 1.4, 1.6] ))
#kT_vals = n.hstack(( n.arange(0.5,8,0.5), [10, 20, 30, 40, 50] ))
kT_vals = n.hstack(([0.1, 0.2], n.arange(0.5,8,0.5), [10, 20, 30, 40, 50] ))

XX_nh, YY_z, ZZ_kt = n.meshgrid(nh_vals, z_vals, kT_vals)

shape_i = XX_nh.shape

matrix_z_nh_kt = itp_frac_obs.reshape(shape_i)

from scipy.interpolate import RegularGridInterpolator
attenuation_3d = RegularGridInterpolator((z_vals, n.log10(nh_vals*1e22), kT_vals), matrix_z_nh_kt)

kT_kev, kT_kev_err, Rspec, R500, R500_err, M500, M500_err, Mgas500, Mgas500_err, R2500,  R2500_err, M2500, M2500_err, Mgas2500, Mgas2500_err, tcool, tcool_err, LXxmm, LXxmm_err = n.loadtxt(os.path.join(SR_dir, 'lovisari_2015_table2.ascii'), unpack=True)

redshift_s18, nh, LX_s18, Rktmax, M500NFWFreeze, M500kTextrp, M500NFWHudson, M500NFWAll, M200NFWFreeze, M200kTextrp, M500PlanckSZ = n.loadtxt(os.path.join(SR_dir, 'schllenberger_2018_tableB2B3.ascii'), unpack=True)

s_mi20, l_mi20, b_mi20, T_mi20, LX_mi20, sig_LX_mi20, f_mi20, NH_mi20, Z_mi20 = n.loadtxt(os.path.join(SR_dir, 'migkas_2020_tableC1.ascii'), unpack=True)

B18_id,    B18_z, B18_R500, B18_LXcin, B18_LXcinbol, B18_TXcin, B18_ZXcin, B18_LXcexbol, B18_LXcex, B18_TXcex, B18_ZXcex, B18_MICM, B18_YXcin, B18_M500 = n.loadtxt(os.path.join(SR_dir, 'bulbul_2018_table1_2.ascii'), unpack=True)

Lo20_planckName, Lo20_z, Lo20_M500, Lo20_Mg500, Lo20_kT, Lo20_kTexc, Lo20_LX, Lo20_LXexc, Lo20_Lbol, Lo20_Lbolexc, Lo20_NT, Lo20_fT, Lo20_Nsb, Lo20_fsb = n.loadtxt(os.path.join(SR_dir, 'lovisari_2020_tableA1.ascii'), unpack=True)

XXL_i = fits.open(os.path.join(os.environ['GIT_AGN_MOCK'],'data/XXL/xxl365gc.fits'))[1].data
XXL = XXL_i[(XXL_i['Class']==1)]

WtG = fits.open(os.path.join(os.environ['GIT_AGN_MOCK'], 'data','WtG','Mantz16_Table2.fits'))[1].data

fig_dir = os.path.join( os.environ['GIT_STACK'], 'figures', 'Ti20_SDSS_stacked_galaxy_profile', 'scaling_relations' )

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoUNIT = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
h = 0.6774
L_box = 1000.0 / h
cosmo = cosmoUNIT

Z_SEL = 'Z3P'
prof_dir = os.path.join( os.environ['GIT_STACK'], 'figures', 'Ti20_SDSS_stacked_galaxy_profile', 'profiles_decomp_AXH' )
MERGE_Mh = Table.read(os.path.join(prof_dir, 'FULL_SUMMARY_Mhalobin_centrals_'+Z_SEL+'.fits'))
t_mh = MERGE_Mh[n.unique(MERGE_Mh['file_ID'], return_index=True)[1]]
t_mh = t_mh [(t_mh['N_gal'] > 1000)&(t_mh['M_min'] >= 12.0)&(t_mh['M_max'] <=14.0)]
t = t_mh


#p_2_CLU = os.path.join( os.environ['HOME'], 'sf_Shared/data/UNIT_fA1i', 'UNIT_fA1i_DIR_eRO_CLU_b8_CM_0_pixS_20.0_M500c_13.0_FX_-14.5_MGAS_Sept2021.fits')
#CLU = Table.read(p_2_CLU)
#s_C = (CLU['redshift_R']>0.01) & (CLU['redshift_R']<0.36) #& (CLU['CLUSTER_FX_soft']>FX_MIN)& (CLU['CLUSTER_FX_soft']<factor*FX_MIN)
#x_clu = CLU['HALO_M500c'][s_C]
#y_clu = CLU['CLUSTER_LX_soft'][s_C] #+ n.log10(frac_op20)
#z_clu = CLU['redshift_R'][s_C]
#xx, yy = x_clu*cosmo.efunc(z_clu), 10**(y_clu-39)/cosmo.efunc(z_clu)

percents = n.arange(0,101,1)
def get_mean_scaling_relation(mass_array_log10, mass_proxy_array):#, percents = percents ):
	m_array = n.arange(13, 15, 0.2)
	mass_mins = m_array[:-1]
	mass_maxs = m_array[1:]
	# define empty columns
	mean_mass_proxy = n.zeros_like(mass_mins)
	std_mass_proxy  = n.zeros_like(mass_mins)
	mean_mass       = n.zeros_like(mass_mins)
	#percentile_values = n.zeros((len(mass_mins), len(percents)))
	# in each bin compute the mean and std of the mass proxy, mass and redshift
	for jj, (m_min, m_max) in enumerate(zip(mass_mins, mass_maxs)):
		selection = ( mass_array_log10 >= m_min ) & ( mass_array_log10 < m_max )
		mproxy_values = mass_proxy_array[ selection ]
		mm_values = mass_array_log10[ selection ]
		mean_mass_proxy[jj] = n.mean( mproxy_values )
		std_mass_proxy [jj] = n.std( mproxy_values )
		#percentile_values[jj] = scoreatpercentile(mproxy_values, percents)
		mean_mass      [jj] = n.mean(mm_values )

	## fit scaling relation for z<1.5
	#good = (mean_mass>0)&(mean_mass_proxy>0)&(std_mass_proxy>0.001)&(redshift_mins < z_max_SR)
	#xdata = mean_mass[good]
	#ydata = mean_mass_proxy[good]
	#yerr = std_mass_proxy[good]/2.
	#xerr = Dlog10M/2.
	#out_m = n.polyfit(xdata, ydata, deg=1, cov=True)#, w=1/yerr)
	#x_model = n.arange(13, 16, 0.1)
	#y_model = n.polyval(out_m[0], x_model)
	#SR.mass_mins  = mass_mins
	#SR.mass_maxs, redshift_mins = mass_maxs, redshift_mins
	#SR.redshift_maxs = redshift_maxs
	#SR.mean_mass_proxy = mean_mass_proxy
	#SR.std_mass_proxy = std_mass_proxy
	#SR.mean_mass = mean_mass
	#SR.mean_redshift = mean_redshift
	#SR.x_model = x_model
	#SR.y_model = y_model
	#SR.out_m = out_m
	#SR.percentile_values = percentile_values
	# output values and boundaries or std
	return mass_mins, mean_mass, mass_maxs, mean_mass_proxy, std_mass_proxy

#mass_mins, mean_mass, mass_maxs, mean_mass_proxy, std_mass_proxy = get_mean_scaling_relation(n.log10(xx), yy)

efeds = Table.read('/home/comparat/sf_Shared/data/erosita/observations/eFEDS_c001_clean/efeds_clusters_full_20210814.fits.gz')
sz = ( efeds['z_final']>0.01 ) & ( efeds['z_final']<0.35 )&(efeds['M500']*cosmo.efunc(efeds['z_final'])>4e13)
efeds = efeds[sz]

fig_out = os.path.join(fig_dir, 'M500c-LX-z.png')
#title_str=r'$M_{500c}>1\times10^{13}M_\odot$'

p.figure(0, (5.5, 5.))
# Co20
#p.fill_between(10**mean_mass, mean_mass_proxy-std_mass_proxy*0.5, mean_mass_proxy+std_mass_proxy*0.5, color='g', alpha=0.6, label='Co 20', rasterized = True)

# popesso 23
z_po23 = 0.1
po23_x_m500  = n.array([7.5e12, 2.2e13])*cosmo.efunc(z_po23)
po23_y_LX500 = n.array([2.8e41, 1.0e42])/1e39/cosmo.efunc(z_po23)
po23_yLO_LX500 = n.array([1e38, 2.2e41])/1e39/cosmo.efunc(z_po23)
po23_yUP_LX500 = n.array([5.5e41, 1.8e42])/1e39/cosmo.efunc(z_po23)
p.errorbar( po23_x_m500 ,
	po23_y_LX500 ,
	yerr= [po23_y_LX500-po23_yLO_LX500 , po23_yUP_LX500-po23_y_LX500],
	marker='v', ls='', mfc='none', label='Po23', markersize=10)

# Bulbul 18
p.plot(B18_M500*1e14*cosmo.efunc(B18_z),
	B18_LXcin * 1e44 / cosmo.efunc(B18_z)/1e39,
	marker='*', ls='', mfc='none', label='Bu19')

# Mantz 16
k_correction_3d_Ma16 = attenuation_3d( n.transpose([WtG['Ma16_z'], n.ones_like(WtG['Ma16_z'])*20.1, WtG['Ma16_kT_keV']]))
p.plot(  WtG['Ma16_Mlen_1e15'] * 1e15 ,
	WtG['Ma16_LX_1e44'] * k_correction_3d_Ma16 * 1e44 / cosmo.efunc(WtG['Ma16_z'] ) /1e39 ,
	marker='s', ls='', mfc='none', label='Ma16')

# Liu 22
p.plot(efeds['M500']*cosmo.efunc(efeds['z_final']), efeds['L500'].T[1]/1e39/cosmo.efunc(efeds['z_final']), marker='^', ls='', mfc='none', label='Li22', rasterized = True)
#p.plot(efeds['M500']*cosmo.efunc(efeds['z_final']), efeds['L500_Bulbul2019'].T[1]/1e39/cosmo.efunc(efeds['z_final']), marker='^', ls='', mfc='none', label='Li22 B', rasterized = True)
#p.plot(efeds['M500']*cosmo.efunc(efeds['z_final']), efeds['L500_Sereno2020'].T[1]/1e39/cosmo.efunc(efeds['z_final']), marker='^', ls='', mfc='none', label='Li22 S', rasterized = True)
#p.plot(efeds['M500']*cosmo.efunc(efeds['z_final']), efeds['L500_Lovisari2015'].T[1]/1e39/cosmo.efunc(efeds['z_final']), marker='^', ls='', mfc='none', label='Li22 L', rasterized = True)

# Lovsari 2020
k_correction_3d_Lo20 = attenuation_3d( n.transpose([Lo20_z, n.ones_like(Lo20_z)*20.1, Lo20_kT]))
p.plot( Lo20_M500*1e14*cosmo.efunc(Lo20_z) ,
	Lo20_LX * k_correction_3d_Lo20 * 1e44 / cosmo.efunc(Lo20_z) /1e39  ,
	marker='o', ls='', mfc='none', label='Lo20')

# Adami 18, XXL
ok = (XXL['Mgas500kpc']*1e11*10**(1.1) * cosmo.efunc(XXL['z'])>4e13)
p.plot( XXL['Mgas500kpc'][ok]*1e11*10**(1.1) * cosmo.efunc(XXL['z'][ok]),
	XXL['LXXL500MT'][ok] * 1e42/cosmo.efunc(XXL['z'][ok]) /1e39,
	marker='o', ls='', label='Ad18', mfc='none')

# Lovisari 2015, groups
k_correction_3d_Lo15 = attenuation_3d( n.transpose([n.ones_like(kT_kev)*0.02, n.ones_like(kT_kev)*20.1, kT_kev]))
ok = (M500 * 1e13 / 0.7>4e13)
p.plot( M500[ok] * 1e13 / 0.7 ,
	LXxmm[ok] * k_correction_3d_Lo15[ok] * 1e43/1e39,
	marker='s', ls='', mfc='none', label='Lo15')

# Schellenberger 18
k_correction_3d_Sc17 = attenuation_3d( n.transpose([redshift_s18, n.ones_like(redshift_s18)*20.1, n.ones_like(redshift_s18)*2]))
ok = (M500NFWFreeze*1e14*cosmo.efunc(redshift_s18) >4e13)
p.plot(M500NFWFreeze[ok]*1e14*cosmo.efunc(redshift_s18[ok]),
	LX_s18[ok] * k_correction_3d_Sc17[ok] * 1e43 / cosmo.efunc(redshift_s18[ok])/1e39,
	marker='o', ls='', label='Sc17', mfc='none')


# ALL

#p.errorbar( t['M500c_mean']*cosmo.efunc(t['redshift_mean']), t['total_LX_R500c']/cosmo.efunc(t['redshift_mean'])/1e39,
			#yerr = [t['total_LX_R500c' ]/cosmo.efunc(t['redshift_mean'])/1e39 - t['total_LX_R500c_lo']/cosmo.efunc(t['redshift_mean'])/1e39, t['total_LX_R500c_up']/cosmo.efunc(t['redshift_mean'])/1e39-t['total_LX_R500c']/cosmo.efunc(t['redshift_mean'])/1e39],
			#xerr = [ t['M500c_mean']*cosmo.efunc(t['redshift_mean'])-t['M500c_Q05']*cosmo.efunc(t['redshift_mean']),t['M500c_Q95']*cosmo.efunc(t['redshift_mean'])-t['M500c_mean']*cosmo.efunc(t['redshift_mean'])],
			#lw=1, ls='', color='darkgreen' )
p.errorbar( t['M500c_mean']*cosmo.efunc(t['redshift_mean']), t['total_LX_R500c']/1e39/cosmo.efunc(t['redshift_mean']),
			yerr = [t['total_LX_R500c' ]/cosmo.efunc(t['redshift_mean'])/1e39-t['total_LX_R500c_lo']/cosmo.efunc(t['redshift_mean'])/1e39, t['total_LX_R500c_up']/cosmo.efunc(t['redshift_mean'])/1e39-t['total_LX_R500c']/cosmo.efunc(t['redshift_mean'])/1e39],
			xerr = t['M500c_std']*cosmo.efunc(t['redshift_mean']),
			lw=2, ls='', color='darkgreen', label='Total emission' )
# Hot gas
#p.errorbar( t['M500c_mean']*cosmo.efunc(t['redshift_mean']), t['hotgas_LX_R500c']/cosmo.efunc(t['redshift_mean'])/1e39,
			#yerr = [t['hotgas_LX_R500c' ]/cosmo.efunc(t['redshift_mean'])/1e39-t['hotgas_LX_R500c_lo']/cosmo.efunc(t['redshift_mean'])/1e39, t['hotgas_LX_R500c_up']/cosmo.efunc(t['redshift_mean'])/1e39-t['hotgas_LX_R500c']/cosmo.efunc(t['redshift_mean'])/1e39],
			#xerr = [ t['M500c_mean']*cosmo.efunc(t['redshift_mean'])-t['M500c_Q05']*cosmo.efunc(t['redshift_mean']),t['M500c_Q95']*cosmo.efunc(t['redshift_mean'])-t['M500c_mean']*cosmo.efunc(t['redshift_mean'])],
			#lw=1, ls='', color='orange' )
p.errorbar( t['M500c_mean']*cosmo.efunc(t['redshift_mean']), t['hotgas_LX_R500c']/cosmo.efunc(t['redshift_mean'])/1e39,
			yerr = [t['hotgas_LX_R500c' ]/cosmo.efunc(t['redshift_mean'])/1e39-t['hotgas_LX_R500c_lo']/cosmo.efunc(t['redshift_mean'])/1e39, t['hotgas_LX_R500c_up']/cosmo.efunc(t['redshift_mean'])/1e39-t['hotgas_LX_R500c']/cosmo.efunc(t['redshift_mean'])/1e39],
			xerr = t['M500c_std']*cosmo.efunc(t['redshift_mean']),
			lw=2, ls='', color='orange', label='Hot gas' )


#fun_ch22 = lambda M500 : n.log(3.36) + n.log(1e43) + (1.44 - 0.07*n.log(1.35/1.35) )*n.log(M500/1.4e14) + 2*n.log(cosmo.efunc(0.35)/cosmo.efunc(0.35))-0.51**n.log(1.35/1.35)
#fun_ch22 = lambda M500 : n.log(3.36) + n.log(1e43) + 1.44*n.log(M500/1.4e14)
#M500s = 10**n.arange(13,15,0.01)*cosmo.efunc(0.35)
#lx = n.e**(fun_ch22(M500s)) /1e39/cosmo.efunc(0.35)
#p.plot(M500s, lx, ls='dashed',color='k')
slope = (44-n.log10(3e41))/(n.log10(4e14)-13)
oo = n.log10(3e41)-slope*13
fun = lambda x : slope * x + oo
x_1 = n.arange(13.2,15,0.01)
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


