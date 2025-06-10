import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as plt
import os, sys
import glob
import numpy as np
from astropy.table import Table, Column, vstack, hstack
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
cosmoUNIT = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
h = 0.6774
L_box = 1000.0 / h
cosmo = cosmoUNIT
sky_map_hdu = Table.read(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS.fits') )
is_eroDE = ( (sky_map_hdu['OWNER']==2)|(sky_map_hdu['OWNER']==0) ) & ( abs(sky_map_hdu['GLAT_CEN']) > 20 ) #& ( sky_map_hdu['DE_CEN'] <= 32 )
SRVMAP_exGAL_eroDE = sky_map_hdu['SRVMAP'][is_eroDE]
benchmark_dir = os.path.join( os.environ['GIT_STMOD_DATA'], 'data/benchmark/xcorr_comparat_2025' )

def get_merge_wth(p_2_Apcf):
	PCF = []
	for el in p_2_Apcf:
		t_i = Table.read(el)
		#print(len(t_i), el)
		#if len(t_i)==73:
		PCF.append(t_i)

	PCF = np.array(PCF)

	Merge = Table()
	Merge['theta_min'] = PCF[-1]['theta_min']
	Merge['theta_max'] = PCF[-1]['theta_max']
	Merge['theta_mid'] = PCF[-1]['theta_mid']
	Merge['wtheta'] = np.zeros_like(PCF['wtheta'].sum(axis=0))
	Merge['N_data'] = PCF['N_data'].sum(axis=0)
	Merge['N2_data'] = PCF['N2_data'].sum(axis=0)
	Merge['N_random'] = PCF['N_random'].sum(axis=0)
	Merge['N2_random'] = PCF['N_random'].sum(axis=0)
	Merge['D1D2_counts'] = PCF['D1D2_counts'].sum(axis=0)
	Merge['D1R2_counts'] = PCF['D1R_counts'].sum(axis=0)
	Merge['D2R1_counts'] = PCF['D2R_counts'].sum(axis=0)
	Merge['R1R2_counts'] = PCF['RR_counts'].sum(axis=0)

	fN1 = Merge['N_random'][0] / Merge['N_data'][0]
	fN2 = Merge['N2_random'][0] / Merge['N2_data'][0]
	cf = np.zeros(len(Merge))
	cf[:] = np.nan
	nonzero = Merge['R1R2_counts'] > 0
	cf[nonzero] = (fN1 * fN2 * Merge['D1D2_counts'][nonzero] -
					fN1 * Merge['D1R2_counts'][nonzero] -
					fN2 * Merge['D2R1_counts'][nonzero] +
					Merge['R1R2_counts'][nonzero]) / Merge['R1R2_counts'][nonzero]
	Merge['wtheta'] = cf
	return Merge

agn_seed = '1' # sys.argv[1] # 1
clu_seed = '1' # sys.argv[2] # 1
LC_dir = 'LCerass'

dir_2pcf = os.path.join(os.environ['UCHUU'], LC_dir, '??????',
                         'GE_e4_merge_AGNseed' + agn_seed.zfill(3) + '_SimBKG_CLUseed' + clu_seed.zfill(3), 'XCORR')

XCORR = {}
XCORR['m10.0'] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALxEVT_wtheta.fits' ) )
XCORR['m10.5'] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_GALxEVT_wtheta.fits' ) )
XCORR['m11.0'] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_GALxEVT_wtheta.fits' ) )

PRED = {}

PRED['m10.0_GxG'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALxEVT_wtheta_prediction_GxG.fits'))
PRED['m10.0_GxB'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALxEVT_wtheta_prediction_GxB.fits'))
PRED['m10.0_GxA'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALxEVT_wtheta_prediction_GxA.fits'))
PRED['m10.0_GxGA' ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALxEVT_wtheta_prediction_GxGA.fits'))
PRED['m10.0_GxGAB'] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALxEVT_wtheta_prediction_GxGAB.fits'))

PRED['m10.0C_GxG'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALCENxEVT_wtheta_prediction_GxG.fits'))
PRED['m10.0C_GxB'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALCENxEVT_wtheta_prediction_GxB.fits'))
PRED['m10.0C_GxA'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALCENxEVT_wtheta_prediction_GxA.fits'))
PRED['m10.0C_GxGA' ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALCENxEVT_wtheta_prediction_GxGA.fits'))
PRED['m10.0C_GxGAB'] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALCENxEVT_wtheta_prediction_GxGAB.fits'))

PRED['m10.0S_GxG'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALSATxEVT_wtheta_prediction_GxG.fits'))
PRED['m10.0S_GxB'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALSATxEVT_wtheta_prediction_GxB.fits'))
PRED['m10.0S_GxA'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALSATxEVT_wtheta_prediction_GxA.fits'))
PRED['m10.0S_GxGA' ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALSATxEVT_wtheta_prediction_GxGA.fits'))
PRED['m10.0S_GxGAB'] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALSATxEVT_wtheta_prediction_GxGAB.fits'))

PRED['m10.5_GxG'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_GALxEVT_wtheta_prediction_GxG.fits'))
PRED['m10.5_GxB'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_GALxEVT_wtheta_prediction_GxB.fits'))
PRED['m10.5_GxA'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_GALxEVT_wtheta_prediction_GxA.fits'))
PRED['m10.5_GxGA' ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_GALxEVT_wtheta_prediction_GxGA.fits'))
PRED['m10.5_GxGAB'] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_GALxEVT_wtheta_prediction_GxGAB.fits'))
PRED['m11.0_GxG'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_GALxEVT_wtheta_prediction_GxG.fits'))
PRED['m11.0_GxB'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_GALxEVT_wtheta_prediction_GxB.fits'))
PRED['m11.0_GxA'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_GALxEVT_wtheta_prediction_GxA.fits'))
PRED['m11.0_GxGA' ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_GALxEVT_wtheta_prediction_GxGA.fits'))
PRED['m11.0_GxGAB'] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_GALxEVT_wtheta_prediction_GxGAB.fits'))

dir_fig = os.path.join( os.environ['GIT_STMOD_DATA'], 'data/validation/validation_GasGal/XCORR' )
os.system('mkdir -p '+dir_fig)

m0 = 10.0
z1 = 0.18
m0_str = str(np.round(m0,1))
z1_str = str(np.round(z1,2))
str_title = r'$z<$'+z1_str+', '+m0_str+r'$<\log_{10}(M*[M_\odot])$'
z_mean = 0.136
BG_val = 6.859*10**(36)*(1+z_mean)**2
cv_r = cosmo.kpc_proper_per_arcmin(z_mean).to(u.kpc/u.deg).value

p_2_2PCF_figure = os.path.join(dir_fig, 'GALxEVT_m'+m0_str+'.png')
plt.figure(33, (6.5,6.))

WTH = XCORR['m10.0']
plt.errorbar(WTH['theta']*cv_r, WTH['wtheta']*BG_val, yerr= WTH['wtheta']*0.05, label='Obs', lw=3, color='black', marker='*')
x_fit = XCORR['m10.0']['theta'][6:]
y_fit = XCORR['m10.0']['wtheta'][6:]

def mix_fit(x_fit, f_A, f_B, f_G):
	N_tot = f_A * PRED['m10.0_GxA']['N_data'] + f_B * PRED['m10.0_GxB']['N_data'] + f_G * PRED['m10.0_GxG']['N_data']
	out = (PRED['m10.0_GxA']['wtheta'] * f_A * PRED['m10.0_GxA']['N_data'] + np.median(PRED['m10.0_GxB']['wtheta']) * f_B *
	 PRED['m10.0_GxB']['N_data'] + PRED['m10.0_GxG']['wtheta'] * f_G * PRED['m10.0_GxG']['N_data']) / N_tot
	itp = interp1d(PRED['m10.0_GxG'  ]['theta_mid'], out)
	return itp(x_fit)

def mix_fit_components(x_fit, f_A, f_B, f_G):
	N_tot = f_A * PRED['m10.0_GxA']['N_data'] + f_B * PRED['m10.0_GxB']['N_data'] + f_G * PRED['m10.0_GxG']['N_data']
	out = np.array([ (PRED['m10.0_GxA']['wtheta'] * f_A * PRED['m10.0_GxA']['N_data'] , np.median(PRED['m10.0_GxB']['wtheta']) * f_B *
	 PRED['m10.0_GxB']['N_data'] , PRED['m10.0_GxG']['wtheta'] * f_G * PRED['m10.0_GxG']['N_data']) ]) / N_tot
	return out

print('fitting 1')
popt, pcov = curve_fit(mix_fit, x_fit, y_fit, bounds=(0, 1))
print('popt', popt)
f_A_nat = PRED['m10.0_GxA']['N_data'][0]/ PRED['m10.0_GxGAB']['N_data'][0]
f_B_nat = PRED['m10.0_GxB']['N_data'][0]/ PRED['m10.0_GxGAB']['N_data'][0]
f_G_nat = PRED['m10.0_GxG']['N_data'][0]/ PRED['m10.0_GxGAB']['N_data'][0]

n_tot_fit = popt[0]*PRED['m10.0_GxA']['N_data'][0] + popt[1]*PRED['m10.0_GxB']['N_data'][0] + popt[2]*PRED['m10.0_GxG']['N_data'][0]
f_A_fit = popt[0]*PRED['m10.0_GxA']['N_data'][0]/ n_tot_fit
f_B_fit = popt[1]*PRED['m10.0_GxB']['N_data'][0]/ n_tot_fit
f_G_fit = popt[2]*PRED['m10.0_GxG']['N_data'][0]/ n_tot_fit
print('AGN',f_A_nat, f_A_fit)
print('BKG',f_B_nat, f_B_fit)
print('GAS',f_G_nat, f_G_fit)

out1, out2, out3 = mix_fit_components(x_fit, *popt)[0]
plt.plot(PRED['m10.0_GxA']['theta_mid']*cv_r, BG_val * (out1), label='GxA', lw=1, color='orange')
plt.plot(PRED['m10.0_GxB']['theta_mid']*cv_r, BG_val * (out2), label='GxB', lw=1, color='green')
plt.plot(PRED['m10.0_GxG']['theta_mid']*cv_r, BG_val * (out3), label='Gxg', lw=1, color='purple')
plt.plot(PRED['m10.0_GxG']['theta_mid']*cv_r, BG_val * (out1+out2+out3), label='Sum', lw=1, color='grey')

def mix_fit_FULL(x_fit, f_AC, f_AS, f_B, f_GC, f_GS):
	N_tot = f_AC * PRED['m10.0C_GxA']['N_data'] + f_AS * PRED['m10.0S_GxA']['N_data'] + f_B * PRED['m10.0_GxB']['N_data'] + f_GC * PRED['m10.0C_GxG']['N_data'] + f_GS * PRED['m10.0S_GxG']['N_data']
	out = (PRED['m10.0C_GxA']['wtheta'] * f_AC * PRED['m10.0C_GxA']['N_data'] +
		   PRED['m10.0S_GxA']['wtheta'] * f_AS * PRED['m10.0S_GxA']['N_data'] +
		   np.median(PRED['m10.0_GxB']['wtheta']) * f_B * PRED['m10.0_GxB']['N_data'] +
		   PRED['m10.0C_GxG']['wtheta'] * f_GC * PRED['m10.0C_GxG']['N_data']+
		   PRED['m10.0C_GxG']['wtheta'] * f_GS * PRED['m10.0S_GxG']['N_data']       ) / N_tot
	itp = interp1d(PRED['m10.0_GxG'  ]['theta_mid'], out)
	return itp(x_fit)

def mix_fit_components_FULL(x_fit, f_AC, f_AS, f_B, f_GC, f_GS):
	N_tot = f_AC * PRED['m10.0C_GxA']['N_data'] + f_AS * PRED['m10.0S_GxA']['N_data'] + f_B * PRED['m10.0_GxB']['N_data'] + f_GC * PRED['m10.0C_GxG']['N_data'] + f_GS * PRED['m10.0S_GxG']['N_data']
	out = (PRED['m10.0C_GxA']['wtheta'] * f_AC * PRED['m10.0C_GxA']['N_data'] ,
		   PRED['m10.0S_GxA']['wtheta'] * f_AS * PRED['m10.0S_GxA']['N_data'] ,
		   np.median(PRED['m10.0_GxB']['wtheta']) * f_B * PRED['m10.0_GxB']['N_data'] ,
		   PRED['m10.0C_GxG']['wtheta'] * f_GC * PRED['m10.0C_GxG']['N_data'],
		   PRED['m10.0C_GxG']['wtheta'] * f_GS * PRED['m10.0S_GxG']['N_data']       ) / N_tot
	return out

print('fitting 2')
popt, pcov = curve_fit(mix_fit_FULL, x_fit, y_fit, bounds=(0, 1))
print('popt more parameters', popt)

out1, out2, out3, out4, out5 = mix_fit_components_FULL(x_fit, *popt)
plt.plot(PRED['m10.0_GxA']['theta_mid']*cv_r, BG_val * (out1), label='GCxA', lw=1, color='orange', ls='dashed')
plt.plot(PRED['m10.0_GxB']['theta_mid']*cv_r, BG_val * (out2), label='GSxA', lw=1, color='orange', ls='dotted')
plt.plot(PRED['m10.0_GxG']['theta_mid']*cv_r, BG_val * (out3), label='GxB', lw=1, color='green', ls='dashed')
plt.plot(PRED['m10.0_GxG']['theta_mid']*cv_r, BG_val * (out4), label='GCxG', lw=1, color='purple', ls='dashed')
plt.plot(PRED['m10.0_GxG']['theta_mid']*cv_r, BG_val * (out5), label='GSxG', lw=1, color='purple', ls='dotted')
plt.plot(PRED['m10.0_GxG']['theta_mid']*cv_r, BG_val * (out1+out2+out3+out4+out5), label='Sum', lw=1, color='grey', ls='dashed')

# plt.plot(x_fit, mix_fit(x_fit, *popt), ls='--', color='darkred', lw=3, zorder=100, label='fit')#: f_A=%5.3f, f_B=%5.3f, f_G=%5.3f' % tuple(popt))
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xscale('log')
plt.yscale('log')
plt.xlim((7, 2000))
plt.ylim((1e33, 2e37))
plt.legend(loc='lower left', ncol=2, fontsize=10, title=str_title)
plt.xlabel(r'$r_p$ (Mpc)',fontsize=18)
plt.ylabel(r'$S_x$ (erg/kpc2/s)',fontsize=18)
plt.title('Correlation Galaxies x 0.5-2 keV events')
plt.tight_layout()
plt.savefig(p_2_2PCF_figure)
plt.clf()
print(p_2_2PCF_figure, 'written')




m0 = 10.5
z1 = 0.26
m0_str = str(np.round(m0,1))
z1_str = str(np.round(z1,2))
str_title = r'$z<$'+z1_str+', '+m0_str+r'$<\log_{10}(M*[M_\odot])$'
z_mean = 0.191
BG_val = 6.772*10**(36)*(1+z_mean)**2
cv_r = cosmo.kpc_proper_per_arcmin(z_mean).to(u.kpc/u.deg).value

p_2_2PCF_figure = os.path.join(dir_fig, 'GALxEVT_m'+m0_str+'.png')
plt.figure(33, (6.5,6.))

WTH = XCORR['m10.5']
plt.errorbar(WTH['theta']*cv_r, WTH['wtheta']*BG_val, yerr= WTH['wtheta']*0.05, label='Obs', lw=3, color='black', marker='*')
x_fit = XCORR['m10.5']['theta'][6:]
y_fit = XCORR['m10.5']['wtheta'][6:]

def mix_fit(x_fit, f_A, f_B, f_G):
	N_tot = f_A * PRED['m10.5_GxA']['N_data'] + f_B * PRED['m10.5_GxB']['N_data'] + f_G * PRED['m10.5_GxG']['N_data']
	out = (PRED['m10.5_GxA']['wtheta'] * f_A * PRED['m10.5_GxA']['N_data'] + np.median(PRED['m10.5_GxB']['wtheta']) * f_B *
	 PRED['m10.5_GxB']['N_data'] + PRED['m10.5_GxG']['wtheta'] * f_G * PRED['m10.5_GxG']['N_data']) / N_tot
	itp = interp1d(PRED['m10.5_GxG'  ]['theta_mid'], out)
	return itp(x_fit)

def mix_fit_components(x_fit, f_A, f_B, f_G):
	N_tot = f_A * PRED['m10.5_GxA']['N_data'] + f_B * PRED['m10.5_GxB']['N_data'] + f_G * PRED['m10.5_GxG']['N_data']
	out = np.array([ (PRED['m10.5_GxA']['wtheta'] * f_A * PRED['m10.5_GxA']['N_data'] , np.median(PRED['m10.5_GxB']['wtheta']) * f_B *
	 PRED['m10.5_GxB']['N_data'] , PRED['m10.5_GxG']['wtheta'] * f_G * PRED['m10.5_GxG']['N_data']) ]) / N_tot
	return out

popt, pcov = curve_fit(mix_fit, x_fit, y_fit, bounds=(0, 1))

f_A_nat = PRED['m10.5_GxA']['N_data'][0]/ PRED['m10.5_GxGAB']['N_data'][0]
f_B_nat = PRED['m10.5_GxB']['N_data'][0]/ PRED['m10.5_GxGAB']['N_data'][0]
f_G_nat = PRED['m10.5_GxG']['N_data'][0]/ PRED['m10.5_GxGAB']['N_data'][0]

n_tot_fit = popt[0]*PRED['m10.5_GxA']['N_data'][0] + popt[1]*PRED['m10.5_GxB']['N_data'][0] + popt[2]*PRED['m10.5_GxG']['N_data'][0]
f_A_fit = popt[0]*PRED['m10.5_GxA']['N_data'][0]/ n_tot_fit
f_B_fit = popt[1]*PRED['m10.5_GxB']['N_data'][0]/ n_tot_fit
f_G_fit = popt[2]*PRED['m10.5_GxG']['N_data'][0]/ n_tot_fit
print('AGN',f_A_nat, f_A_fit)
print('BKG',f_B_nat, f_B_fit)
print('GAS',f_G_nat, f_G_fit)

out1, out2, out3 = mix_fit_components(x_fit, *popt)[0]
plt.plot(PRED['m10.5_GxA']['theta_mid']*cv_r, BG_val * (out1), label='GxA', lw=1, color='orange')
plt.plot(PRED['m10.5_GxB']['theta_mid']*cv_r, BG_val * (out2), label='GxB', lw=1, color='green')
plt.plot(PRED['m10.5_GxG']['theta_mid']*cv_r, BG_val * (out3), label='Gxg', lw=1, color='purple')
plt.plot(PRED['m10.5_GxG']['theta_mid']*cv_r, BG_val * (out1+out2+out3), label='Sum', lw=1, color='grey')

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xscale('log')
plt.yscale('log')
plt.xlim((7, 2000))
plt.ylim((1e33, 2e37))
plt.legend(loc='lower left', ncol=2, fontsize=10, title=str_title)
plt.xlabel(r'$r_p$ (Mpc)',fontsize=18)
plt.ylabel(r'$S_x$ (erg/kpc2/s)',fontsize=18)
plt.title('Correlation Galaxies x 0.5-2 keV events')
plt.tight_layout()
plt.savefig(p_2_2PCF_figure)
plt.clf()
print(p_2_2PCF_figure, 'written')





m0 = 11.0
z1 = 0.35
m0_str = str(np.round(m0,1))
z1_str = str(np.round(z1,2))
z_mean = 0.252
BG_val = 6.693*10**(36)*(1+z_mean)**2
cv_r = cosmo.kpc_proper_per_arcmin(z_mean).to(u.kpc/u.deg).value

p_2_2PCF_figure = os.path.join(dir_fig, 'GALxEVT_m'+m0_str+'.png')
plt.figure(33, (6.5,6.))

WTH = XCORR['m11.0']
plt.errorbar(WTH['theta']*cv_r, WTH['wtheta']*BG_val, yerr= WTH['wtheta']*0.05, label='Obs', lw=3, color='black', marker='*')
x_fit = XCORR['m11.0']['theta'][6:]
y_fit = XCORR['m11.0']['wtheta'][6:]

def mix_fit(x_fit, f_A, f_B, f_G):
	N_tot = f_A * PRED['m11.0_GxA']['N_data'] + f_B * PRED['m11.0_GxB']['N_data'] + f_G * PRED['m11.0_GxG']['N_data']
	out = (PRED['m11.0_GxA']['wtheta'] * f_A * PRED['m11.0_GxA']['N_data'] + np.median(PRED['m11.0_GxB']['wtheta']) * f_B *
	 PRED['m11.0_GxB']['N_data'] + PRED['m11.0_GxG']['wtheta'] * f_G * PRED['m11.0_GxG']['N_data']) / N_tot
	itp = interp1d(PRED['m11.0_GxG'  ]['theta_mid'], out)
	return itp(x_fit)

def mix_fit_components(x_fit, f_A, f_B, f_G):
	N_tot = f_A * PRED['m11.0_GxA']['N_data'] + f_B * PRED['m11.0_GxB']['N_data'] + f_G * PRED['m11.0_GxG']['N_data']
	out = np.array([ (PRED['m11.0_GxA']['wtheta'] * f_A * PRED['m11.0_GxA']['N_data'] , np.median(PRED['m11.0_GxB']['wtheta']) * f_B *
	 PRED['m11.0_GxB']['N_data'] , PRED['m11.0_GxG']['wtheta'] * f_G * PRED['m11.0_GxG']['N_data']) ]) / N_tot
	return out

popt, pcov = curve_fit(mix_fit, x_fit, y_fit, bounds=(0, 1))

f_A_nat = PRED['m11.0_GxA']['N_data'][0]/ PRED['m11.0_GxGAB']['N_data'][0]
f_B_nat = PRED['m11.0_GxB']['N_data'][0]/ PRED['m11.0_GxGAB']['N_data'][0]
f_G_nat = PRED['m11.0_GxG']['N_data'][0]/ PRED['m11.0_GxGAB']['N_data'][0]

n_tot_fit = popt[0]*PRED['m11.0_GxA']['N_data'][0] + popt[1]*PRED['m11.0_GxB']['N_data'][0] + popt[2]*PRED['m11.0_GxG']['N_data'][0]
f_A_fit = popt[0]*PRED['m11.0_GxA']['N_data'][0]/ n_tot_fit
f_B_fit = popt[1]*PRED['m11.0_GxB']['N_data'][0]/ n_tot_fit
f_G_fit = popt[2]*PRED['m11.0_GxG']['N_data'][0]/ n_tot_fit
print('AGN',f_A_nat, f_A_fit)
print('BKG',f_B_nat, f_B_fit)
print('GAS',f_G_nat, f_G_fit)

out1, out2, out3 = mix_fit_components(x_fit, *popt)[0]
plt.plot(PRED['m11.0_GxA']['theta_mid']*cv_r, BG_val * (out1), label='GxA', lw=1, color='orange')
plt.plot(PRED['m11.0_GxB']['theta_mid']*cv_r, BG_val * (out2), label='GxB', lw=1, color='green')
plt.plot(PRED['m11.0_GxG']['theta_mid']*cv_r, BG_val * (out3), label='Gxg', lw=1, color='purple')
plt.plot(PRED['m11.0_GxG']['theta_mid']*cv_r, BG_val * (out1+out2+out3), label='Sum', lw=1, color='grey')

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xscale('log')
plt.yscale('log')
plt.xlim((7, 2000))
plt.ylim((1e33, 2e37))
plt.legend(loc='lower left', ncol=2, fontsize=10, title=str_title)
plt.xlabel(r'$r_p$ (Mpc)',fontsize=18)
plt.ylabel(r'$S_x$ (erg/kpc2/s)',fontsize=18)
plt.title('Correlation Galaxies x 0.5-2 keV events')
plt.tight_layout()
plt.savefig(p_2_2PCF_figure)
plt.clf()
print(p_2_2PCF_figure, 'written')

