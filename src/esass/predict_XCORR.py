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

m0 = 10.0
z1 = 0.18
m0_str = str(np.round(m0,1))
z1_str = str(np.round(z1,1))
p_2_all_xcorr = {}
p_2_all_srvval = {}

basename = 'GAL_m'+m0_str+'_evGAS_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxG'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxG'].sort()
p_2_all_srvval['GxG'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxG']])
Merge = get_merge_wth( p_2_all_xcorr['GxG'][np.isin(p_2_all_srvval['GxG'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALxEVT_wtheta_prediction_GxG.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)

basename = 'GAL_m'+m0_str+'_evBKG_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxB'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxB'].sort()
p_2_all_srvval['GxB'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxB']])
Merge = get_merge_wth( p_2_all_xcorr['GxB'][np.isin(p_2_all_srvval['GxB'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALxEVT_wtheta_prediction_GxB.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)

basename = 'GAL_m'+m0_str+'_evAGN_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxA'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxA'].sort()
p_2_all_srvval['GxA'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxA']])
Merge = get_merge_wth( p_2_all_xcorr['GxA'][np.isin(p_2_all_srvval['GxA'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALxEVT_wtheta_prediction_GxA.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)

basename = 'GAL_m' + m0_str + '_evAGNevCLU_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxGA'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxGA'].sort()
p_2_all_srvval['GxGA'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxGA']])
Merge = get_merge_wth( p_2_all_xcorr['GxGA'][np.isin(p_2_all_srvval['GxGA'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALxEVT_wtheta_prediction_GxGA.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)

basename = 'GAL_m' + m0_str + '_evAGNevCLUevBKG_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxGAB'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxGAB'].sort()
p_2_all_srvval['GxGAB'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxGAB']])
Merge = get_merge_wth( p_2_all_xcorr['GxGAB'][np.isin(p_2_all_srvval['GxGAB'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALxEVT_wtheta_prediction_GxGAB.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)


m0 = 10.5
z1 = 0.26
m0_str = str(np.round(m0,1))
z1_str = str(np.round(z1,1))
p_2_all_xcorr = {}
p_2_all_srvval = {}

basename = 'GAL_m'+m0_str+'_evGAS_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxG'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxG'].sort()
p_2_all_srvval['GxG'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxG']])
Merge = get_merge_wth( p_2_all_xcorr['GxG'][np.isin(p_2_all_srvval['GxG'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_GALxEVT_wtheta_prediction_GxG.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)

basename = 'GAL_m'+m0_str+'_evBKG_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxB'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxB'].sort()
p_2_all_srvval['GxB'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxB']])
Merge = get_merge_wth( p_2_all_xcorr['GxB'][np.isin(p_2_all_srvval['GxB'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_GALxEVT_wtheta_prediction_GxB.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)

basename = 'GAL_m'+m0_str+'_evAGN_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxA'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxA'].sort()
p_2_all_srvval['GxA'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxA']])
Merge = get_merge_wth( p_2_all_xcorr['GxA'][np.isin(p_2_all_srvval['GxA'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_GALxEVT_wtheta_prediction_GxA.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)

basename = 'GAL_m' + m0_str + '_evAGNevCLU_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxGA'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxGA'].sort()
p_2_all_srvval['GxGA'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxGA']])
Merge = get_merge_wth( p_2_all_xcorr['GxGA'][np.isin(p_2_all_srvval['GxGA'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_GALxEVT_wtheta_prediction_GxGA.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)

basename = 'GAL_m' + m0_str + '_evAGNevCLUevBKG_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxGAB'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxGAB'].sort()
p_2_all_srvval['GxGAB'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxGAB']])
Merge = get_merge_wth( p_2_all_xcorr['GxGAB'][np.isin(p_2_all_srvval['GxGAB'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_GALxEVT_wtheta_prediction_GxGAB.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)



m0 = 11.0
z1 = 0.35
m0_str = str(np.round(m0,1))
z1_str = str(np.round(z1,1))
p_2_all_xcorr = {}
p_2_all_srvval = {}

basename = 'GAL_m'+m0_str+'_evGAS_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxG'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxG'].sort()
p_2_all_srvval['GxG'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxG']])
Merge = get_merge_wth( p_2_all_xcorr['GxG'][np.isin(p_2_all_srvval['GxG'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_GALxEVT_wtheta_prediction_GxG.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)

basename = 'GAL_m'+m0_str+'_evBKG_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxB'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxB'].sort()
p_2_all_srvval['GxB'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxB']])
Merge = get_merge_wth( p_2_all_xcorr['GxB'][np.isin(p_2_all_srvval['GxB'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_GALxEVT_wtheta_prediction_GxB.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)

basename = 'GAL_m'+m0_str+'_evAGN_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxA'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxA'].sort()
p_2_all_srvval['GxA'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxA']])
Merge = get_merge_wth( p_2_all_xcorr['GxA'][np.isin(p_2_all_srvval['GxA'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_GALxEVT_wtheta_prediction_GxA.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)

basename = 'GAL_m' + m0_str + '_evAGNevCLU_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxGA'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxGA'].sort()
p_2_all_srvval['GxGA'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxGA']])
Merge = get_merge_wth( p_2_all_xcorr['GxGA'][np.isin(p_2_all_srvval['GxGA'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_GALxEVT_wtheta_prediction_GxGA.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)

basename = 'GAL_m' + m0_str + '_evAGNevCLUevBKG_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxGAB'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxGAB'].sort()
p_2_all_srvval['GxGAB'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxGAB']])
Merge = get_merge_wth( p_2_all_xcorr['GxGAB'][np.isin(p_2_all_srvval['GxGAB'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_GALxEVT_wtheta_prediction_GxGAB.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)


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


# mix = lambda f_A, f_B, f_G : (PRED['m11.0_GxA']['wtheta']*f_A*PRED['m11.0_GxA']['N_data']+  PRED['m11.0_GxB']['wtheta']*f_B*PRED['m11.0_GxB']['N_data'] + PRED['m11.0_GxG']['wtheta']*f_G*PRED['m11.0_GxG']['N_data'])/PRED['m11.0_GxGAB']['N_data']
# mix (1,1,1)

m0 = 10.0
z1 = 0.18
m0_str = str(np.round(m0,1))
z1_str = str(np.round(z1,2))
str_title = r'$z<$'+z1_str+', '+m0_str+r'$<\log_{10}(M*[M_\odot])$'

p_2_2PCF_figure = os.path.join(dir_fig, 'GALxEVT_m'+m0_str+'.png')
plt.figure(33, (6.5,6.))

WTH = XCORR['m10.0']
plt.errorbar(WTH['theta'], WTH['wtheta'], yerr= WTH['wtheta']*0.05, label='Obs', lw=3, color='black', marker='*')
x_fit = XCORR['m10.0']['theta'][2:]
y_fit = XCORR['m10.0']['wtheta'][2:]
plt.plot(x_fit, y_fit, ls='', color='red', marker='*', zorder=80)

def mix_fit(x_fit, f_A, f_B, f_G):
	N_tot = f_A * PRED['m10.0_GxA']['N_data'] + f_B * PRED['m10.0_GxB']['N_data'] + f_G * PRED['m10.0_GxG']['N_data']
	out = (PRED['m10.0_GxA']['wtheta'] * f_A * PRED['m10.0_GxA']['N_data'] + PRED['m10.0_GxB']['wtheta'] * f_B *
	 PRED['m10.0_GxB']['N_data'] + PRED['m10.0_GxG']['wtheta'] * f_G * PRED['m10.0_GxG']['N_data']) / N_tot
	itp = interp1d(PRED['m10.0_GxG'  ]['theta_mid'], out)
	return itp(x_fit)

popt, pcov = curve_fit(mix_fit, x_fit, y_fit, bounds=(0, 1))

mix_f = lambda f_A, f_B : (PRED['m10.0_GxA']['wtheta']*f_A*PRED['m10.0_GxA']['N_data']+  PRED['m10.0_GxB']['wtheta']*f_B*PRED['m10.0_GxB']['N_data'] + PRED['m10.0_GxG']['wtheta']*f_A*PRED['m10.0_GxG']['N_data'])/PRED['m10.0_GxGAB']['N_data']
mix_N = lambda n_A, n_B, n_G : (PRED['m10.0_GxA']['wtheta']*n_A+  PRED['m10.0_GxB']['wtheta']*n_B + PRED['m10.0_GxG']['wtheta']*n_G)/PRED['m10.0_GxGAB']['N_data']
f_A = 4
#f_B = 0.91
n_A = f_A*(PRED['m10.0_GxA']['N_data'][0]+PRED['m10.0_GxG']['N_data'][0])
n_B = PRED['m10.0_GxGAB']['N_data'][0] - n_A
f_B = n_B/PRED['m10.0_GxB']['N_data'][0]
#f_B*PRED['m10.0_GxB']['N_data'][0]
#n_G = PRED['m10.0_GxGAB']['N_data'][0]-n_B-n_A
print(f_A, f_B)
out = mix_f(f_A, f_B)
# print(out/PRED['m10.0_GxGAB']['wtheta'])


WTH = PRED['m10.0_GxG'  ]
plt.plot(WTH['theta_mid'], WTH['wtheta'], label='GxG', lw=1)
WTH = PRED['m10.0_GxB'  ]
plt.plot(WTH['theta_mid'], WTH['wtheta'], label='GxB', lw=1)
WTH = PRED['m10.0_GxA'  ]
plt.plot(WTH['theta_mid'], WTH['wtheta'], label='GxA', lw=1)
WTH = PRED['m10.0_GxGA' ]
plt.plot(WTH['theta_mid'], WTH['wtheta'], label='GxGA', lw=1)
WTH = PRED['m10.0_GxGAB']

plt.plot(WTH['theta_mid'], WTH['wtheta'], label='GxGAB', lw=1)
#
# mix_f = lambda f_A, f_B : (PRED['m10.0_GxA']['wtheta']*f_A*PRED['m10.0_GxA']['N_data']+  PRED['m10.0_GxB']['wtheta']*f_B*PRED['m10.0_GxB']['N_data'] + PRED['m10.0_GxG']['wtheta']*f_A*PRED['m10.0_GxG']['N_data'])/PRED['m10.0_GxGAB']['N_data']
# f_A = 2
# n_A = f_A*(PRED['m10.0_GxA']['N_data'][0]+PRED['m10.0_GxG']['N_data'][0])
# n_B = PRED['m10.0_GxGAB']['N_data'][0] - n_A
# f_B = n_B/PRED['m10.0_GxB']['N_data'][0]
# out = mix_f(f_A, f_B)
# plt.plot(WTH['theta_mid'], out, label='GxGA('+str(np.round(f_A,2))+')B('+str(np.round(f_B,2))+')', lw=2)
# f_A = 3
# n_A = f_A*(PRED['m10.0_GxA']['N_data'][0]+PRED['m10.0_GxG']['N_data'][0])
# n_B = PRED['m10.0_GxGAB']['N_data'][0] - n_A
# f_B = n_B/PRED['m10.0_GxB']['N_data'][0]
# out = mix_f(f_A, f_B)
# plt.plot(WTH['theta_mid'], out, label='GxGA('+str(np.round(f_A,2))+')B('+str(np.round(f_B,2))+')', lw=2)
# f_A = 4
# n_A = f_A*(PRED['m10.0_GxA']['N_data'][0]+PRED['m10.0_GxG']['N_data'][0])
# n_B = PRED['m10.0_GxGAB']['N_data'][0] - n_A
# f_B = n_B/PRED['m10.0_GxB']['N_data'][0]
# out = mix_f(f_A, f_B)
# plt.plot(WTH['theta_mid'], out, label='GxGA('+str(np.round(f_A,2))+')B('+str(np.round(f_B,2))+')', lw=2)

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

plt.plot(x_fit, mix_fit(x_fit, *popt), ls='--', color='darkred', lw=3, zorder=100, label='fit')#: f_A=%5.3f, f_B=%5.3f, f_G=%5.3f' % tuple(popt))
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xscale('log')
plt.yscale('log')
#plt.ylim((1e-3, 2e3))
#plt.xlim((8e-4, 3))
plt.legend(loc='lower left', ncol=2, fontsize=10, title=str_title)
plt.xlabel(r'$\theta$ (deg)',fontsize=18)
plt.ylabel(r'$w(\theta)$',fontsize=18)
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

p_2_2PCF_figure = os.path.join(dir_fig, 'GALxEVT_m'+m0_str+'.png')
plt.figure(33, (6.5,6.))

WTH = XCORR['m10.5']
plt.errorbar(WTH['theta'], WTH['wtheta'], yerr= WTH['wtheta']*0.05, label='Obs', lw=3, color='black', marker='*')

x_fit = XCORR['m10.5']['theta'][2:]
y_fit = XCORR['m10.5']['wtheta'][2:]
plt.plot(x_fit, y_fit, ls='', color='red', marker='*', zorder=80)

def mix_fit(x_fit, f_A, f_B, f_G):
	N_tot = f_A * PRED['m10.5_GxA']['N_data'] + f_B * PRED['m10.5_GxB']['N_data'] + f_G * PRED['m10.5_GxG']['N_data']
	out = (PRED['m10.5_GxA']['wtheta'] * f_A * PRED['m10.5_GxA']['N_data'] + PRED['m10.5_GxB']['wtheta'] * f_B *
	 PRED['m10.5_GxB']['N_data'] + PRED['m10.5_GxG']['wtheta'] * f_G * PRED['m10.5_GxG']['N_data']) / N_tot
	itp = interp1d(PRED['m10.5_GxG'  ]['theta_mid'], out)
	return itp(x_fit)

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

plt.plot(x_fit, mix_fit(x_fit, *popt), ls='--', color='darkred', lw=3, zorder=100, label='fit')#: f_A=%5.3f, f_B=%5.3f, f_G=%5.3f' % tuple(popt))

WTH = PRED['m10.5_GxG'  ]
plt.plot(WTH['theta_mid'], WTH['wtheta'], label='GxG', lw=1)
WTH = PRED['m10.5_GxB'  ]
plt.plot(WTH['theta_mid'], WTH['wtheta'], label='GxB', lw=1)
WTH = PRED['m10.5_GxA'  ]
plt.plot(WTH['theta_mid'], WTH['wtheta'], label='GxA', lw=1)
WTH = PRED['m10.5_GxGA' ]
plt.plot(WTH['theta_mid'], WTH['wtheta'], label='GxGA', lw=1)
WTH = PRED['m10.5_GxGAB']
plt.plot(WTH['theta_mid'], WTH['wtheta'], label='GxGAB', lw=1)

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xscale('log')
plt.yscale('log')
#plt.ylim((1e-3, 2e3))
#plt.xlim((8e-4, 3))
plt.legend(loc='lower left', ncol=2, fontsize=10, title=str_title)
plt.xlabel(r'$\theta$ (deg)',fontsize=18)
plt.ylabel(r'$w(\theta)$',fontsize=18)
plt.title('Correlation Galaxies x 0.5-2 keV events')
plt.tight_layout()
plt.savefig(p_2_2PCF_figure)
plt.clf()
print(p_2_2PCF_figure, 'written')





m0 = 11.0
z1 = 0.35
m0_str = str(np.round(m0,1))
z1_str = str(np.round(z1,2))
str_title = r'$z<$'+z1_str+', '+m0_str+r'$<\log_{10}(M*[M_\odot])$'


p_2_2PCF_figure = os.path.join(dir_fig, 'GALxEVT_m'+m0_str+'.png')
plt.figure(33, (6.5,6.))

WTH = XCORR['m11.0']
plt.errorbar(WTH['theta'], WTH['wtheta'], yerr= WTH['wtheta']*0.05, label='Obs', lw=3, color='black', marker='*')

x_fit = XCORR['m11.0']['theta'][2:]
y_fit = XCORR['m11.0']['wtheta'][2:]
plt.plot(x_fit, y_fit, ls='', color='red', marker='*', zorder=80)

def mix_fit(x_fit, f_A, f_B, f_G):
	N_tot = f_A * PRED['m11.0_GxA']['N_data'] + f_B * PRED['m11.0_GxB']['N_data'] + f_G * PRED['m11.0_GxG']['N_data']
	out = (PRED['m11.0_GxA']['wtheta'] * f_A * PRED['m11.0_GxA']['N_data'] + PRED['m11.0_GxB']['wtheta'] * f_B *
	 PRED['m11.0_GxB']['N_data'] + PRED['m11.0_GxG']['wtheta'] * f_G * PRED['m11.0_GxG']['N_data']) / N_tot
	itp = interp1d(PRED['m11.0_GxG'  ]['theta_mid'], out)
	return itp(x_fit)

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

plt.plot(x_fit, mix_fit(x_fit, *popt), ls='--', color='darkred', lw=3, zorder=100, label='fit')#: f_A=%5.3f, f_B=%5.3f, f_G=%5.3f' % tuple(popt))

WTH = PRED['m11.0_GxG'  ]
plt.plot(WTH['theta_mid'], WTH['wtheta'], label='GxG', lw=1)
WTH = PRED['m11.0_GxB'  ]
plt.plot(WTH['theta_mid'], WTH['wtheta'], label='GxB', lw=1)
WTH = PRED['m11.0_GxA'  ]
plt.plot(WTH['theta_mid'], WTH['wtheta'], label='GxA', lw=1)
WTH = PRED['m11.0_GxGA' ]
plt.plot(WTH['theta_mid'], WTH['wtheta'], label='GxGA', lw=1)
WTH = PRED['m11.0_GxGAB']
plt.plot(WTH['theta_mid'], WTH['wtheta'], label='GxGAB', lw=1)
# yerr =  abs(WTH['wtheta'] * (0.01**2 + WTH['D1D2_counts']**-1. + WTH['D1R2_counts']**-1. + WTH['D2R1_counts']**-1. + WTH['R1R2_counts']**-1. )**0.5)

# mix_f = lambda f_A, f_B : (PRED['m11.0_GxA']['wtheta']*f_A*PRED['m11.0_GxA']['N_data']+  PRED['m11.0_GxB']['wtheta']*f_B*PRED['m11.0_GxB']['N_data'] + PRED['m11.0_GxG']['wtheta']*f_A*PRED['m11.0_GxG']['N_data'])/PRED['m11.0_GxGAB']['N_data']
# f_A = 1.2
# n_A = f_A*(PRED['m11.0_GxA']['N_data'][0]+PRED['m11.0_GxG']['N_data'][0])
# n_B = PRED['m11.0_GxGAB']['N_data'][0] - n_A
# f_B = n_B/PRED['m11.0_GxB']['N_data'][0]
# out = mix_f(f_A, f_B)
# plt.plot(WTH['theta_mid'], out, label='GxGA('+str(np.round(f_A,2))+')B('+str(np.round(f_B,2))+')', lw=2)
# f_A = 1.5
# n_A = f_A*(PRED['m11.0_GxA']['N_data'][0]+PRED['m11.0_GxG']['N_data'][0])
# n_B = PRED['m11.0_GxGAB']['N_data'][0] - n_A
# f_B = n_B/PRED['m11.0_GxB']['N_data'][0]
# out = mix_f(f_A, f_B)
# plt.plot(WTH['theta_mid'], out, label='GxGA('+str(np.round(f_A,2))+')B('+str(np.round(f_B,2))+')', lw=2)
# f_A = 2
# n_A = f_A*(PRED['m11.0_GxA']['N_data'][0]+PRED['m11.0_GxG']['N_data'][0])
# n_B = PRED['m11.0_GxGAB']['N_data'][0] - n_A
# f_B = n_B/PRED['m11.0_GxB']['N_data'][0]
# out = mix_f(f_A, f_B)
# plt.plot(WTH['theta_mid'], out, label='GxGA('+str(np.round(f_A,2))+')B('+str(np.round(f_B,2))+')', lw=2)

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xscale('log')
plt.yscale('log')
#plt.ylim((1e-3, 2e3))
#plt.xlim((8e-4, 3))
plt.legend(loc='lower left', ncol=2, fontsize=10, title=str_title)
plt.xlabel(r'$\theta$ (deg)',fontsize=18)
plt.ylabel(r'$w(\theta)$',fontsize=18)
plt.title('Correlation Galaxies x 0.5-2 keV events')
plt.tight_layout()
plt.savefig(p_2_2PCF_figure)
plt.clf()
print(p_2_2PCF_figure, 'written')

