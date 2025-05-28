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
cosmoUNIT = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
h = 0.6774
L_box = 1000.0 / h
cosmo = cosmoUNIT
sky_map_hdu = Table.read(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS.fits') )
is_eroDE = ( (sky_map_hdu['OWNER']==2)|(sky_map_hdu['OWNER']==0) ) & ( abs(sky_map_hdu['GLAT_CEN']) > 20 ) & ( sky_map_hdu['DE_CEN'] <= 32 )
SRVMAP_exGAL_eroDE = sky_map_hdu['SRVMAP'][is_eroDE]

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
benchmark_dir = os.path.join( os.environ['GIT_STMOD_DATA'], 'data/benchmark/xcorr_comparat_2025' )
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


m0 = 10.0
z1 = 0.18
m0_str = str(np.round(m0,1))
z1_str = str(np.round(z1,2))
str_title = r'$z<$'+z1_str+', '+m0_str+r'$<\log_{10}(M*[M_\odot])$'

p_2_2PCF_figure = os.path.join(dir_fig, 'GALxEVT_m'+m0_str+'.png')
plt.figure(33, (6.5,6.))

WTH = XCORR['m10.0']
plt.errorbar(WTH['theta'], WTH['wtheta'], yerr= WTH['wtheta']*0.05, label='Obs', lw=3, color='black', marker='*')

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
# yerr =  abs(WTH['wtheta'] * (0.01**2 + WTH['D1D2_counts']**-1. + WTH['D1R2_counts']**-1. + WTH['D2R1_counts']**-1. + WTH['R1R2_counts']**-1. )**0.5)

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
# yerr =  abs(WTH['wtheta'] * (0.01**2 + WTH['D1D2_counts']**-1. + WTH['D1R2_counts']**-1. + WTH['D2R1_counts']**-1. + WTH['R1R2_counts']**-1. )**0.5)

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

