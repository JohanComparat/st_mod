# mamba activate /home/idies/workspace/erosim/software/JC/env/clustering

import os
import glob
import numpy as np
from astropy.table import Table
os.environ['UCHUU']='/home/idies/workspace/erosim/Uchuu'
os.environ['GIT_STMOD']='/home/idies/workspace/erosim/software/st_mod'
os.environ['GIT_STMOD_DATA']='/home/idies/workspace/erosim/software/st_mod_data'

sky_map_hdu = Table.read(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS.fits') )
is_eroDE = ( (sky_map_hdu['OWNER']==2)|(sky_map_hdu['OWNER']==0) ) & ( abs(sky_map_hdu['GLAT_CEN']) > 20 ) #& ( sky_map_hdu['DE_CEN'] <= 32 )
SRVMAP_exGAL_eroDE = sky_map_hdu['SRVMAP'][is_eroDE]
is_eroDE30 = ( (sky_map_hdu['OWNER']==2)|(sky_map_hdu['OWNER']==0) ) & ( abs(sky_map_hdu['GLAT_CEN']) > 30 ) #& ( sky_map_hdu['DE_CEN'] <= 32 )
SRVMAP_exGAL_eroDE30 = sky_map_hdu['SRVMAP'][is_eroDE30]
benchmark_dir = os.path.join( os.environ['GIT_STMOD_DATA'], 'data/benchmark/xcorr_comparat_2025', 'XCORR_SHELLzlt04' )
os.system('mkdir -p ' + benchmark_dir)
# benchmark_dir = os.path.join( os.environ['GIT_STMOD_DATA'], 'data/benchmark/xcorr_comparat_2025', 'XCORR' )
# os.system('mkdir -p ' + benchmark_dir)


def get_merge_wth(p_2_Apcf):
	print('merging', len(p_2_Apcf), 'tiles')
	PCF = []
	for el in p_2_Apcf:
		print(el)
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
clu_seed = '19' # sys.argv[2] # 1
LC_dir = 'LCerass'

dir_2pcf = os.path.join(os.environ['UCHUU'], LC_dir, '??????',
                         'GE_e4_merge_AGNseed' + agn_seed.zfill(3) + '_SimBKG_CLUseed' + clu_seed.zfill(3), 'XCORR_SHELLzlt04')
# dir_2pcf = os.path.join(os.environ['UCHUU'], LC_dir, '??????',
#                          'GE_e4_merge_AGNseed' + agn_seed.zfill(3) + '_SimBKG_CLUseed' + clu_seed.zfill(3), 'XCORR')

m0 = 10.0
z1 = 0.18
m0_str = str(np.round(m0,1))
z1_str = str(np.round(z1,1))
p_2_all_xcorr = {}
p_2_all_srvval = {}

basename = 'GALCEN_m'+m0_str+'_evGAS_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxG'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxG'].sort()
p_2_all_srvval['GxG'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxG']])
Merge = get_merge_wth( p_2_all_xcorr['GxG'][np.isin(p_2_all_srvval['GxG'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALCENxEVT_wtheta_prediction_GxG.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)

basename = 'GALCEN_m'+m0_str+'_evBKG_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxB'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxB'].sort()

# MED_WTH_GxB = np.array([np.median(Table.read(t)['wtheta']) for t in p_2_all_xcorr['GxB']])
# import matplotlib
# matplotlib.use('Agg')
# matplotlib.rcParams.update({'font.size': 14})
# import matplotlib.pyplot as plt
# dir_fig = os.path.join( os.environ['GIT_STMOD_DATA'], 'data/validation/validation_GasGal/XCORR_SHELL' )

# p_2_2PCF_figure = os.path.join(dir_fig, 'BG_med_hist.png')
# plt.figure(11, (6.5,5.5))
# plt.hist(np.log10(MED_WTH_GxB[MED_WTH_GxB>0]), bins=20,  histtype='step', color='k', lw=2, label='pos')
# plt.hist(np.log10(-MED_WTH_GxB[MED_WTH_GxB<0]), bins=20,  histtype='step', color='b', lw=2, label='neg')
# plt.xticks(fontsize=14)
# plt.yticks(fontsize=14)
# # plt.xscale('log')
# plt.yscale('log')
# # plt.xlim((5, 5000))
# # plt.ylim((3e33, 3e38))
# plt.xlabel(r'log10 median GxBG ',fontsize=18)
# plt.ylabel(r'$N$ ',fontsize=18)
# plt.legend(fontsize=14)
# plt.tight_layout()
# plt.savefig(p_2_2PCF_figure)#, transparent = True)
# plt.clf()
# print(p_2_2PCF_figure, 'written')

# p_2_2PCF_figure = os.path.join(dir_fig, 'BG_med_hist-full.png')
# plt.figure(11, (6.5,5.5))
# plt.hist(MED_WTH_GxB,  bins=40, histtype='step', color='k', lw=4, label='all', density=True)
# plt.hist(MED_WTH_GxB[np.isin(p_2_all_srvval['GxB'], SRVMAP_exGAL_eroDE)],  bins=40, histtype='step', color='r', lw=2, label='glat>20', density=True)
# plt.hist(MED_WTH_GxB[np.isin(p_2_all_srvval['GxB'], SRVMAP_exGAL_eroDE30)],  bins=40, histtype='step', color='b', lw=2, label='glat>30', density=True)
# plt.xticks(fontsize=14)
# plt.yticks(fontsize=14)
# # plt.xscale('log')
# plt.yscale('log')
# # plt.xlim((5, 5000))
# # plt.ylim((3e33, 3e38))
# plt.xlabel(r'median GxBG ',fontsize=18)
# plt.ylabel(r'$N$ ',fontsize=18)
# plt.legend(fontsize=14)
# plt.tight_layout()
# plt.savefig(p_2_2PCF_figure)#, transparent = True)
# plt.clf()
# print(p_2_2PCF_figure, 'written')
# print(np.std(MED_WTH_GxB), np.std(MED_WTH_GxB[np.isin(p_2_all_srvval['GxB'], SRVMAP_exGAL_eroDE)]), np.std(MED_WTH_GxB[np.isin(p_2_all_srvval['GxB'], SRVMAP_exGAL_eroDE30)]))

p_2_all_srvval['GxB'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxB']])
Merge = get_merge_wth( p_2_all_xcorr['GxB'][np.isin(p_2_all_srvval['GxB'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALCENxEVT_wtheta_prediction_GxB.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)

basename = 'GALCEN_m'+m0_str+'_evAGN_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxA'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxA'].sort()
p_2_all_srvval['GxA'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxA']])
Merge = get_merge_wth( p_2_all_xcorr['GxA'][np.isin(p_2_all_srvval['GxA'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALCENxEVT_wtheta_prediction_GxA.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)

basename = 'GALCEN_m' + m0_str + '_evAGNevCLU_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxGA'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxGA'].sort()
p_2_all_srvval['GxGA'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxGA']])
Merge = get_merge_wth( p_2_all_xcorr['GxGA'][np.isin(p_2_all_srvval['GxGA'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALCENxEVT_wtheta_prediction_GxGA.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)

basename = 'GALCEN_m' + m0_str + '_evAGNevCLUevBKG_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxGAB'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxGAB'].sort()
p_2_all_srvval['GxGAB'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxGAB']])
Merge = get_merge_wth( p_2_all_xcorr['GxGAB'][np.isin(p_2_all_srvval['GxGAB'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALCENxEVT_wtheta_prediction_GxGAB.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)


m0 = 10.5
z1 = 0.26
m0_str = str(np.round(m0,1))
z1_str = str(np.round(z1,1))
p_2_all_xcorr = {}
p_2_all_srvval = {}

basename = 'GALCEN_m'+m0_str+'_evGAS_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxG'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxG'].sort()
p_2_all_srvval['GxG'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxG']])
Merge = get_merge_wth( p_2_all_xcorr['GxG'][np.isin(p_2_all_srvval['GxG'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_GALCENxEVT_wtheta_prediction_GxG.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)

basename = 'GALCEN_m'+m0_str+'_evBKG_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxB'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxB'].sort()
p_2_all_srvval['GxB'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxB']])
Merge = get_merge_wth( p_2_all_xcorr['GxB'][np.isin(p_2_all_srvval['GxB'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_GALCENxEVT_wtheta_prediction_GxB.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)

basename = 'GALCEN_m'+m0_str+'_evAGN_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxA'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxA'].sort()
p_2_all_srvval['GxA'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxA']])
Merge = get_merge_wth( p_2_all_xcorr['GxA'][np.isin(p_2_all_srvval['GxA'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_GALCENxEVT_wtheta_prediction_GxA.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)

basename = 'GALCEN_m' + m0_str + '_evAGNevCLU_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxGA'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxGA'].sort()
p_2_all_srvval['GxGA'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxGA']])
Merge = get_merge_wth( p_2_all_xcorr['GxGA'][np.isin(p_2_all_srvval['GxGA'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_GALCENxEVT_wtheta_prediction_GxGA.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)

basename = 'GALCEN_m' + m0_str + '_evAGNevCLUevBKG_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxGAB'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxGAB'].sort()
p_2_all_srvval['GxGAB'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxGAB']])
Merge = get_merge_wth( p_2_all_xcorr['GxGAB'][np.isin(p_2_all_srvval['GxGAB'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_GALCENxEVT_wtheta_prediction_GxGAB.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)



m0 = 11.0
z1 = 0.35
m0_str = str(np.round(m0,1))
z1_str = str(np.round(z1,1))
p_2_all_xcorr = {}
p_2_all_srvval = {}

basename = 'GALCEN_m'+m0_str+'_evGAS_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxG'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxG'].sort()
p_2_all_srvval['GxG'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxG']])
Merge = get_merge_wth( p_2_all_xcorr['GxG'][np.isin(p_2_all_srvval['GxG'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_GALCENxEVT_wtheta_prediction_GxG.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)

basename = 'GALCEN_m'+m0_str+'_evBKG_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxB'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxB'].sort()
p_2_all_srvval['GxB'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxB']])
Merge = get_merge_wth( p_2_all_xcorr['GxB'][np.isin(p_2_all_srvval['GxB'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_GALCENxEVT_wtheta_prediction_GxB.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)

basename = 'GALCEN_m'+m0_str+'_evAGN_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxA'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxA'].sort()
p_2_all_srvval['GxA'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxA']])
Merge = get_merge_wth( p_2_all_xcorr['GxA'][np.isin(p_2_all_srvval['GxA'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_GALCENxEVT_wtheta_prediction_GxA.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)

basename = 'GALCEN_m' + m0_str + '_evAGNevCLU_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxGA'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxGA'].sort()
p_2_all_srvval['GxGA'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxGA']])
Merge = get_merge_wth( p_2_all_xcorr['GxGA'][np.isin(p_2_all_srvval['GxGA'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_GALCENxEVT_wtheta_prediction_GxGA.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)

basename = 'GALCEN_m' + m0_str + '_evAGNevCLUevBKG_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxGAB'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxGAB'].sort()
p_2_all_srvval['GxGAB'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxGAB']])
Merge = get_merge_wth( p_2_all_xcorr['GxGAB'][np.isin(p_2_all_srvval['GxGAB'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_GALCENxEVT_wtheta_prediction_GxGAB.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)



m0 = 10.0
z1 = 0.18
m0_str = str(np.round(m0,1))
z1_str = str(np.round(z1,1))
p_2_all_xcorr = {}
p_2_all_srvval = {}

basename = 'GALSAT_m'+m0_str+'_evGAS_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxG'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxG'].sort()
p_2_all_srvval['GxG'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxG']])
Merge = get_merge_wth( p_2_all_xcorr['GxG'][np.isin(p_2_all_srvval['GxG'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALSATxEVT_wtheta_prediction_GxG.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)

basename = 'GALSAT_m'+m0_str+'_evBKG_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxB'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxB'].sort()
p_2_all_srvval['GxB'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxB']])
Merge = get_merge_wth( p_2_all_xcorr['GxB'][np.isin(p_2_all_srvval['GxB'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALSATxEVT_wtheta_prediction_GxB.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)

basename = 'GALSAT_m'+m0_str+'_evAGN_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxA'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxA'].sort()
p_2_all_srvval['GxA'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxA']])
Merge = get_merge_wth( p_2_all_xcorr['GxA'][np.isin(p_2_all_srvval['GxA'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALSATxEVT_wtheta_prediction_GxA.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)

basename = 'GALSAT_m' + m0_str + '_evAGNevCLU_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxGA'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxGA'].sort()
p_2_all_srvval['GxGA'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxGA']])
Merge = get_merge_wth( p_2_all_xcorr['GxGA'][np.isin(p_2_all_srvval['GxGA'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALSATxEVT_wtheta_prediction_GxGA.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)

basename = 'GALSAT_m' + m0_str + '_evAGNevCLUevBKG_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxGAB'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxGAB'].sort()
p_2_all_srvval['GxGAB'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxGAB']])
Merge = get_merge_wth( p_2_all_xcorr['GxGAB'][np.isin(p_2_all_srvval['GxGAB'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALSATxEVT_wtheta_prediction_GxGAB.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)


m0 = 10.5
z1 = 0.26
m0_str = str(np.round(m0,1))
z1_str = str(np.round(z1,1))
p_2_all_xcorr = {}
p_2_all_srvval = {}

basename = 'GALSAT_m'+m0_str+'_evGAS_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxG'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxG'].sort()
p_2_all_srvval['GxG'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxG']])
Merge = get_merge_wth( p_2_all_xcorr['GxG'][np.isin(p_2_all_srvval['GxG'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_GALSATxEVT_wtheta_prediction_GxG.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)

basename = 'GALSAT_m'+m0_str+'_evBKG_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxB'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxB'].sort()
p_2_all_srvval['GxB'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxB']])
Merge = get_merge_wth( p_2_all_xcorr['GxB'][np.isin(p_2_all_srvval['GxB'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_GALSATxEVT_wtheta_prediction_GxB.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)

basename = 'GALSAT_m'+m0_str+'_evAGN_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxA'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxA'].sort()
p_2_all_srvval['GxA'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxA']])
Merge = get_merge_wth( p_2_all_xcorr['GxA'][np.isin(p_2_all_srvval['GxA'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_GALSATxEVT_wtheta_prediction_GxA.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)

basename = 'GALSAT_m' + m0_str + '_evAGNevCLU_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxGA'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxGA'].sort()
p_2_all_srvval['GxGA'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxGA']])
Merge = get_merge_wth( p_2_all_xcorr['GxGA'][np.isin(p_2_all_srvval['GxGA'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_GALSATxEVT_wtheta_prediction_GxGA.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)

basename = 'GALSAT_m' + m0_str + '_evAGNevCLUevBKG_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxGAB'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxGAB'].sort()
p_2_all_srvval['GxGAB'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxGAB']])
Merge = get_merge_wth( p_2_all_xcorr['GxGAB'][np.isin(p_2_all_srvval['GxGAB'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_GALSATxEVT_wtheta_prediction_GxGAB.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)



m0 = 11.0
z1 = 0.35
m0_str = str(np.round(m0,1))
z1_str = str(np.round(z1,1))
p_2_all_xcorr = {}
p_2_all_srvval = {}

basename = 'GALSAT_m'+m0_str+'_evGAS_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxG'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxG'].sort()
p_2_all_srvval['GxG'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxG']])
Merge = get_merge_wth( p_2_all_xcorr['GxG'][np.isin(p_2_all_srvval['GxG'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_GALSATxEVT_wtheta_prediction_GxG.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)

basename = 'GALSAT_m'+m0_str+'_evBKG_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxB'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxB'].sort()
p_2_all_srvval['GxB'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxB']])
Merge = get_merge_wth( p_2_all_xcorr['GxB'][np.isin(p_2_all_srvval['GxB'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_GALSATxEVT_wtheta_prediction_GxB.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)

basename = 'GALSAT_m'+m0_str+'_evAGN_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxA'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxA'].sort()
p_2_all_srvval['GxA'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxA']])
Merge = get_merge_wth( p_2_all_xcorr['GxA'][np.isin(p_2_all_srvval['GxA'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_GALSATxEVT_wtheta_prediction_GxA.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)

basename = 'GALSAT_m' + m0_str + '_evAGNevCLU_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxGA'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxGA'].sort()
p_2_all_srvval['GxGA'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxGA']])
Merge = get_merge_wth( p_2_all_xcorr['GxGA'][np.isin(p_2_all_srvval['GxGA'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_GALSATxEVT_wtheta_prediction_GxGA.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)

basename = 'GALSAT_m' + m0_str + '_evAGNevCLUevBKG_CROSSCORR_05E20.wtheta.2pcf.fits'
p_2_all_xcorr['GxGAB'] = np.array( glob.glob( os.path.join(dir_2pcf, basename  ) ) )
p_2_all_xcorr['GxGAB'].sort()
p_2_all_srvval['GxGAB'] = np.array([int(el.split('/')[7]) for el in p_2_all_xcorr['GxGAB']])
Merge = get_merge_wth( p_2_all_xcorr['GxGAB'][np.isin(p_2_all_srvval['GxGAB'], SRVMAP_exGAL_eroDE)] )
p_2_2PCF_GALxEVTc030singleRRDR10 = os.path.join(benchmark_dir, 'LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_GALSATxEVT_wtheta_prediction_GxGAB.fits')
Merge.write(p_2_2PCF_GALxEVTc030singleRRDR10, overwrite = True)
print(p_2_2PCF_GALxEVTc030singleRRDR10)




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

