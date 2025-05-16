import os, sys, glob
import numpy as n
import astropy.io.fits as fits
#import healpy
import time
t0 = time.time()
from astropy.table import Table, vstack
import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as plt

clu_seed =  sys.argv[1] # 1

nl = lambda sel : len(sel.nonzero()[0])

sky_map_hdu = Table.read(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS.fits') )

LC_dir = 'LCerass'
top_dir = os.path.join(os.environ['UCHUU'], LC_dir)

def get_srvmap(ra, dec):
    return sky_map_hdu['SRVMAP'].value[(sky_map_hdu['RA_MIN']<ra ) & ( sky_map_hdu['RA_MAX'] >= ra ) & ( sky_map_hdu['DE_MIN']<dec ) & ( sky_map_hdu['DE_MAX'] >= dec)]

for sky_tile in sky_map_hdu[(sky_map_hdu['OWNER']==2)|(sky_map_hdu['OWNER']==0)]:
	print(sky_tile)
	sky_tile_id = str(sky_tile['SRVMAP'])
	str_field = str(sky_tile['SRVMAP']).zfill(6)
	evt_list = np.array(glob.glob(os.path.join(os.environ['UCHUU'], LC_dir, str_field, 's4_c030', '*_Image_c030.fits.gz' ) ) )
	log_dir = os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'logs-erass8')
	os.system('mkdir -p '+log_dir)
	#esass_dir = os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'sim_evt_e4_merge')
	esass_dir = os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'GE_e4_merge_SimBKG_CLUseed'+clu_seed.zfill(3))
	os.system('mkdir -p '+esass_dir)

	path_2_event_file = os.path.join(esass_dir, 'evt_'+str_field+'.fits')
	path_2_simeventCLU_file = os.path.join(esass_dir, 'simCLUevt_'+str_field+'.fits')
	path_2_simeventBKG_file = os.path.join(esass_dir, 'simBKGevt_'+str_field+'.fits')
	if len(evt_list)==0 or os.path.isfile(path_2_event_file):
		print('continue', len(evt_list)==0, os.path.isfile(path_2_event_file))
		continue
	bg_dir      = os.path.join( os.environ['UCHUU'], LC_dir, str_field, 'pBG2' ) # 'evt_particle_???.fits' )
	BG_evt_files = n.array( glob.glob( os.path.join( bg_dir, '*.fits' ) ) )
	if len(BG_evt_files)==0:
		print('continue', len(BG_evt_files), 'no BG files')
		continue
	hdul_raw = fits.open(evt_list[0])
	texps = np.array([ np.sum(hdul_raw['GTI1'].data['STOP']-hdul_raw['GTI1'].data['START'])
			, np.sum(hdul_raw['GTI2'].data['STOP']-hdul_raw['GTI2'].data['START'])
			, np.sum(hdul_raw['GTI3'].data['STOP']-hdul_raw['GTI3'].data['START'])
			, np.sum(hdul_raw['GTI4'].data['STOP']-hdul_raw['GTI4'].data['START'])
			, np.sum(hdul_raw['GTI5'].data['STOP']-hdul_raw['GTI5'].data['START'])
			, np.sum(hdul_raw['GTI6'].data['STOP']-hdul_raw['GTI6'].data['START'])
			, np.sum(hdul_raw['GTI7'].data['STOP']-hdul_raw['GTI7'].data['START']) ])
	N_ev_OBS = len(hdul_raw['EVENTS'].data)
	hdul = fits.open(evt_list[0])
	hdul['EVENTS'].data['RA'][hdul['EVENTS'].data['RA']==0]=1e-6
	SRV_ev = np.array([get_srvmap(e0,e1) for e0,e1 in zip(hdul['EVENTS'].data['RA'], hdul['EVENTS'].data['DEC']) ])
	to_replace = ( np.hstack((SRV_ev)) == sky_tile['SRVMAP'] )
	ids_to_replace = np.arange(len(to_replace))[to_replace]
	N_ev_OBS = len(ids_to_replace)
	cluster_dir = os.path.join(os.environ['UCHUU'], LC_dir, str_field, "Att_eRASS8_sixte_v27_SEED_"+str(clu_seed).zfill(3) +"_events_cluster" )

	bg_all = []
	for el in BG_evt_files:
		tt0 = Table.read(el)
		tt0['TM_NR'] = int(os.path.basename(el).split('_')[-2][-1])
		tt0.keep_columns(['RA', 'DEC','RAWX', 'RAWY', 'PHA', 'SIGNAL', 'TM_NR'])
		bg_all.append(tt0)
	bg_all = vstack((bg_all))

	N_evs = []
	for NCCD, tEXP in zip(n.arange(7)+1, texps):
		CL_evt_files = n.array( glob.glob( os.path.join( cluster_dir, 't0erass_ccd' + str(NCCD) + '_evt.fits' ) ) )
		hdu_C = fits.open(CL_evt_files[0])
		texp_C = np.sum(hdu_C[2].data['STOP']-hdu_C[2].data['START'])
		frac_C = tEXP/texp_C
		N_ev_C = len(hdu_C[1].data)

		bg_tm = bg_all[bg_all['TM_NR']==NCCD]
		N_BG = len(bg_tm)
		N_evs.append([ N_ev_C, N_BG])

	f_BKG = 0.98
	f_CLU = 0.02
	NC, NB = np.transpose(N_evs).sum(axis=1)
	N_ev_B = int(1. * N_ev_OBS/7)+1
	N_ev_C = int(f_CLU * N_ev_OBS/7)+1

	if N_ev_B*7>len(bg_all):
		print('continue', 'not enough BG events', len(bg_all), 'when ', N_ev_B*7, 'are needed')
		continue

	data_C = []
	data_B = []

	for NCCD, tEXP in zip(n.arange(7)+1, texps):
		CL_evt_files = n.array( glob.glob( os.path.join( cluster_dir, 't0erass_ccd' + str(NCCD) + '_evt.fits' ) ) )

		hdu_C = fits.open(CL_evt_files[0])
		if N_ev_C<=len(hdu_C[1].data):
			id_C = np.random.choice(np.arange(len(hdu_C[1].data)), size = N_ev_C, replace = False)
		else:
			id_C = np.arange(len(hdu_C[1].data))
		data_C.append( Table(hdu_C['EVENTS'].data[id_C]) )

		bg_tm = bg_all[bg_all['TM_NR']==NCCD]
		if N_ev_B<len(bg_tm):
			id_B = np.random.choice(np.arange(len(bg_tm)), size = N_ev_B, replace = False)
			data_B.append( bg_tm[id_B] )
		else:
			print('continue', 'not enough BG events', len(bg_tm), 'when ', N_ev_B, 'are needed')
			continue


	data_C = vstack((data_C))
	data_B = vstack((data_B))
	data_B = data_B[:N_ev_OBS-len(data_C)]

	fi_up = ['RA', 'DEC','RAWX', 'RAWY', 'PHA']#, 'X', 'Y']
	for fn in fi_up:
		hdul['EVENTS'].data[fn][ids_to_replace] = np.hstack((data_C[fn], data_B[fn]))[:N_ev_OBS]

	fn='SIGNAL'
	hdul['EVENTS'].data['PI'][ids_to_replace] = 1000.*np.hstack((data_C[fn], data_B[fn]))[:N_ev_OBS]

	hdul.writeto(path_2_event_file, overwrite=True)
	print(path_2_event_file)
	print('='*100)
	print('='*100)
	data_C.write(path_2_simeventCLU_file, overwrite = True)
	data_B.write(path_2_simeventBKG_file, overwrite = True)
	print(path_2_simeventCLU_file, len(data_C))
	print(path_2_simeventBKG_file, len(data_B))
	N_sim =  len(data_C) + len(data_B)
	print(f_BKG, np.round(len(data_B)/N_sim,4))
	print(f_CLU, np.round(len(data_C)/N_sim,4))




