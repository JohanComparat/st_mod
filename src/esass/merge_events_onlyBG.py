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
	esass_dir = os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'GE_e4_merge_SimBKG')
	os.system('mkdir -p '+esass_dir)

	path_2_event_file = os.path.join(esass_dir, 'evt_'+str_field+'.fits')
	path_2_simeventAGN_file = os.path.join(esass_dir, 'simAGNevt_'+str_field+'.fits')
	path_2_simeventCLU_file = os.path.join(esass_dir, 'simCLUevt_'+str_field+'.fits')
	path_2_simeventBKG_file = os.path.join(esass_dir, 'simBKGevt_'+str_field+'.fits')
	if len(evt_list)==0 or os.path.isfile(path_2_event_file):
		print('continue', len(evt_list)==0, os.path.isfile(path_2_event_file))
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
	SRV_ev = np.array([get_srvmap(e0,e1) for e0,e1 in zip(hdul['EVENTS'].data['RA'], hdul['EVENTS'].data['DEC']) ])
	to_replace = ( np.hstack((SRV_ev)) == sky_tile['SRVMAP'] )
	ids_to_replace = np.arange(len(to_replace))[to_replace]
	N_ev_OBS = len(ids_to_replace)
	bg_dir      = os.path.join( os.environ['UCHUU'], LC_dir, str_field, 'pBG' ) # 'evt_particle_???.fits' )

	BG_evt_files = n.array( glob.glob( os.path.join( bg_dir, 'evt_particle_???.fits' ) ) )
	bg_all = vstack(([Table.read(el) for el in BG_evt_files]))

	N_ev_B = int(N_ev_OBS/7)+5
	data_B = []

	for NCCD, tEXP in zip(n.arange(7)+1, texps):
		bg_tm = bg_all[bg_all['TM_NR']==NCCD]
		if N_ev_B<len(bg_tm):
			id_B = np.random.choice(np.arange(len(bg_tm)), size = N_ev_B, replace = False)
			data_B.append( bg_tm[id_B] )
		else:
			print('continue', 'not enough BG events', len(bg_tm), 'when ', N_ev_B, 'are needed')
			continue

	data_B = vstack((data_B))
	data_B = data_B[:N_ev_OBS]
	#data_S = vstack((data_S))

	fi_up = ['RA', 'DEC','RAWX', 'RAWY', 'PHA']#, 'X', 'Y']
	for fn in fi_up:
		hdul['EVENTS'].data[fn][ids_to_replace] = data_B[fn]

	hdul['EVENTS'].data['PI'][ids_to_replace] = data_B['ENERGY']

	hdul.writeto(path_2_event_file, overwrite=True)
	print(path_2_event_file)
	print('='*100)
	print('='*100)
	data_B.write(path_2_simeventBKG_file, overwrite = True)
	print(path_2_simeventBKG_file, len(data_B))


