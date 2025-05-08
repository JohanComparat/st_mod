"""

Pre-rpocess sixte event files to input eSASS

Adapted from bash to python

Based on the Script from Teng Liu

"""
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

agn_seed = sys.argv[1] # 1

sky_map_hdu = Table.read(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS.fits') )

LC_dir = 'LCerass'
top_dir = os.path.join(os.environ['UCHUU'], LC_dir)

def get_srvmap(ra, dec):
    return sky_map_hdu['SRVMAP'].value[(sky_map_hdu['RA_MIN']<ra ) & ( sky_map_hdu['RA_MAX'] >= ra ) & ( sky_map_hdu['DE_MIN']<dec ) & ( sky_map_hdu['DE_MAX'] >= dec)]

for sky_tile in sky_map_hdu[(sky_map_hdu['OWNER']==2)|(sky_map_hdu['OWNER']==0)][:1]:
	print(sky_tile)
	sky_tile_id = str(sky_tile['SRVMAP'])
	str_field = str(sky_tile['SRVMAP']).zfill(6)
	evt_list = np.array(glob.glob(os.path.join(os.environ['UCHUU'], LC_dir, str_field, 's4_c030', '*_Image_c030.fits.gz' ) ) )
	log_dir = os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'logs-erass8')
	os.system('mkdir -p '+log_dir)
	esass_dir = os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'GE_e4_merge_AGNseed'+agn_seed.zfill(3)+'_SimBKG')
	os.system('mkdir -p '+esass_dir)

	path_2_event_file = os.path.join(esass_dir, 'evt_'+str_field+'.fits')
	path_2_simeventAGN_file = os.path.join(esass_dir, 'simAGNevt_'+str_field+'.fits')
	path_2_simeventBKG_file = os.path.join(esass_dir, 'simBKGevt_'+str_field+'.fits')
	#if len(evt_list)==0 or os.path.isfile(path_2_event_file):
		#continue
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
	agn_dir = os.path.join(os.environ['UCHUU'], LC_dir, str_field, "eRASS8_SEED_"+agn_seed.zfill(3) +"_events_AGN_2025_04" )
	bg_dir      = os.path.join( os.environ['UCHUU'], LC_dir, str_field, 'pBG' )

	BG_evt_files = n.array( glob.glob( os.path.join( bg_dir, 'evt_particle_???.fits' ) ) )
	bg_all = vstack(([Table.read(el) for el in BG_evt_files]))
	bg_all['TIME'].min()

	N_evs = []
	for NCCD, tEXP in zip(n.arange(7)+1, texps):
		agn_evt_files = n.array( glob.glob( os.path.join( agn_dir, 't0erass_ccd' + str(NCCD) + '_evt.fits' ) ) )

		hdu_A = fits.open(agn_evt_files[0])
		texp_A = np.sum(hdu_A[2].data['STOP']-hdu_A[2].data['START'])
		frac_A = tEXP/texp_A
		N_ev_A = len(hdu_A[1].data)

		bg_tm = bg_all[bg_all['TM_NR']==NCCD]
		N_BG = len(bg_tm)
		N_evs.append([N_ev_A, N_BG])

	f_AGN = 0.14
	f_BKG = 0.90
	NA, NB = np.transpose(N_evs).sum(axis=1)
	N_ev_A = int(f_AGN * N_ev_OBS/7)+1
	N_ev_B = int(0.95 * N_ev_OBS/7)+1

	data_A = []
	data_B = []

	for NCCD, tEXP in zip(n.arange(7)+1, texps):
		agn_evt_files = n.array( glob.glob( os.path.join( agn_dir, 't0erass_ccd' + str(NCCD) + '_evt.fits' ) ) )

		hdu_A = fits.open(agn_evt_files[0])
		id_A = np.random.choice(np.arange(len(hdu_A[1].data)), size = N_ev_A, replace = False)
		data_A.append( Table(hdu_A['EVENTS'].data[id_A]) )

		bg_tm = bg_all[bg_all['TM_NR']==NCCD]
		id_B = np.random.choice(np.arange(len(bg_tm)), size = N_ev_B, replace = False)
		data_B.append( bg_tm[id_B] )


	data_A = vstack((data_A))
	data_B = vstack((data_B))
	data_B = data_B[:N_ev_OBS-len(data_A)]

	fi_up = ['RA', 'DEC','RAWX', 'RAWY', 'PHA']#, 'X', 'Y']
	for fn in fi_up:
		hdul['EVENTS'].data[fn][ids_to_replace] = np.hstack((data_A[fn], data_B[fn]))[:N_ev_OBS]

	fn='SIGNAL'
	hdul['EVENTS'].data['PI'][ids_to_replace] = 1000.*np.hstack((data_A[fn], data_B['ENERGY']/1000.))[:N_ev_OBS]

	hdul.writeto(path_2_event_file, overwrite=True)
	print(path_2_event_file)
	print('='*100)
	print('='*100)
	data_A.write(path_2_simeventAGN_file, overwrite = True)
	data_B.write(path_2_simeventBKG_file, overwrite = True)
	print(path_2_simeventAGN_file, len(data_A))
	print(path_2_simeventBKG_file, len(data_B))
	N_sim = len(data_A) + len(data_B)
	print(f_AGN, np.round(len(data_A)/N_sim,4))
	print(f_BKG, np.round(len(data_B)/N_sim,4))



