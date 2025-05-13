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
agn_seed = sys.argv[1] # 1
clu_seed = sys.argv[2] # 1

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
	esass_dir = os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'GE_e4_merge_AGNseed'+agn_seed.zfill(3)+'_SimBKG_CLUseed'+clu_seed.zfill(3))
	os.system('mkdir -p '+esass_dir)

	path_2_event_file = os.path.join(esass_dir, 'evt_'+str_field+'.fits')
	path_2_simeventAGN_file = os.path.join(esass_dir, 'simAGNevt_'+str_field+'.fits')
	path_2_simeventCLU_file = os.path.join(esass_dir, 'simCLUevt_'+str_field+'.fits')
	path_2_simeventBKG_file = os.path.join(esass_dir, 'simBKGevt_'+str_field+'.fits')
	if len(evt_list)==0 or os.path.isfile(path_2_event_file):
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
	agn_dir = os.path.join(os.environ['UCHUU'], LC_dir, str_field, "eRASS8_SEED_"+str(agn_seed).zfill(3) +"_events_AGN_2025_04" )
	cluster_dir = os.path.join(os.environ['UCHUU'], LC_dir, str_field, "Att_eRASS8_sixte_v27_SEED_"+str(clu_seed).zfill(3) +"_events_cluster" )
	stars_dir   = os.path.join( os.environ['UCHUU'], LC_dir, str_field, 'stars')#, 'simulated_photons_ccd?.fits' )
	bg_dir      = os.path.join( os.environ['UCHUU'], LC_dir, str_field, 'pBG' ) # 'evt_particle_???.fits' )
	#RA  = str(sky_tile['RA_CEN'])
	#DEC = str(sky_tile['DE_CEN'])

	BG_evt_files = n.array( glob.glob( os.path.join( bg_dir, 'evt_particle_???.fits' ) ) )
	bg_all = vstack(([Table.read(el) for el in BG_evt_files]))
	bg_all['TIME'].min()

	N_evs = []
	for NCCD, tEXP in zip(n.arange(7)+1, texps):
		agn_evt_files = n.array( glob.glob( os.path.join( agn_dir, 't0erass_ccd' + str(NCCD) + '_evt.fits' ) ) )
		CL_evt_files = n.array( glob.glob( os.path.join( cluster_dir, 't0erass_ccd' + str(NCCD) + '_evt.fits' ) ) )
		ST_evt_files = n.array( glob.glob( os.path.join( stars_dir, 'simulated_photons_ccd' + str(NCCD) + '.fits' ) ) )

		hdu_A = fits.open(agn_evt_files[0])
		texp_A = np.sum(hdu_A[2].data['STOP']-hdu_A[2].data['START'])
		frac_A = tEXP/texp_A
		N_ev_A = len(hdu_A[1].data)

		hdu_C = fits.open(CL_evt_files[0])
		texp_C = np.sum(hdu_C[2].data['STOP']-hdu_C[2].data['START'])
		frac_C = tEXP/texp_C
		N_ev_C = len(hdu_C[1].data)

		#hdu_S = fits.open(ST_evt_files[0])
		#N_ev_S = len(hdu_S[1].data)

		bg_tm = bg_all[bg_all['TM_NR']==NCCD]
		N_BG = len(bg_tm)
		N_evs.append([N_ev_A, N_ev_C, N_BG])

	f_AGN = 0.123
	f_BKG = 0.858
	f_CLU = 0.019
	NA, NC, NB = np.transpose(N_evs).sum(axis=1)
	N_ev_A = int(f_AGN * N_ev_OBS/7)+1
	N_ev_B = int(0.95 * N_ev_OBS/7)+1
	N_ev_C = int(f_CLU * N_ev_OBS/7)+1
	#frac_all = N_ev_OBS / np.sum(N_evs)+0.01

	data_A = []
	data_C = []
	#data_S = []
	data_B = []

	for NCCD, tEXP in zip(n.arange(7)+1, texps):
		agn_evt_files = n.array( glob.glob( os.path.join( agn_dir, 't0erass_ccd' + str(NCCD) + '_evt.fits' ) ) )
		CL_evt_files = n.array( glob.glob( os.path.join( cluster_dir, 't0erass_ccd' + str(NCCD) + '_evt.fits' ) ) )
		ST_evt_files = n.array( glob.glob( os.path.join( stars_dir, 'simulated_photons_ccd' + str(NCCD) + '.fits' ) ) )
		hdu_A = fits.open(agn_evt_files[0])
		#N_ev_A = int(len(hdu_A[1].data) * frac_all) + 20
		id_A = np.random.choice(np.arange(len(hdu_A[1].data)), size = N_ev_A, replace = False)
		data_A.append( Table(hdu_A['EVENTS'].data[id_A]) )

		hdu_C = fits.open(CL_evt_files[0])
		#N_ev_C = int(len(hdu_C[1].data) * frac_all) + 20
		if N_ev_C<=len(hdu_C[1].data):
			id_C = np.random.choice(np.arange(len(hdu_C[1].data)), size = N_ev_C, replace = False)
		else:
			id_C = np.arange(len(hdu_C[1].data))
		data_C.append( Table(hdu_C['EVENTS'].data[id_C]) )

		#hdu_S = fits.open(ST_evt_files[0])
		#N_ev_S = int(len(hdu_S[1].data) * frac_all) +1
		#id_S = np.random.choice(np.arange(len(hdu_S[1].data)), size = N_ev_S, replace = False)
		#data_S.append( Table(hdu_S[1].data[id_S]) )

		bg_tm = bg_all[bg_all['TM_NR']==NCCD]
		#N_ev_B = int(len(bg_tm) * frac_all) + 100
		if N_ev_B<len(bg_tm):
			id_B = np.random.choice(np.arange(len(bg_tm)), size = N_ev_B, replace = False)
			data_B.append( bg_tm[id_B] )
		else:
			print('continue', 'not enough BG events', len(bg_tm), 'when ', N_ev_B, 'are needed')
			continue

	data_A = vstack((data_A))
	data_C = vstack((data_C))
	data_B = vstack((data_B))
	data_B = data_B[:N_ev_OBS-(len(data_A)-len(data_C))]
	#data_S = vstack((data_S))

	fi_up = ['RA', 'DEC','RAWX', 'RAWY', 'PHA']#, 'X', 'Y']
	for fn in fi_up:
		hdul['EVENTS'].data[fn][ids_to_replace] = np.hstack((data_C[fn], data_A[fn], data_B[fn]))[:N_ev_OBS]

	fn='SIGNAL'
	hdul['EVENTS'].data['PI'][ids_to_replace] = 1000.*np.hstack((data_C[fn], data_A[fn], data_B['ENERGY']/1000.))[:N_ev_OBS]

	hdul.writeto(path_2_event_file, overwrite=True)
	print(path_2_event_file)
	print('='*100)
	print('='*100)
	data_A.write(path_2_simeventAGN_file, overwrite = True)
	data_C.write(path_2_simeventCLU_file, overwrite = True)
	data_B.write(path_2_simeventBKG_file, overwrite = True)
	print(path_2_simeventAGN_file, len(data_A))
	print(path_2_simeventCLU_file, len(data_C))
	print(path_2_simeventBKG_file, len(data_B))
	N_sim = len(data_A) + len(data_C) + len(data_B)
	print(f_AGN, np.round(len(data_A)/N_sim,4))
	print(f_BKG, np.round(len(data_B)/N_sim,4))
	print(f_CLU, np.round(len(data_C)/N_sim,4))

	#fig_out = os.path.join('test-ra-dec.png' )
	#plt.figure(0, (6.,6.))
	#plt.plot(hdul['EVENTS'].data['RA'], hdul['EVENTS'].data['DEC'], 'k,')
	##plt.plot(data_A['RA'], data_A['DEC'], 'b,')
	##plt.plot(data_C['RA'], data_C['DEC'], 'r,')
	##plt.plot(data_S['RA'], data_S['DEC'], 'g,')
	##plt.plot(data_B['RA'], data_B['DEC'], 'b,')
	#plt.savefig(fig_out)
	##plt.clf()
	##fig_out = os.path.join('test-ra-dec-BG.png' )
	##plt.figure(0, (6.,6.))
	##plt.plot(hdul['EVENTS'].data['RA'], hdul['EVENTS'].data['DEC'], 'k,')
	##plt.plot(data_B['RA'], data_B['DEC'], 'b,')
	##plt.savefig(fig_out)
	##plt.clf()


### eRASS:4 take half of the eRASS8 events

##In [61]: hdu_C[1].data.columns
##Out[61]:
##ColDefs(
    ##name = 'TIME'; format = 'D'; unit = 's'
    ##name = 'FRAME'; format = 'J'
    ##name = 'PHA'; format = 'J'; unit = 'ADU'
    ##name = 'PI'; format = 'J'; unit = 'ADU'
    ##name = 'SIGNAL'; format = 'E'; unit = 'keV'
    ##name = 'RAWX'; format = 'I'; unit = 'pixel'
    ##name = 'RAWY'; format = 'I'; unit = 'pixel'
    ##name = 'RA'; format = 'D'; unit = 'deg'
    ##name = 'DEC'; format = 'D'; unit = 'deg'
    ##name = 'PH_ID'; format = '2J'
    ##name = 'SRC_ID'; format = '2J'
    ##name = 'TYPE'; format = 'I'
    ##name = 'NPIXELS'; format = 'J'
    ##name = 'PILEUP'; format = 'I'
    ##name = 'SIGNALS'; format = '9E'; unit = 'keV'
    ##name = 'PHAS'; format = '9J'; unit = 'ADU'
##)

## Stars have a stupid format and lack information ... needs to be re-generated
##ColDefs(
    ##name = 'RA'; format = 'D'; unit = 'deg'
    ##name = 'DEC'; format = 'D'; unit = 'deg'
    ##name = 'SRC_ID_1'; format = 'J'
    ##name = 'TIME'; format = 'D'; unit = 's'
    ##name = 'SIGNAL'; format = 'E'; unit = 'keV'
    ##name = 'SRVMAP'; format = 'J'


##In [72]: bg_all.info()
##<Table length=473336>
    ##name     dtype  unit
##----------- ------- ----
       ##TIME float64    s
         ##RA float64  deg
        ##DEC float64  deg
          ##X float64
          ##Y float64
     ##ENERGY float32
  ##EV_WEIGHT float32
       ##RAWX   int16
       ##RAWY   int16
       ##SUBX float64
       ##SUBY float64
        ##PHA   int16
    ##PAT_TYP   int16
    ##PAT_INF   uint8
      ##TM_NR   uint8
       ##FLAG   int32
  ##FRAMETIME float64    s
 ##RECORDTIME float64
##HEALPIX_VAL   int64  deg
     ##SRVMAP   int32

##In [62]: hdul['EVENTS'].columns
##Out[62]:
##ColDefs(
    ##name = 'TIME'; format = 'D'; unit = 's'
    ##name = 'RA'; format = 'D'; unit = 'deg'; bscale = 1e-06; bzero = 0.0
    ##name = 'DEC'; format = 'D'; unit = 'deg'; bscale = 1e-06; bzero = 0.0
    ##name = 'X'; format = 'D'; coord_type = 'RA---SIN'; coord_unit = 'deg'; coord_ref_point = 0.0; coord_ref_value = 120.659341; coord_inc = -1.3888889095849e-05
    ##name = 'Y'; format = 'D'; coord_type = 'DEC--SIN'; coord_unit = 'deg'; coord_ref_point = 0.0; coord_ref_value = 42.008289; coord_inc = 1.3888889095849e-05
    ##name = 'PI'; format = 'E'
    ##name = 'EV_WEIGHT'; format = 'E'
    ##name = 'RAWX'; format = 'I'
    ##name = 'RAWY'; format = 'I'
    ##name = 'SUBX'; format = 'D'; bscale = 0.00666666667; bzero = -0.843333333
    ##name = 'SUBY'; format = 'D'; bscale = 0.00666666667; bzero = -0.843333333
    ##name = 'PHA'; format = 'I'
    ##name = 'PAT_TYP'; format = 'I'
    ##name = 'PAT_INF'; format = 'B'
    ##name = 'TM_NR'; format = 'B'
    ##name = 'FLAG'; format = 'J'
    ##name = 'FRAMETIME'; format = 'D'; unit = 's'
    ##name = 'RECORDTIME'; format = 'D'
##)




