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

nl = lambda sel : len(sel.nonzero()[0])

sky_map_hdu = Table.read(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS.fits') )

LC_dir = 'LCerass'
top_dir = os.path.join(os.environ['UCHUU'], LC_dir)


for sky_tile in sky_map_hdu[(sky_map_hdu['OWNER']==2)|(sky_map_hdu['OWNER']==0)][:1]:#[604:605] :

	sky_tile_id = str(sky_tile['SRVMAP'])
	str_field = str(sky_tile['SRVMAP']).zfill(6)

	agn_seed = 1
	agn_dir = os.path.join(os.environ['UCHUU'], LC_dir, str_field, "eRASS8_SEED_"+str(agn_seed).zfill(3) +"_events_AGN_2025_04" )

	clu_seed = 1
	cluster_dir = os.path.join(os.environ['UCHUU'], LC_dir, str_field, "Att_eRASS8_sixte_v27_SEED_"+str(clu_seed).zfill(3) +"_events_cluster" )

	stars_dir   = os.path.join( os.environ['UCHUU'], LC_dir, str_field, 'stars')#, 'simulated_photons_ccd?.fits' )
	bg_dir      = os.path.join( os.environ['UCHUU'], LC_dir, str_field, 'pBG' ) # 'evt_particle_???.fits' )

	RA  = str(sky_tile['RA_CEN'])
	DEC = str(sky_tile['DE_CEN'])

	log_dir = os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'logs-erass8')
	os.system('mkdir -p '+log_dir)

	esass_dir = os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'sim_evt_AGN_CLU')

	Attitude_File = "/home/idies/workspace/erosim/erosita_attitude/eRASS_4yr_epc85_att.fits"

	path_2_tmp_file = os.path.join(esass_dir, 'tmp_'+str_field+'.fits')
	path_2_tmp2_file = os.path.join(esass_dir, 'tmp2_'+str_field+'.fits')
	path_2_event_file = os.path.join(esass_dir, 'evt_'+str_field+'.fits')

	# concatenates events
	path_2_list = os.path.join(esass_dir, 'FTCOPY.list')
	cmd_list = "ls " + os.path.join(esass_dir, "*_FTCOPY.fits") + " > " + path_2_list
	print( cmd_list )
	os.system( cmd_list )

	# evtool command (requires eSASS to be sourced)
	# evtool eventfiles="@tmp_${hid}.list" outfile=tmp_${hid}.fits pattern=15 clobber=yes
	cmd_evtool = """evtool eventfiles='@"""+path_2_list+"""' outfile="""+ path_2_tmp_file+" pattern=15 clobber=yes "
	#cmd_evtool_erass1 = """evtool eventfiles='@"""+path_2_list+"""' GTI="617943605 649479605" outfile="""+ path_2_tmp_file+" pattern=15 clobber=yes "
	print(cmd_evtool)
	os.system(cmd_evtool)

	# gives the tile center position
	# radec2xy file=tmp_${hid}.fits ra0=${ra_cen} dec0=${dec_cen}
	cmd_radec = "radec2xy file="+path_2_tmp_file+" ra0="+RA+" dec0="+DEC
	print(cmd_radec)
	os.system(cmd_radec)

	# ftcopy tmp_${hid}.fits"[(RAWX-192.5)**2+(RAWY-192.5)**2<=191.5**2]" evt_${hid}.fits clobber=yes
	cmd_copy = "ftcopy "+path_2_tmp_file+""""[(RAWX-192.5)**2+(RAWY-192.5)**2<=191.5**2]" """+path_2_event_file+" clobber=yes"
	print(cmd_copy)
	os.system(cmd_copy)

	# fparkey SIXTE evt_${hid}.fits[1] CREATOR add=yes
	cmd_parkey = "fparkey SIXTE "+path_2_event_file+"[1] CREATOR add=yes"
	print(cmd_parkey)
	os.system(cmd_parkey)

	## deleting unused files

	path_rm_list = os.path.join(esass_dir, 'RM.list')

	to_remove = "ls " + os.path.join(esass_dir, "?_t*erass_ccd?_evt.fits") + " > " + path_rm_list
	print( to_remove )
	os.system( to_remove )
	cmd_list_rm = n.loadtxt(path_rm_list, dtype='str')
	for element in cmd_list_rm:
		os.remove(element)

	os.remove(path_rm_list)

	to_remove = "ls " + os.path.join(esass_dir, "?_t*erass_ccd?_evt_FTCOPY.fits") + " > " + path_rm_list
	print( to_remove )
	os.system( to_remove )
	cmd_list_rm = n.loadtxt(path_rm_list, dtype='str')
	for element in cmd_list_rm:
		os.remove(element)


	os.remove(path_rm_list)
	os.remove(path_2_list)

	os.remove(path_2_tmp_file)

