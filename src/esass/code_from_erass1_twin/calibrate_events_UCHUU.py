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


for sky_tile in sky_map_hdu[(sky_map_hdu['OWNER']==2)|(sky_map_hdu['OWNER']==0)]:

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
	os.system('mkdir -p '+esass_dir)

	Attitude_File = "/home/idies/workspace/erosim/erosita_attitude/eRASS_4yr_epc85_att.fits"

	path_2_tmp_file = os.path.join(esass_dir, 'tmp_'+str_field+'.fits')
	path_2_tmp2_file = os.path.join(esass_dir, 'tmp2_'+str_field+'.fits')
	path_2_event_file = os.path.join(esass_dir, 'evt_'+str_field+'.fits')

	def command_ero_cal_event( nccd = 1, in_file='', out_file='', Attitude_File = Attitude_File ):
		"""
		ero_calevents Projection=AIT Attitude=${ATTITUDE} clobber=yes EvtFile=${SIMU_DIR}/erass_ccd${Nccd}_evt.fits eroEvtFile=${SIMU_DIR}/cal_erass_ccd${Nccd}_evt.fits CCDNR=${Nccd} RA=${ra_cen} Dec=${dec_cen}
		"""
		part1 = "ero_calevents Projection=AIT Attitud="+Attitude_File
		part2 = " EvtFile="+in_file
		part3 = " eroEvtFile="+out_file
		part4 = " CCDNR="+str(nccd)+" RA="+RA+" DEC="+DEC+" clobber=yes"
		full_command = part1 + part2 + part3 + part4
		return full_command

	def command_fparkey_15( p_2_file ):
		"""
		fparkey 15 ${SIMU_DIR}/cal_erass_ccd${Nccd}_evt.fits[1] PAT_SEL add=yes
		"""
		return "fparkey 15 "+ p_2_file +"[1] PAT_SEL add=yes"


	def command_ftcopy( in_file='', out_file='' ):
		"""
		ftcopy ${SIMU_DIR}/cal_erass_ccd${Nccd}_evt.fits[1]"[col *,FRAMETIME=TIME,RECORDTIME=TIME]" ${SIMU_DIR}/evt_erass_ccd${Nccd}.fits clobber=yes

		"""
		part1 = """ftcopy '""" + in_file
		part2 = """[col *,RECORDTIME=TIME,FRAMETIME=TIME]'"""
		part3 = " " + out_file
		part4 = " clobber=yes"
		full_command = part1 + part2 + part3 + part4
		return full_command


	#if os.path.isfile(path_2_tmp_file)==False:
	for NCCD in n.arange(7)+1:
		#NCCD = 1
		print("#==========================================================")
		print("# NCCD=",NCCD)
		# a: AGN
		agn_evt_files = n.array( glob.glob( os.path.join( agn_dir, 't0erass_ccd' + str(NCCD) + '_evt.fits' ) ) )
		agn_evt_OUT = n.array( [ os.path.join(esass_dir, 'a_' + os.path.basename(el) ) for el in agn_evt_files ] )
		for a_in, a_out in zip(agn_evt_files,agn_evt_OUT):
			#a_in = agn_evt_files[0]
			#a_out = agn_evt_OUT[0]
			cmd_a = command_ero_cal_event(NCCD, a_in, a_out)
			print(cmd_a)
			os.system(cmd_a)
			cmd_a_1 = command_fparkey_15( a_out)
			print(cmd_a_1)
			os.system(cmd_a_1)
			a_out_ftcopy = a_out[:-5] + '_FTCOPY.fits'
			cmd_a_2 = command_ftcopy( a_out,  a_out_ftcopy)
			print(cmd_a_2)
			os.system(cmd_a_2)

		# c: Clusters
		CL_evt_files = n.array( glob.glob( os.path.join( cluster_dir, 't0erass_ccd' + str(NCCD) + '_evt.fits' ) ) )
		CL_evt_OUT = n.array( [ os.path.join(esass_dir, 'c_' + os.path.basename(el) ) for el in CL_evt_files ] )
		for c_in, c_out in zip(CL_evt_files, CL_evt_OUT):
			#c_in = CL_evt_files[0]
			#c_out = CL_evt_OUT[0]
			cmd_c = command_ero_cal_event(NCCD, c_in, c_out)
			print(cmd_c)
			os.system(cmd_c)
			cmd_c_1 = command_fparkey_15( c_out)
			print(cmd_c_1)
			os.system(cmd_c_1)
			c_out_ftcopy = c_out[:-5] + '_FTCOPY.fits'
			cmd_c_2 = command_ftcopy( c_out,  c_out_ftcopy)
			print(cmd_c_2)
			os.system(cmd_c_2)
		"""
		# d: Stars
		ST_evt_files = n.array( glob.glob( os.path.join( stars_dir, 'simulated_photons_ccd' + str(NCCD) + '.fits' ) ) )
		ST_evt_OUT = n.array( [ os.path.join(esass_dir, 'd_' + os.path.basename(el) ) for el in ST_evt_files ] )
		for d_in, d_out in zip(ST_evt_files, ST_evt_OUT):
			#d_in = ST_evt_files[0]
			#d_out = ST_evt_OUT[0]
			cmd_d = command_ero_cal_event(NCCD, d_in, d_out)
			print(cmd_d)
			os.system(cmd_d)
			cmd_d_1 = command_fparkey_15( d_out)
			print(cmd_d_1)
			os.system(cmd_d_1)
			d_out_ftcopy = d_out[:-5] + '_FTCOPY.fits'
			cmd_d_2 = command_ftcopy( d_out,  d_out_ftcopy)
			print(cmd_d_2)
			os.system(cmd_d_2)
		"""
	"""
	# b: Background
	BG_evt_files = n.array( glob.glob( os.path.join( bg_dir, 'evt_particle_???.fits' ) ) )
	BG_evt_OUT = n.array( [ os.path.join(esass_dir, 'b_' + os.path.basename(el) ) for el in BG_evt_files ] )
	for b_in, b_ou in zip(BG_evt_files, BG_evt_OUT):
		b_out_ftcopy = b_ou[:-5] + '_FTCOPY.fits'
		cmd_b_2 = command_ftcopy( b_in,  b_out_ftcopy)
		print(cmd_b_2)
		os.system(cmd_b_2)

	# concatenates events
	path_2_list = os.path.join(esass_dir, 'FTCOPY.list')
	cmd_list = "ls " + os.path.join(esass_dir, "*_FTCOPY.fits") + " > " + path_2_list
	print( cmd_list )
	os.system( cmd_list )
	"""
