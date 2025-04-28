"""
Simulate eROSITA event data from sky model using SIXTE.
A. Malyali, 2019. amalyali@mpe.mpg.de
modification by J.C:
 - more flexible with directories and multiple simput files
 - merge of the 3 gti files obtained using 'mgtime'

 ero_vis GTIfile=/data26s/comparat/simulations/erosim/eRASS8_events_cluster_UNIT_fA1i_DIR_June2020_cat_CLU_SIMPUT_b8_CM_0_pixS_20.0/000/erass_c_000000_N_2.gti Simput=/data24s/comparat/simulation/UNIT/ROCKSTAR_HALOS/fixedAmp_InvPhase_001/cat_CLU_SIMPUT_b8_CM_0_pixS_20.0/c_000000_N_2.fit Exposure=55200.000000 Attitude=/data40s/erosim/eRASS/eRASS_4yr_epc85_att.fits TSTART=629332205.000000 dt=0.5 visibility_range=1.0 clobber=yes

"""
import subprocess
import os
import errno
import sys, glob
import healpy
import numpy as n
from astropy.io import fits
pix_ids = n.arange(healpy.nside2npix(8) )
ra_cen_s, dec_cen_s = healpy.pix2ang(8, pix_ids, nest=True, lonlat=True)

class Simulator:
	"""
	SIXTE simulator for eROSITA observations.
	1. Compute GTI file for given simput
	2. Simulate eROSITA observations of simput, using GTI to speed things up.
	"""

	def __init__(self, block_id, with_bkg_par, t_start, exposure, seed, simput, data_dir, ra_cen, dec_cen, p_2_att_file):
		# def __init__(self, with_bkg_par, t_start, exposure, seed, simput,
		# simput2, simput3):
		"""
		:param with_bkg_par: Simulate with particle background.
		:param t_start: Start time of simulation. Input units of [s]
		:param exposure: Length of time to simulate for after t_start
		:param seed: Seed for random number generator.
		:param simput: Simput files (ie. the sky model)
		"""
		self._block_id = 't'+str(block_id)
		self._with_bkg_par = bool(with_bkg_par)
		self._t_start = float(t_start)  # secs
		self._exposure = float(exposure)
		self._seed = int(seed)
		self._simput = simput
		self._N_simputs = len(simput)
		self._data_dir = data_dir
		self._ra_cen = ra_cen
		self._dec_cen = dec_cen
		self._p_2_att_file = p_2_att_file

	def make_event_directory(self):
		"""
		Check for whether directory exists and create if not.
		"""
		try:
			os.makedirs(self._data_dir)
		except OSError as e:
			print('directory already exists')

	def compute_gti(self, simput_ero, gti_file ):
		"""
		Compute the GTI (good time interval) file given the simput file.
		Use this as input to SIXTE call for reducing computational expense.
		"""
		cmd = ["ero_vis",
				"GTIfile=%s" % gti_file,
				"Simput=%s" % simput_ero,
				"Exposure=%f" % self._exposure,
				"Attitude=%s" % self._p_2_att_file,
				"TSTART=%f" % self._t_start,
				#"RA=%s" % self._ra_cen,
				#"Dec=%s" % self._dec_cen,
				"dt=0.5",
				"visibility_range=1.0",
				"clobber=yes"
				]

		command = " ".join(cmd)
		print('=========================================')
		print(command)
		print('=========================================')
		os.system(command)
		hd=fits.open(gti_file)
		number_of_lines = hd[1].header['NAXIS2']
		print('file N lines', number_of_lines)
		if number_of_lines==0:
			cmd_del = 'rm '+gti_file
			print('=========================================')
			print(cmd_del)
			print('=========================================')
			os.system(cmd_del)

	def run_sixte(self):
		"""
		Launch erosim from python.
		"""
		print('N simput files', self._N_simputs)
		prefix = os.path.join(self._data_dir, self._block_id+"erass_" )
		if self._N_simputs==1:
			cmd = ["erosim",
					"Simput=%s" % self._simput[0],
					"Prefix=%s" % prefix,
					"Attitude=%s" % self._p_2_att_file,
					"RA=%s" % self._ra_cen,
					"Dec=%s" % self._dec_cen,
					"GTIFile=%s/erass.gti" % self._data_dir,
					"TSTART=%s" % self._t_start,
					"Exposure=%s" % self._exposure,
					"MJDREF=51543.875",
					"dt=0.5",
					"Seed=%s" % self._seed,
					"clobber=yes",
					"chatter=3"
					]

		if self._N_simputs>=2:
			cmd = ["erosim", "Simput=%s" % self._simput[0] ]
			for jj in n.arange(len(self._simput))[1:]:
				cmd.append("Simput"+str(int(jj+1))+"=%s" % self._simput[jj])
			cmd_end = ["Prefix=%s" % prefix,
					"Attitude=%s" % self._p_2_att_file,
					"RA=%s" % self._ra_cen,
					"Dec=%s" % self._dec_cen,
					"GTIFile=%s/erass.gti" % self._data_dir,
					"TSTART=%s" % self._t_start,
					"Exposure=%s" % self._exposure,
					"MJDREF=51543.875",
					"dt=0.5",
					"Seed=%s" % self._seed,
					"clobber=yes",
					"chatter=3"
					]
			for el in cmd_end:
				cmd.append(el)

		if self._with_bkg_par is True:
			cmd.append("Background=yes")
		else:
			cmd.append("Background=no")
		command = " ".join(cmd)
		print('=========================================')
		print(command)
		print('=========================================')
		os.system(command)

	def run_all(self):
		"""
		Run SIXTE simulation of eRASS 8
		"""
		print('make event directory')
		self.make_event_directory()
		print('compute gti with ero_vis')
		if self._N_simputs==1:
			path_to_gti = os.path.join(self._data_dir, "erass.gti")
			self.compute_gti(self._simput[0], path_to_gti)

		if self._N_simputs>=2:
			for el in self._simput:
				gti_file = os.path.join(self._data_dir, "erass_" + os.path.basename(el)[:-4] + ".gti")
				self.compute_gti(el, gti_file)
			print('merges gti files')
			path_to_gti_list = os.path.join(self._data_dir, "gti.list")
			#print(path_to_gti_list)
			path_to_gti = os.path.join(self._data_dir, "erass.gti")
			#print(path_to_gti)
			gti_files = glob.glob(self._data_dir + "/*.gti")
			if len(gti_files)>0:
				command_list = "ls " + self._data_dir + "/*.gti > " + path_to_gti_list
				print('=========================================')
				print(command_list)
				print('=========================================')
				os.system(command_list)
				# ls /data40s/erosim/eRASS/sixte/000/erass_SIMPUT_000000*.gti >
				# gti.list
				command_merge = "mgtime ingtis=@" + path_to_gti_list + " outgti=" + path_to_gti + " merge=OR"
				print('=========================================')
				print(command_merge)
				print('=========================================')
				os.system(command_merge)
		print('is there a gti file ?', os.path.isfile(path_to_gti))
		if os.path.isfile(path_to_gti):
			print('run sixte with erosim')
			self.run_sixte()
			list_2_delete = n.array(glob.glob(os.path.join(self._data_dir,  "*.gti" )))
			print('============================================')
			print('============================================')
			print('deleting', list_2_delete)
			print('============================================')
			print('============================================')
			del_list = n.array([ os.remove(el) for el in  list_2_delete])

if __name__ == '__main__':
	#erass_option = "eRASS8" # sys.argv[4]
	#print(erass_option)
	seed = 42
	tile_id = sys.argv[1]  # '355'
	erass_option = sys.argv[2] # "eRASS1" # sys.argv[4]
	ra_cen = ra_cen_s[int(tile_id)]
	dec_cen = dec_cen_s[int(tile_id)]
	data_dir = os.path.join(os.environ['HOME'], "workspace/erosim/eRASS8_events_BG_FG/" + tile_id )
	print(data_dir)
	topdir = os.path.join(os.environ['HOME'], 'workspace/erosim/simput/bkg_erosita_simput_full_sky/')
	simput_files = n.array([ os.path.join(topdir, 'SIMPUT_000' + tile_id + '.fits') ])
	print(topdir  )
	print( tile_id, simput_files)

	if erass_option == "eRASS8":
		p_2_att_file = "/home/idies/workspace/erosim/erosita_attitude/eRASS_4yr_epc85_att.fits"
		t_starts = n.array([ 617943605 ])
		t_stops = n.array([ 744174005 ])

	for jj, (t0, t1) in enumerate( zip( t_starts, t_stops ) ):
		print('+++++++++++++++++++++++++++++++++++++++++++++++++')
		print('+++++++++++++++++++++++++++++++++++++++++++++++++')
		print('STEP', jj, (t0, t1))
		print('+++++++++++++++++++++++++++++++++++++++++++++++++')
		print('+++++++++++++++++++++++++++++++++++++++++++++++++')
		bkg = 0
		t_start = t0 # 617943605.0
		exposure = t1 - t0 # 31536000 * 4 # = 4 years  # 31536000 = 1year # 15750000 = 1/2 year
		# Launch...
		# 3 files
		#Simulator(bkg, t_start, exposure, seed, simput_file_1, simput_file, simput_file_2).run_all()
		# 2 files
		Simulator(
			jj,
			bkg,
			t_start,
			exposure,
			seed,
			simput_files,
			data_dir,
			ra_cen,
			dec_cen,
			p_2_att_file).run_all()


