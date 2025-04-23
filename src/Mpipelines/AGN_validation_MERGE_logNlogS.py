"""
What it does
------------

Computes and tabulates the logNlogS, logNlogR, 

Creates figures with the tabulate data

References
----------

Command to run
--------------

python3 003_2_agn_compute_XLF_logNlogS_R.py environmentVAR ftyp

arguments
---------

environmentVAR: environment variable linking to the directory where files are e.g. "MD10"
It will then work in the directory : $environmentVAR/hlists/fits/

ftyp: 'all' or 'sat', type of AGN populating a central halo or satellite haloes.

Dependencies
------------

import time, os, sys, numpy, scipy, astropy, h5py, matplotlib

"""
import time
t0 = time.time()
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as plt
import sys, os, glob
from scipy.integrate import quad

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

#import extinction

from scipy.interpolate import interp1d

from astropy.table import Table

import numpy as np

print('MERGE logNlogS files')
print('------------------------------------------------')
print('------------------------------------------------')

LC_dirs = np.array([ 'FullSky'])#, 'LC1800', 'LC0060', 'LC0002' ])#[::-1]

all_z_dirs = np.array([ 'z0p02',
						'z0p05',
						'z0p09',
						'z0p14',
						'z0p19',
						'z0p25',
						'z0p30',
						'z0p36',
						'z0p43',
						'z0p49',
						'z0p56',
						'z0p63',
						'z0p70',
						'z0p78',
						'z0p86',
						'z0p94',
						'z1p03',
						'z1p12',
						'z1p22',
						'z1p32',
						'z1p43',
						'z1p54',
						'z1p65',
						'z1p77',
						'z1p90',
						'z2p03',
						'z2p17',
						'z2p31',
						'z2p46',
						'z2p62',
						'z2p78',
						'z2p95',
						'z3p13',
						'z3p32',
						'z3p61',
						'z3p93',
						'z4p27',
						'z4p63',
						'z5p15',
						'z5p73' ])
all_suffixes = np.array([ 	#'sigma_0.8_fsat_0.0'  ,
							'sigma_0.8_fsat_8.0'  ,
							#'sigma_0.8_fsat_20.0'  ,
							#'sigma_0.6_fsat_0.0'  ,
							#'sigma_0.6_fsat_8.0'  ,
							#'sigma_0.4_fsat_0.0'  ,
							#'sigma_0.4_fsat_8.0'  ,
							#'sigma_1.0_fsat_10.0' ])
							])

def get_lnls(LC_dir = 'LC0002', suffix = 'sigma_0.8_fsat_0.0', z_dir=''):
	DATA = Table()
	if z_dir=='':
		logNlogS_files = np.array(glob.glob(os.path.join(os.environ['UCHUU'], LC_dir, 'z?p??', 'replication_*', 'logNlogS_AGN_list_'+suffix+'.fits')))
	else:
		logNlogS_files = np.array(glob.glob(os.path.join(os.environ['UCHUU'], LC_dir, z_dir, 'replication_*', 'logNlogS_AGN_list_'+suffix+'.fits')))
	logNlogS_soft = []
	logNlogS_hard = []
	#print(logNlogS_files)
	if len(logNlogS_files)>=1:
		for p_2_logNlogS in logNlogS_files:
			t_out = Table.read(p_2_logNlogS)
			#print(t_out)
			logNlogS_hard.append(t_out['N_hard'] )#/ t_out['area'])
			logNlogS_soft.append(t_out['N_soft'] )#/ t_out['area'])

		logNlogS_hard = np.transpose(logNlogS_hard)
		logNlogS_soft = np.transpose(logNlogS_soft)
		DATA['FX_lo'] =	t_out['FX_lo']
		DATA['FX_hi'] =	t_out['FX_hi']
		DATA['dNdlog10S_hard'] =	np.sum(logNlogS_hard, axis=1)
		DATA['dNdlog10S_soft'] =	np.sum(logNlogS_soft, axis=1)
		DATA['dNdlog10S_CM_soft'] = np.cumsum(DATA['dNdlog10S_soft'][::-1])[::-1]
		DATA['dNdlog10S_CM_hard'] = np.cumsum(DATA['dNdlog10S_hard'][::-1])[::-1]
		return DATA

for suffix in all_suffixes: #= 'sigma_0.8_fsat_0.0'
	print('='*100)
	print(suffix)
	for LC_dir in LC_dirs:
		print('*'*100)
		print(LC_dir)
		for z_dir in all_z_dirs :
			print(z_dir)
			DATA  = get_lnls(LC_dir = LC_dir, suffix=suffix, z_dir = z_dir)
			try:
				DATA.write(os.path.join(os.environ['UCHUU'], LC_dir, z_dir, 'logNlogS_AGN_list_'+suffix+'.fits'), overwrite = True)
				print(os.path.join(os.environ['UCHUU'], LC_dir, z_dir, 'logNlogS_AGN_list_'+suffix+'.fits'), 'written')
			except(AttributeError):
				print('no data')
				continue
