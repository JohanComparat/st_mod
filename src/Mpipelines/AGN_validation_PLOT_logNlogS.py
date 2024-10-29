"""
Plots the logNlogS
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

from scipy.interpolate import interp1d

from astropy.table import Table

import numpy as np

print('plots logNlogS')
print('------------------------------------------------')
print('------------------------------------------------')
validation_dir           = os.path.join(os.environ['GIT_STMOD'], 'data', 'validation','validation_AGN')
validation_dir_lNlS = os.path.join(validation_dir, 'XrayLogNlogS')
os.system('mkdir -p ' + validation_dir_lNlS            )

LC_dirs = np.array([ 'FullSky', 'LC1800', 'LC0060', 'LC0002' ])[::-1]
area = {}
area['FullSky'] = 129600. / np.pi
area['LC1800'] = 1800.
area['LC0060'] = 60.
area['LC0002'] = 3.42

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
all_suffixes = np.array([ 	'sigma_0.8_fsat_0.0'  ,
							'sigma_0.8_fsat_8.0'  ,
							'sigma_0.8_fsat_20.0'  ,
							'sigma_0.6_fsat_0.0'  ,
							'sigma_0.6_fsat_8.0'  ,
							'sigma_0.4_fsat_0.0'  ,
							'sigma_0.4_fsat_8.0'  ,
							'sigma_1.0_fsat_10.0' ])

def get_lnls(LC_dir = 'LC0002', suffix = 'sigma_0.8_fsat_0.0', z_dir=''):
	print(LC_dir, suffix, z_dir)
	DATA = Table()
	if z_dir=='':
		logNlogS_files = np.array(glob.glob(os.path.join(os.environ['UCHUU'], LC_dir, 'z?p??', 'logNlogS_AGN_list_'+suffix+'.fits')))
	else:
		logNlogS_files = np.array(glob.glob(os.path.join(os.environ['UCHUU'], LC_dir, z_dir, 'logNlogS_AGN_list_'+suffix+'.fits')))
	logNlogS_soft = []
	logNlogS_hard = []
	print(logNlogS_files)
	for p_2_logNlogS in logNlogS_files:
		t_out = Table.read(p_2_logNlogS)
		logNlogS_hard.append(t_out['dNdlog10S_hard'])
		logNlogS_soft.append(t_out['dNdlog10S_soft'])

	logNlogS_hard = np.transpose(logNlogS_hard)
	logNlogS_soft = np.transpose(logNlogS_soft)
	#print(logNlogS_hard.shape)
	#print(logNlogS_soft.shape)
	#print(np.sum(logNlogS_hard, axis=1))
	#print(np.sum(logNlogS_soft, axis=1))
	DATA['FX_lo'] =	t_out['FX_lo']
	DATA['FX_hi'] =	t_out['FX_hi']
	DATA['dNdlog10S_hard'] =	np.sum(logNlogS_hard, axis=1)
	DATA['dNdlog10S_soft'] =	np.sum(logNlogS_soft, axis=1)
	#print(np.cumsum(DATA['dNdlog10S_hard'][::-1])[::-1])
	#print(np.cumsum(DATA['dNdlog10S_soft'][::-1])[::-1])
	DATA['dNdlog10S_CM_soft'] = np.cumsum(DATA['dNdlog10S_soft'][::-1])[::-1]
	DATA['dNdlog10S_CM_hard'] = np.cumsum(DATA['dNdlog10S_hard'][::-1])[::-1]
	return DATA

for LC_dir in LC_dirs:
	DATA = {}
	for suffix in all_suffixes:
		DATA[suffix]  = get_lnls(LC_dir = LC_dir, suffix = 'sigma_0.8_fsat_0.0' )

	plt.figure(1, (6, 6))
	for suffix in all_suffixes:
		plt.plot((DATA[suffix]['FX_lo']+DATA[suffix]['FX_hi'])/2., DATA[suffix]['dNdlog10S_CM_soft']/area[LC_dir], lw=2, ls='dashed', label=suffix)

	# Georgakakis 2008
	path_2_logNlogS_data = os.path.join(
		os.environ["GIT_STMOD"],
		'data/validation/validation_AGN/literature_data',
		'logNlogS_Georgakakis_08_AGN.data')
	x_data, y_data = np.loadtxt(path_2_logNlogS_data, unpack=True)
	plt.plot(np.log10(x_data),y_data, lw=3, ls='dotted', color='g', label='G08')

	# Merloni 2012
	path_2_logNlogS_data = os.path.join(
		os.environ["GIT_STMOD"],
		'data/validation/validation_AGN/literature_data',
		'logNlogS_Merloni_12_AGN.data')
	x_data, y_data = np.loadtxt(path_2_logNlogS_data, unpack=True)
	plt.plot(np.log10(x_data),y_data, lw=3, ls='dotted', color='r', label='M12')

	# Mateos 2008
	path_2_logNlogS_data = os.path.join(
		os.environ["GIT_STMOD"],
		'data/validation/validation_AGN/literature_data',
		'logNlogS_Mateos_08_AGN.data')
	x_data, y_data, err = np.loadtxt(path_2_logNlogS_data, unpack=True)
	plt.plot(x_data,y_data, lw=3, ls='dotted', color='b', label='M08')

	plt.xlabel('log10(F_X[0.5-2 keV])')
	plt.ylabel('N(>F_X) [/deg2]')
	plt.legend(frameon=False, loc=0)
	plt.yscale('log')
	plt.xlim((-19, -11.5))
	plt.ylim((1e-2, 1e5))
	plt.grid()
	plt.tight_layout()
	plt.savefig(os.path.join(validation_dir_lNlS, LC_dir+"_logN_logS_soft_AGN.png"))
	plt.clf()
	print(os.path.join(validation_dir_lNlS, LC_dir+"_logN_logS_soft_AGN.png"), 'written')


	plt.figure(1, (6, 6))

	ref_line = interp1d( (DATA[suffix]['FX_lo']+DATA[suffix]['FX_hi'])/2., DATA[suffix]['dNdlog10S_CM_soft']/area[LC_dir] )

	# Georgakakis 2008
	path_2_logNlogS_data = os.path.join(
		os.environ["GIT_STMOD"],
		'data/validation/validation_AGN/literature_data',
		'logNlogS_Georgakakis_08_AGN.data')
	x_data, y_data = np.loadtxt(path_2_logNlogS_data, unpack=True)
	plt.plot(np.log10(x_data),y_data/ref_line(np.log10(x_data)), lw=3, ls='dotted', color='g', label='G08')

	# Merloni 2012
	path_2_logNlogS_data = os.path.join(
		os.environ["GIT_STMOD"],
		'data/validation/validation_AGN/literature_data',
		'logNlogS_Merloni_12_AGN.data')
	x_data, y_data = np.loadtxt(path_2_logNlogS_data, unpack=True)
	plt.plot(np.log10(x_data),y_data/ref_line(np.log10(x_data)), lw=3, ls='dotted', color='r', label='M12')

	# Mateos 2008
	path_2_logNlogS_data = os.path.join(
		os.environ["GIT_STMOD"],
		'data/validation/validation_AGN/literature_data/logNlogS_Mateos_08_AGN.data')
	x_data, y_data, err = np.loadtxt(path_2_logNlogS_data, unpack=True)
	plt.plot(x_data,y_data/ref_line(x_data), lw=3, ls='dotted', color='b', label='M08')

	plt.xlabel('log10(F_X[0.5-2 keV])')
	plt.ylabel('N(>F_X)/(N Uchuu AGN,>F_X) '+suffix)
	plt.legend(frameon=False, loc=0)
	plt.xlim((-19, -11.5))
	plt.ylim((0.3, 1.7))
	plt.tight_layout()
	plt.grid()
	plt.savefig(os.path.join(validation_dir_lNlS, LC_dir+"_logN_logS_soft_AGN_ratio.png"))
	plt.clf()


	plt.figure(1, (6, 6))
	for suffix in all_suffixes:
		plt.plot((DATA[suffix]['FX_lo']+DATA[suffix]['FX_hi'])/2., DATA[suffix]['dNdlog10S_CM_hard']/area[LC_dir], lw=2, ls='dashed', label=suffix)

	# Merloni 2012
	path_2_logNlogS_data = os.path.join(
		os.environ["GIT_STMOD"],
		'data/validation/validation_AGN/literature_data',
		'logNlogS_Merloni_12_AGN_hard.data')
	x_data, y_data = np.loadtxt(path_2_logNlogS_data, unpack=True)
	plt.plot(np.log10(x_data),y_data, lw=3, ls='dotted', color='r', label='M12')

	plt.xlabel('log10(F_X[2-10 keV])')
	plt.ylabel('log10(>F_X) [/deg2]')
	plt.legend(frameon=False, loc=0)
	plt.yscale('log')
	plt.xlim((-19, -11.5))
	#plt.ylim((-2, 5))
	#lt p.title('Mocks')
	plt.grid()
	plt.savefig(os.path.join(validation_dir_lNlS, LC_dir+"_logN_logS_hard_AGN.png"))
	plt.clf()


for jj, z_dir in enumerate(all_z_dirs):#[:22]:
	print(z_dir)
	DATA = {}
	suffix = all_suffixes[0]
	if jj<22:
		all_LC_dirs = LC_dirs[::-1]
	elif jj>=22 and jj<26:
		all_LC_dirs = LC_dirs[::-1][1:]
	elif jj>=26 and jj<35:
		all_LC_dirs = LC_dirs[::-1][2:]
	elif jj>35 :
		all_LC_dirs = LC_dirs[::-1][3:]

	for LC_dir in all_LC_dirs:
		DATA[LC_dir]  = get_lnls(LC_dir = LC_dir, suffix = suffix, z_dir=z_dir )

	plt.figure(1, (6, 6))
	for LC_dir in all_LC_dirs:
		plt.plot((DATA[LC_dir]['FX_lo']+DATA[LC_dir]['FX_hi'])/2., DATA[LC_dir]['dNdlog10S_CM_soft']/area[LC_dir], lw=2, ls='dashed', label=LC_dir)

	# Georgakakis 2008
	path_2_logNlogS_data = os.path.join(
		os.environ["GIT_STMOD"],
		'data/validation/validation_AGN/literature_data',
		'logNlogS_Georgakakis_08_AGN.data')
	x_data, y_data = np.loadtxt(path_2_logNlogS_data, unpack=True)
	plt.plot(np.log10(x_data), y_data, lw=3, ls='dotted', color='g', label='G08')

	# Merloni 2012
	path_2_logNlogS_data = os.path.join(
		os.environ["GIT_STMOD"],
		'data/validation/validation_AGN/literature_data',
		'logNlogS_Merloni_12_AGN.data')
	x_data, y_data = np.loadtxt(path_2_logNlogS_data, unpack=True)
	plt.plot(np.log10(x_data),y_data, lw=3, ls='dotted', color='r', label='M12')

	# Mateos 2008
	path_2_logNlogS_data = os.path.join(
		os.environ["GIT_STMOD"],
		'data/validation/validation_AGN/literature_data',
		'logNlogS_Mateos_08_AGN.data')
	x_data, y_data, err = np.loadtxt(path_2_logNlogS_data, unpack=True)
	plt.plot(x_data,y_data, lw=3, ls='dotted', color='b', label='M08')

	plt.xlabel('log10(F_X[0.5-2 keV])')
	plt.ylabel('log10(>F_X) [/deg2]')
	plt.legend(frameon=False, loc=0)
	plt.yscale('log')
	plt.xlim((-19, -11.5))
	#plt.ylim((-2, 5))
	#lt p.title('Mocks')
	plt.grid()
	plt.savefig(os.path.join(validation_dir_lNlS, z_dir+"_logN_logS_soft_AGN.png"))
	plt.clf()
	print(os.path.join(validation_dir_lNlS, z_dir+"_logN_logS_soft_AGN.png"), 'written')

