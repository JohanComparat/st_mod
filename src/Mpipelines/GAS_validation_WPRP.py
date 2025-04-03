import time
t0 = time.time()
import numpy as np
import numpy as n
import healpy as hp
from astropy.table import Table, Column, vstack, hstack
import sys, os, glob

import astropy.units as u
import astropy.constants as cc
import astropy.io.fits as fits
from scipy.stats import norm
speed_light = cc.c.to(u.km/u.s).value
from Corrfunc.utils import convert_3d_counts_to_cf
from Corrfunc.mocks.DDtheta_mocks import DDtheta_mocks
import Corrfunc
from Corrfunc.mocks.DDrppi_mocks import DDrppi_mocks
from Corrfunc.io import read_catalog
from Corrfunc.utils import convert_rp_pi_counts_to_wp
from Corrfunc.theory.DD import DD

from scipy.interpolate import interp1d
from astropy.cosmology import FlatLambdaCDM
cosmoUCHUU = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
h = 0.6774
cosmo = cosmoUCHUU
zs = np.arange(0.0000001, 7.1, 0.001)
dm_itp = interp1d(zs, cosmo.distmod(zs).value)

# cosmology setup
z_array = n.arange(0.001, 7.5, 0.001)
d_C = cosmo.comoving_distance(z_array)
dc_mpc = (d_C).value
dc_interpolation = interp1d(z_array, dc_mpc)
z_interpolation = interp1d(dc_mpc, z_array)

def get_xyz(RARA, DECDEC, ZZZZ):
	phi   = ( RARA   - 180 ) * n.pi / 180.
	theta = ( DECDEC + 90 ) * n.pi / 180.
	rr    = dc_interpolation( ZZZZ )
	xx = rr * n.cos( phi   ) * n.sin( theta )
	yy = rr * n.sin( phi   ) * n.sin( theta )
	zz = rr * n.cos( theta )
	return n.array(list(xx)), n.array(list(yy)), n.array(list(zz))

def get_radec(xn, yn, zn):
	rr = (xn**2 + yn**2 + zn**2)**0.5
	# angular coordinates
	theta = n.arccos(zn / rr) * 180 / n.pi
	phi = n.arctan2(yn, xn) * 180 / n.pi
	ra = phi + 180.
	dec = theta - 90.
	return ra, dec, z_interpolation(rr)

z_dir = sys.argv[1]
#LC_dir = sys.argv[2] # 'FullSky'
LC_dir='FullSky'

validation_dir       = os.path.join(os.environ['GIT_STMOD_DATA'], 'data', 'validation','validation_GAS')
validation_dir_WPRP = os.path.join(validation_dir, 'WPRP')
os.system('mkdir -p ' + os.path.join(validation_dir_WPRP, z_dir) )

sys.path.append( os.path.join(os.environ['GIT_STMOD'], 'src') )
from models import GAL as GG
from models import GAS as MG

C_GAL = GG.GAL(z_dir, LC_dir=LC_dir)
C_GAS = MG.GAS(z_dir, b_HS=0.8, logM500c_min=11., logFX_min=-18, LC_dir=LC_dir)

def tabulate_wprp_clustering_noW(RA, DEC, Z, rand_RA , rand_DEC, rand_Z, out_file='test.fits', pimax = 100.0 ):
	"""
	wprp direct estimate
	path_2_data : path to the catalogue to correlate
	path_2_random
	"""
	#
	CZ = Z * speed_light
	rand_CZ = rand_Z * speed_light
	N = len(RA)
	rand_N = len(rand_RA)
	#print(N, rand_N, out_file, time.time()-t0)
	bins = 10**np.arange(-1.6, 1.81, 0.2)
	nbins = len(bins)-1
	#print('bins', bins, bins.shape)
	x = (bins[1:]+bins[:-1])/2.
	cosmology = 2
	nthreads = 16
	# Auto pairs counts in DD
	autocorr=1
	DD_counts = DDrppi_mocks(autocorr, cosmology, nthreads, pimax, bins,
								np.array(list(RA)).astype('float'),
								np.array(list(DEC)).astype('float'),
								np.array(list(CZ)).astype('float') )#, is_comoving_dist=True)
	#print('DD',DD_counts['npairs'], DD_counts['npairs'].shape)
	# Auto pairs counts in DR
	autocorr=0
	DR_counts = DDrppi_mocks(autocorr, cosmology, nthreads, pimax, bins,
							np.array(list(RA)).astype('float'),
							np.array(list(DEC)).astype('float'),
							np.array(list(CZ)).astype('float'),
							RA2=rand_RA.astype('float'),
							DEC2=rand_DEC.astype('float'),
							CZ2=rand_CZ.astype('float'))
	#print('DR',DR_counts['npairs'])
	# Auto pairs counts in RR
	autocorr=1
	RR_counts = DDrppi_mocks(autocorr, cosmology, nthreads, pimax, bins,
							rand_RA.astype('float'),
							rand_DEC.astype('float'),
							rand_CZ.astype('float'))
	#print('RR',RR_counts['npairs'])
	# All the pair counts are done, get the angular correlation function
	wp = convert_rp_pi_counts_to_wp(N, N, rand_N, rand_N, DD_counts, DR_counts, DR_counts, RR_counts, nbins, pimax)
	t = Table()
	t.add_column(Column(data = bins[:-1], name='rp_min', unit='Mpc'  ) )
	t.add_column(Column(data = bins[1:], name='rp_max', unit='Mpc'  ) )
	t.add_column(Column(data = x, name='rp_mid', unit='Mpc'  ) )
	t.add_column(Column(data = wp, name='wprp', unit=''  ) )
	t.add_column(Column(data = np.ones_like(x) * N, name='N_data', unit=''  ) )
	t.add_column(Column(data = np.ones_like(x) * rand_N, name='N_random', unit=''  ) )
	t.add_column(Column(data = np.ones_like(x) * pimax, name='pimax', unit=''  ) )
	# repeat when removing random 10%
	print(out_file)
	t.write(out_file, overwrite=True, format='fits')
	print(out_file, time.time()-t0, 's')

def tabulate_wprp_clustering_Cross(RA, DEC, Z, rand_RA , rand_DEC, rand_Z, RA2, DEC2, Z2, rand_RA2 , rand_DEC2, rand_Z2, out_file='test.fits', pimax = 100.0 ):
	CZ = Z * speed_light
	rand_CZ = rand_Z * speed_light
	CZ2 = Z2 * speed_light
	rand_CZ2 = rand_Z2 * speed_light
	N = len(RA)
	rand_N = len(rand_RA)
	N2 = len(RA2)
	rand_N2 = len(rand_RA2)
	#print(N, rand_N, out_file, time.time()-t0)
	bins = 10**np.arange(-2.0, 1.81, 0.05)
	nbins = len(bins)-1
	#print('bins', bins, bins.shape)
	x = (bins[1:]+bins[:-1])/2.
	cosmology = 2
	nthreads = 16
	# Auto pairs counts in DD
	autocorr=0
	DD_counts = DDrppi_mocks(autocorr, cosmology, nthreads, pimax, bins,
								np.array(list(RA)).astype('float'),
								np.array(list(DEC)).astype('float'),
								np.array(list(CZ)).astype('float'),
								RA2= np.array(list(RA2)).astype('float'),
								DEC2=np.array(list(DEC2)).astype('float'),
								CZ2= np.array(list(CZ2)).astype('float')
								)#, is_comoving_dist=True)
	# Auto pairs counts in DR
	D1R_counts = DDrppi_mocks(autocorr, cosmology, nthreads, pimax, bins,
							np.array(list(RA)).astype('float'),
							np.array(list(DEC)).astype('float'),
							np.array(list(CZ)).astype('float'),
							RA2 =rand_RA2.astype('float'),
							DEC2=rand_DEC2.astype('float'),
							CZ2 =rand_CZ2.astype('float'))
	D2R_counts = DDrppi_mocks(autocorr, cosmology, nthreads, pimax, bins,
							np.array(list(RA2)).astype('float'),
							np.array(list(DEC2)).astype('float'),
							np.array(list(CZ2)).astype('float'),
							RA2 =rand_RA.astype('float'),
							DEC2=rand_DEC.astype('float'),
							CZ2 =rand_CZ.astype('float'))
	# Auto pairs counts in RR
	RR_counts = DDrppi_mocks(autocorr, cosmology, nthreads, pimax, bins,
							rand_RA.astype('float'),
							rand_DEC.astype('float'),
							rand_CZ.astype('float'),
							RA2 =rand_RA2.astype('float'),
							DEC2=rand_DEC2.astype('float'),
							CZ2 =rand_CZ2.astype('float'))
	#print('RR',RR_counts['npairs'])
	# All the pair counts are done, get the angular correlation function
	wp = convert_rp_pi_counts_to_wp(N, N2, rand_N, rand_N2, DD_counts, D1R_counts, D2R_counts, RR_counts, nbins, pimax)
	t = Table()
	t.add_column(Column(data = bins[:-1], name='rp_min', unit='Mpc'  ) )
	t.add_column(Column(data = bins[1:], name='rp_max', unit='Mpc'  ) )
	t.add_column(Column(data = x, name='rp_mid', unit='Mpc'  ) )
	t.add_column(Column(data = wp, name='wprp', unit=''  ) )
	t.add_column(Column(data = np.ones_like(x) * N, name='N_data', unit=''  ) )
	t.add_column(Column(data = np.ones_like(x) * rand_N, name='N_random', unit=''  ) )
	t.add_column(Column(data = np.ones_like(x) * pimax, name='pimax', unit=''  ) )
	print(out_file)
	t.write(out_file, overwrite=True, format='fits')
	print(out_file, time.time()-t0, 's')

LX_mins = np.arange(42.5, 44.1, 0.1)
Ms_mins = np.arange(10, 11.75, 0.25)[::-1]

for meta in C_GAL.LC_MetaData[:1]:#[(enough_area)&(small_difference_minmax_1)&(small_difference_minmax_2)]:
	#
	print(meta)
	# retrieve the resulting catalogues and meta data
	str_replic = 'replication_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz'])
	p_2_catalogue = os.path.join(os.environ['UCHUU'], LC_dir, z_dir, 'replication_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']), 'glist.fits')
	p_2_catal_MAG = os.path.join(os.environ['UCHUU'], LC_dir, z_dir, 'replication_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']), 'Kmatch_mags.fits')
	GAL = Table.read(p_2_catalogue)
	p_2_catal_GAS_b08 = os.path.join(os.environ['UCHUU'], LC_dir, z_dir, 'replication_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']), 'Xgas_bHS0.8.fits')
	XGA = Table.read(p_2_catal_GAS_b08)

	for Ms_min in Ms_mins:
		s1 = (np.log10(GAL['obs_sm'])>Ms_min)
		data = Table()
		data['RA'] = GAL['RA'][s1]
		data['DEC'] = GAL['DEC'][s1]
		data['Z'] = GAL['redshift_S'][s1]
		M_str = 'Ms_gt_'+str(np.round(Ms_min,2))
		z0 = n.min(data['Z'] )
		z1 = n.max(data['Z'] )
		N_data = len(data)
		print(N_data, 'galaxies')
		x_low, x_high = n.min(GAL['x'][s1]), n.max(GAL['x'][s1])
		y_low, y_high = n.min(GAL['y'][s1]), n.max(GAL['y'][s1])
		z_low, z_high = n.min(GAL['z'][s1]), n.max(GAL['z'][s1])
		size = int(N_data*40)
		uu = n.random.uniform(size=size)
		dec_fs = n.arccos(1 - 2 * uu) * 180 / n.pi - 90.
		ra_fs  = n.random.uniform(size=size) * 2 * n.pi * 180 / n.pi
		z_fs = n.tile(data['Z'], int(size*1./N_data)+1)[:size]
		RX, RY, RZ = get_xyz(ra_fs, dec_fs, z_fs)
		selectionR = ( RX >= x_low ) & ( RX <= x_high ) & ( RY >= y_low ) & ( RY <= y_high ) & ( RZ >= z_low ) & ( RZ <= z_high )
		rand = Table()
		rand['RA'] = ra_fs[selectionR]
		rand['DEC'] = dec_fs[selectionR]
		rand['Z'] = z_fs[selectionR]
		N_RD = len(rand['RA'])
		print(N_RD, 'randoms')
		p_2_2PCF_G = os.path.join(validation_dir_WPRP, z_dir, M_str+'_'+str_replic+'-wprp-pimax100-2pcf.fits' )
		print(p_2_2PCF_G)
		if os.path.isfile(p_2_2PCF_G)==False:
			tabulate_wprp_clustering_noW(data['RA'], data['DEC'], data['Z'], rand['RA'] , rand['DEC'], rand['Z'],  out_file=p_2_2PCF_G )

		for LX_min in LX_mins:
			s2 = ( XGA['CLUSTER_LX_soft_RF_R500c'] > LX_min )
			CC = Table()
			CC['RA'] = XGA['RA'][s2]
			CC['DEC'] = XGA['DEC'][s2]
			CC['Z'] = XGA['redshift_S'][s2]
			N_CLU = len(CC)
			print(N_CLU, 'clusters')
			LX_str = 'HaloLX_gt_'+str(np.round(LX_min,1))
			x_low, x_high = n.min(GAL['x'][s1]), n.max(GAL['x'][s1])
			y_low, y_high = n.min(GAL['y'][s1]), n.max(GAL['y'][s1])
			z_low, z_high = n.min(GAL['z'][s1]), n.max(GAL['z'][s1])
			size = int(N_CLU*40)
			uu = n.random.uniform(size=size)
			dec_fs = n.arccos(1 - 2 * uu) * 180 / n.pi - 90.
			ra_fs  = n.random.uniform(size=size) * 2 * n.pi * 180 / n.pi
			z_fs = n.tile(CC['Z'], int(size*1./N_CLU)+1)[:size]
			RX, RY, RZ = get_xyz(ra_fs, dec_fs, z_fs)
			selectionR = ( RX >= x_low ) & ( RX <= x_high ) & ( RY >= y_low ) & ( RY <= y_high ) & ( RZ >= z_low ) & ( RZ <= z_high )
			CCrand = Table()
			CCrand['RA'] = ra_fs[selectionR]
			CCrand['DEC'] = dec_fs[selectionR]
			CCrand['Z'] = z_fs[selectionR]
			N_CCRD = len(CCrand['RA'])
			print(N_CCRD, 'cluster randoms')
			p_2_2PCF_X = os.path.join(validation_dir_WPRP, z_dir, M_str+'_'+LX_str+'_'+str_replic+'-wprp-pimax100-2pcf.fits' )
			print(p_2_2PCF_X)
			if os.path.isfile(p_2_2PCF_X)==False and len(CC['RA'])>10:
				tabulate_wprp_clustering_Cross(
								data['RA'], data['DEC'], data['Z'],
								rand['RA'] , rand['DEC'], rand['Z'],
								CC['RA'], CC['DEC'], CC['Z'],
								CCrand['RA'] , CCrand['DEC'], CCrand['Z'],
								out_file=p_2_2PCF_X )

	#
	# star-forming
	#
	for Ms_min in Ms_mins:
		s1 = (np.log10(GAL['obs_sm'])>Ms_min)&(np.log10(GAL['obs_sfr']/GAL['obs_sm'])>=-11)
		data = Table()
		data['RA'] = GAL['RA'][s1]
		data['DEC'] = GAL['DEC'][s1]
		data['Z'] = GAL['redshift_S'][s1]
		M_str = 'Ms_gt_'+str(np.round(Ms_min,2))
		z0 = n.min(data['Z'] )
		z1 = n.max(data['Z'] )
		N_data = len(data)
		print(N_data, 'galaxies')
		x_low, x_high = n.min(GAL['x'][s1]), n.max(GAL['x'][s1])
		y_low, y_high = n.min(GAL['y'][s1]), n.max(GAL['y'][s1])
		z_low, z_high = n.min(GAL['z'][s1]), n.max(GAL['z'][s1])
		size = int(N_data*40)
		uu = n.random.uniform(size=size)
		dec_fs = n.arccos(1 - 2 * uu) * 180 / n.pi - 90.
		ra_fs  = n.random.uniform(size=size) * 2 * n.pi * 180 / n.pi
		z_fs = n.tile(data['Z'], int(size*1./N_data)+1)[:size]
		RX, RY, RZ = get_xyz(ra_fs, dec_fs, z_fs)
		selectionR = ( RX >= x_low ) & ( RX <= x_high ) & ( RY >= y_low ) & ( RY <= y_high ) & ( RZ >= z_low ) & ( RZ <= z_high )
		rand = Table()
		rand['RA'] = ra_fs[selectionR]
		rand['DEC'] = dec_fs[selectionR]
		rand['Z'] = z_fs[selectionR]
		N_RD = len(rand['RA'])
		print(N_RD, 'randoms')
		p_2_2PCF_G = os.path.join(validation_dir_WPRP, z_dir, M_str+'_BC_'+str_replic+'-wprp-pimax100-2pcf.fits' )
		print(p_2_2PCF_G)
		if os.path.isfile(p_2_2PCF_G)==False:
			tabulate_wprp_clustering_noW(data['RA'], data['DEC'], data['Z'], rand['RA'] , rand['DEC'], rand['Z'],  out_file=p_2_2PCF_G )

		for LX_min in LX_mins:
			s2 = ( XGA['CLUSTER_LX_soft_RF_R500c'] > LX_min )
			CC = Table()
			CC['RA'] = XGA['RA'][s2]
			CC['DEC'] = XGA['DEC'][s2]
			CC['Z'] = XGA['redshift_S'][s2]
			N_CLU = len(CC)
			print(N_CLU, 'clusters')
			LX_str = 'HaloLX_gt_'+str(np.round(LX_min,1))
			x_low, x_high = n.min(GAL['x'][s1]), n.max(GAL['x'][s1])
			y_low, y_high = n.min(GAL['y'][s1]), n.max(GAL['y'][s1])
			z_low, z_high = n.min(GAL['z'][s1]), n.max(GAL['z'][s1])
			size = int(N_CLU*40)
			uu = n.random.uniform(size=size)
			dec_fs = n.arccos(1 - 2 * uu) * 180 / n.pi - 90.
			ra_fs  = n.random.uniform(size=size) * 2 * n.pi * 180 / n.pi
			z_fs = n.tile(CC['Z'], int(size*1./N_CLU)+1)[:size]
			RX, RY, RZ = get_xyz(ra_fs, dec_fs, z_fs)
			selectionR = ( RX >= x_low ) & ( RX <= x_high ) & ( RY >= y_low ) & ( RY <= y_high ) & ( RZ >= z_low ) & ( RZ <= z_high )
			CCrand = Table()
			CCrand['RA'] = ra_fs[selectionR]
			CCrand['DEC'] = dec_fs[selectionR]
			CCrand['Z'] = z_fs[selectionR]
			N_CCRD = len(CCrand['RA'])
			print(N_CCRD, 'cluster randoms')
			p_2_2PCF_X = os.path.join(validation_dir_WPRP, z_dir, M_str+'_BC_'+LX_str+'_'+str_replic+'-wprp-pimax100-2pcf.fits' )
			print(p_2_2PCF_X)
			if os.path.isfile(p_2_2PCF_X)==False and len(CC['RA'])>10:
				tabulate_wprp_clustering_Cross(
								data['RA'], data['DEC'], data['Z'],
								rand['RA'] , rand['DEC'], rand['Z'],
								CC['RA'], CC['DEC'], CC['Z'],
								CCrand['RA'] , CCrand['DEC'], CCrand['Z'],
								out_file=p_2_2PCF_X )

	#
	# quiescent
	#
	for Ms_min in Ms_mins:
		s1 = (np.log10(GAL['obs_sm'])>Ms_min)&(np.log10(GAL['obs_sfr']/GAL['obs_sm'])<-11)
		data = Table()
		data['RA'] = GAL['RA'][s1]
		data['DEC'] = GAL['DEC'][s1]
		data['Z'] = GAL['redshift_S'][s1]
		M_str = 'Ms_gt_'+str(np.round(Ms_min,2))
		z0 = n.min(data['Z'] )
		z1 = n.max(data['Z'] )
		N_data = len(data)
		print(N_data, 'galaxies')
		x_low, x_high = n.min(GAL['x'][s1]), n.max(GAL['x'][s1])
		y_low, y_high = n.min(GAL['y'][s1]), n.max(GAL['y'][s1])
		z_low, z_high = n.min(GAL['z'][s1]), n.max(GAL['z'][s1])
		size = int(N_data*40)
		uu = n.random.uniform(size=size)
		dec_fs = n.arccos(1 - 2 * uu) * 180 / n.pi - 90.
		ra_fs  = n.random.uniform(size=size) * 2 * n.pi * 180 / n.pi
		z_fs = n.tile(data['Z'], int(size*1./N_data)+1)[:size]
		RX, RY, RZ = get_xyz(ra_fs, dec_fs, z_fs)
		selectionR = ( RX >= x_low ) & ( RX <= x_high ) & ( RY >= y_low ) & ( RY <= y_high ) & ( RZ >= z_low ) & ( RZ <= z_high )
		rand = Table()
		rand['RA'] = ra_fs[selectionR]
		rand['DEC'] = dec_fs[selectionR]
		rand['Z'] = z_fs[selectionR]
		N_RD = len(rand['RA'])
		print(N_RD, 'randoms')
		p_2_2PCF_G = os.path.join(validation_dir_WPRP, z_dir, M_str+'_RS_'+str_replic+'-wprp-pimax100-2pcf.fits' )
		print(p_2_2PCF_G)
		if os.path.isfile(p_2_2PCF_G)==False:
			tabulate_wprp_clustering_noW(data['RA'], data['DEC'], data['Z'], rand['RA'] , rand['DEC'], rand['Z'],  out_file=p_2_2PCF_G )

		for LX_min in LX_mins:
			s2 = ( XGA['CLUSTER_LX_soft_RF_R500c'] > LX_min )
			CC = Table()
			CC['RA'] = XGA['RA'][s2]
			CC['DEC'] = XGA['DEC'][s2]
			CC['Z'] = XGA['redshift_S'][s2]
			N_CLU = len(CC)
			print(N_CLU, 'clusters')
			LX_str = 'HaloLX_gt_'+str(np.round(LX_min,1))
			x_low, x_high = n.min(GAL['x'][s1]), n.max(GAL['x'][s1])
			y_low, y_high = n.min(GAL['y'][s1]), n.max(GAL['y'][s1])
			z_low, z_high = n.min(GAL['z'][s1]), n.max(GAL['z'][s1])
			size = int(N_CLU*40)
			uu = n.random.uniform(size=size)
			dec_fs = n.arccos(1 - 2 * uu) * 180 / n.pi - 90.
			ra_fs  = n.random.uniform(size=size) * 2 * n.pi * 180 / n.pi
			z_fs = n.tile(CC['Z'], int(size*1./N_CLU)+1)[:size]
			RX, RY, RZ = get_xyz(ra_fs, dec_fs, z_fs)
			selectionR = ( RX >= x_low ) & ( RX <= x_high ) & ( RY >= y_low ) & ( RY <= y_high ) & ( RZ >= z_low ) & ( RZ <= z_high )
			CCrand = Table()
			CCrand['RA'] = ra_fs[selectionR]
			CCrand['DEC'] = dec_fs[selectionR]
			CCrand['Z'] = z_fs[selectionR]
			N_CCRD = len(CCrand['RA'])
			print(N_CCRD, 'cluster randoms')
			p_2_2PCF_X = os.path.join(validation_dir_WPRP, z_dir, M_str+'_RS_'+LX_str+'_'+str_replic+'-wprp-pimax100-2pcf.fits' )
			print(p_2_2PCF_X)
			if os.path.isfile(p_2_2PCF_X)==False and len(CC['RA'])>10:
				tabulate_wprp_clustering_Cross(
								data['RA'], data['DEC'], data['Z'],
								rand['RA'] , rand['DEC'], rand['Z'],
								CC['RA'], CC['DEC'], CC['Z'],
								CCrand['RA'] , CCrand['DEC'], CCrand['Z'],
								out_file=p_2_2PCF_X )
