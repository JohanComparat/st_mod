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




validation_dir       = os.path.join(os.environ['GIT_STMOD'], 'data', 'validation','validation_GAL')
validation_dir_WPRP = os.path.join(validation_dir, 'WPRP')
os.system('mkdir -p ' + validation_dir_WPRP       )

sys.path.append( os.path.join(os.environ['GIT_STMOD'], 'src') )
from models import GAL as GG

z_dir = sys.argv[1]
#LC_dir = sys.argv[2] # 'FullSky'
C_GAL = GG.GAL(z_dir, LC_dir='FullSky')
LC_dir='FullSky'

def tabulate_wprp_clustering_noW(RA, DEC, Z, rand_RA , rand_DEC, rand_Z, out_file='test.fits', CV_frac=0.01, pimax = 100.0, N_JK=20 ):
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
	t.add_column(Column(data = np.ones_like(x) * CV_frac, name='CV_frac', unit=''  ) )
	# repeat when removing random 10%
	wprp_JK = np.zeros((N_JK, len(wp)))
	for jj in np.arange(N_JK):
		s1=(np.random.random(len(RA))<0.9)
		ra1, dec1, cz1 = np.array(list(RA)).astype('float')[s1], np.array(list(DEC)).astype('float')[s1], np.array(list(CZ)).astype('float')[s1]
		# Auto pairs counts in DD
		autocorr=1
		DD_counts = DDrppi_mocks(autocorr, cosmology, nthreads, pimax, bins,  ra1, dec1, cz1)
		# Auto pairs counts in DR
		autocorr=0
		DR_counts = DDrppi_mocks(autocorr, cosmology, nthreads, pimax, bins, ra1, dec1, cz1, RA2=rand_RA.astype('float'), DEC2=rand_DEC.astype('float'), CZ2=rand_CZ.astype('float'))
		# Auto pairs counts in RR
		autocorr=1
		RR_counts = DDrppi_mocks(autocorr, cosmology, nthreads, pimax, bins, rand_RA.astype('float'), rand_DEC.astype('float'), rand_CZ.astype('float'))
		# All the pair counts are done, get the angular correlation function
		wprp_JK[jj] = convert_rp_pi_counts_to_wp(N, N, rand_N, rand_N, DD_counts, DR_counts, DR_counts, RR_counts, nbins, pimax)
	#print ( "wprp_JK.mean(axis=0)", wprp_JK.mean(axis=0).shape, wprp_JK.mean(axis=0),  wprp_JK.mean(axis=0) / wp  )
	#print ( "wprp_JK.std(axis=0) ", wprp_JK.std(axis=0) .shape, wprp_JK.std(axis=0),  wprp_JK.std(axis=0) / wp  )
	#print ( "wp                ", wp                .shape, wp                  )
	#t['wprp_JK'] = wprp_JK
	t.add_column(Column(data = wprp_JK.mean(axis=0), name='wprp_JK_mean', unit=''  ) )
	t.add_column(Column(data = wprp_JK.std(axis=0), name='wprp_JK_std', unit=''  ) )
	print(out_file)
	t.write(out_file, overwrite=True, format='fits')
	print(out_file, time.time()-t0, 's')


#enough_area = (C_GAL.LC_MetaData['area_DC_max']>=0.5*np.max(C_GAL.LC_MetaData['area_DC_max'])) & (C_GAL.LC_MetaData['area_DC_max']>0)
#C_GAL.LC_MetaData['mean_area'] = (C_GAL.LC_MetaData['area_DC_min'] + C_GAL.LC_MetaData['area_DC_max']) / 2.
#small_difference_minmax_1 = ( C_GAL.LC_MetaData['area_DC_min'] / C_GAL.LC_MetaData['mean_area'] >= 0.8 ) & ( C_GAL.LC_MetaData['area_DC_min'] / C_GAL.LC_MetaData['mean_area'] <= 1.2 )
#small_difference_minmax_2 = ( C_GAL.LC_MetaData['area_DC_max'] / C_GAL.LC_MetaData['mean_area'] >= 0.8 ) & ( C_GAL.LC_MetaData['area_DC_max'] / C_GAL.LC_MetaData['mean_area'] <= 1.2 )
#for meta in C_GAL.LC_MetaData:#[(enough_area)&(small_difference_minmax_1)&(small_difference_minmax_2)]:
	##
	#print(meta)
	## retrieve the resulting catalogues and meta data
	#str_replic = 'replication_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz'])
	#p_2_catalogue = os.path.join(os.environ['UCHUU'], LC_dir, z_dir, 'replication_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']), 'glist.fits')
	#p_2_catal_MAG = os.path.join(os.environ['UCHUU'], LC_dir, z_dir, 'replication_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']), 'Kmatch_mags.fits')
	#GAL = Table.read(p_2_catalogue)
	#MAG = Table.read(p_2_catal_MAG)
	##z_min, z_max = np.min(GAL['redshift_S']), np.max(GAL['redshift_S'])
	##print('z_min, z_max=', z_min, z_max)
	##volume_mock = (( cosmo.comoving_volume(z_max) - cosmo.comoving_volume(z_min) ) * meta['mean_area'] * np.pi / 129600.).value
	##z_mean = np.mean(GAL['redshift_S'])
	#z_mins = n.array([0.05])#, 0.1, 0.2, 0.3, 0.4])
	#z_maxs = n.array([0.5 ])#, 0.2, 0.3, 0.4, 0.5])
	#magr_maxs = n.array([17.77, 18, 19, 19.5, 20, 21])
	#for r_max in magr_maxs:
		#s_Z = ( MAG['rmag']<=r_max )
		##s1 = ( n.random.random(len(data_i)) < 4e5 / len(data_i) )
		#data = Table()
		#data['RA'] = GAL['RA'][s_Z]
		#data['DEC'] = GAL['DEC'][s_Z]
		#data['Z'] = GAL['redshift_S'][s_Z]
		#N_data = len(data)
		#x_low, x_high = n.min(GAL['x']), n.max(GAL['x'])
		#y_low, y_high = n.min(GAL['y']), n.max(GAL['y'])
		#z_low, z_high = n.min(GAL['z']), n.max(GAL['z'])

		#size = int(N_data*40)
		#uu = n.random.uniform(size=size)
		#dec_fs = n.arccos(1 - 2 * uu) * 180 / n.pi - 90.
		#ra_fs  = n.random.uniform(size=size) * 2 * n.pi * 180 / n.pi
		#z_fs = n.tile(data['Z'], int(size*1./N_data)+1)[:size]
		#RX, RY, RZ = get_xyz(ra_fs, dec_fs, z_fs)
		#selectionR = ( RX >= x_low ) & ( RX <= x_high ) & ( RY >= y_low ) & ( RY <= y_high ) & ( RZ >= z_low ) & ( RZ <= z_high )
		#rand = Table()
		#rand['RA'] = ra_fs[selectionR]
		#rand['DEC'] = dec_fs[selectionR]
		#rand['Z'] = z_fs[selectionR]
		#N_RD = len(rand['RA'])
		###for jj, p2_rand in enumerate(p2_rands):
		##print(z_min, z_max, N_data,'data and ', N_RD,'randoms opened')
		##p_2_data = os.path.join(validation_dir_WPRP, 'LS_DR10_8R20_'+str(z_min) +'_zPHOT_' + str(z_max) +'-Nmax-R0-DATA.fits' )
		##data.write(p_2_data, format='fits', overwrite = True)
		##p_2_rand = os.path.join(validation_dir_WPRP, 'LS_DR10_8R20_'+str(z_min) +'_zPHOT_' + str(z_max) +'-Nmax-R0-RAND.fits' )
		##rand.write(p_2_rand, format='fits', overwrite = True)
		##
		#p_2_2PCF = os.path.join(validation_dir_WPRP, z_dir+'_'+str_replic+'_rmag_' + str(r_max) +'-wprp-pimax60-2pcf.fits' )
		#if os.path.isfile(p_2_2PCF)==False:
			#tabulate_wprp_clustering_noW(data['RA'], data['DEC'], data['Z'], rand['RA'] , rand['DEC'], rand['Z'],  out_file=p_2_2PCF, CV_frac=0.01, pimax = 60.0, N_JK = 20 )
		#p_2_2PCF = os.path.join(validation_dir_WPRP, z_dir+'_'+str_replic+'_rmag_' + str(r_max) +'-wprp-pimax40-2pcf.fits' )
		#if os.path.isfile(p_2_2PCF)==False:
			#tabulate_wprp_clustering_noW(data['RA'], data['DEC'], data['Z'], rand['RA'] , rand['DEC'], rand['Z'],  out_file=p_2_2PCF, CV_frac=0.01, pimax = 40.0, N_JK = 20 )

MsBds= n.array([9
		,9.5
		,10
		,10.5
		,10.75
		,11.0
		,11.25
		,11.5 , 12.0 ])
MsMins = MsBds[:-1]
MsMaxs = MsBds[1:]

z2psf = n.array([ 0.07
		, 0.09
		, 0.12
		, 0.17
		, 0.21
		, 0.33
		, 0.5
		, 0.5  ])

for meta in C_GAL.LC_MetaData:#[(enough_area)&(small_difference_minmax_1)&(small_difference_minmax_2)]:
	#
	print(meta)
	# retrieve the resulting catalogues and meta data
	str_replic = 'replication_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz'])
	p_2_catalogue = os.path.join(os.environ['UCHUU'], LC_dir, z_dir, 'replication_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']), 'glist.fits')
	p_2_catal_MAG = os.path.join(os.environ['UCHUU'], LC_dir, z_dir, 'replication_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']), 'Kmatch_mags.fits')
	GAL = Table.read(p_2_catalogue)
	#MAG = Table.read(p_2_catal_MAG)
	#z_min, z_max = np.min(GAL['redshift_S']), np.max(GAL['redshift_S'])
	#print('z_min, z_max=', z_min, z_max)
	#volume_mock = (( cosmo.comoving_volume(z_max) - cosmo.comoving_volume(z_min) ) * meta['mean_area'] * np.pi / 129600.).value
	#z_mean = np.mean(GAL['redshift_S'])
	#z_mins = n.array([0.05])#, 0.1, 0.2, 0.3, 0.4])
	#z_maxs = n.array([0.5 ])#, 0.2, 0.3, 0.4, 0.5])
	#magr_maxs = n.array([17.77, 18, 19, 19.5, 20, 21])
	for m0, m1, z_max in zip(MsMins, MsMaxs, z2psf):
		M_str = 'Ms_gt_'+str(m0)
		p_2_2PCF = os.path.join(validation_dir_WPRP, z_dir+'_'+M_str+'-wprp-pimax40-2pcf.fits' )
		if os.path.isfile(p_2_2PCF)==False:
			s_Z = ( n.log10(GAL['obs_sm'])>=m0 )
			z0 = n.min(GAL['redshift_S'][s_Z])
			z1 = n.max(GAL['redshift_S'][s_Z])
			data = Table()
			data['RA'] = GAL['RA'][s_Z]
			data['DEC'] = GAL['DEC'][s_Z]
			data['Z'] = GAL['redshift_S'][s_Z]
			N_data = len(data)
			if N_data>100:
				x_low, x_high = n.min(GAL['x']), n.max(GAL['x'])
				y_low, y_high = n.min(GAL['y']), n.max(GAL['y'])
				z_low, z_high = n.min(GAL['z']), n.max(GAL['z'])

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

				tabulate_wprp_clustering_noW(data['RA'], data['DEC'], data['Z'], rand['RA'] , rand['DEC'], rand['Z'],  out_file=p_2_2PCF, CV_frac=0.01, pimax = 40.0, N_JK = 10 )


		z_errs = n.array([0.005, 0.01, 0.02, 0.03, 0.05])
		for z_err in z_errs:
			out_str = 'Ms_gt_'+str(m0)+'_dz_'+str(int(1000*z_err)).zfill(4)
			p_2_2PCF = os.path.join(validation_dir_WPRP, z_dir+'_'+out_str+'-wprp-pimax40-2pcf.fits' )
			if os.path.isfile(p_2_2PCF)==False:
				s_Z = ( n.log10(GAL['obs_sm'])>=m0 )
				z0 = n.min(GAL['redshift_S'][s_Z])
				z1 = n.max(GAL['redshift_S'][s_Z])
				z_true = GAL['redshift_S'][s_Z]
				new_z = norm.rvs(loc = z_true, scale = (1+z_true)*z_err, size=len(z_true))
				s_Z2 = ( new_z>=z0 ) & ( new_z <=z1 )
				data = Table()
				data['RA'] = GAL['RA'][s_Z][s_Z2]
				data['DEC'] = GAL['DEC'][s_Z][s_Z2]
				data['Z'] = new_z[s_Z2]
				N_data = len(data)
				if N_data>100:

					x_low, x_high = n.min(GAL['x']), n.max(GAL['x'])
					y_low, y_high = n.min(GAL['y']), n.max(GAL['y'])
					z_low, z_high = n.min(GAL['z']), n.max(GAL['z'])

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

					tabulate_wprp_clustering_noW(data['RA'], data['DEC'], data['Z'], rand['RA'] , rand['DEC'], rand['Z'],  out_file=p_2_2PCF, CV_frac=0.01, pimax = 40.0, N_JK = 10 )
