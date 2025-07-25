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

validation_dir       = os.path.join(os.environ['GIT_STMOD_DATA'], 'data', 'validation','validation_AGN')
validation_dir_WPRP = os.path.join(validation_dir, 'WPRP')
os.system('mkdir -p ' + validation_dir_WPRP       )

sys.path.append( os.path.join(os.environ['GIT_STMOD'], 'src') )
from models import AGN as GG

z_dir = sys.argv[1] # 'z0p09'
#LC_dir = sys.argv[2] # 'FullSky'
#C_AGN = GG.GAL(z_dir, LC_dir='FullSky')
LC_dir='FullSky'
C_AGN = GG.AGN(z_dir, LC_dir=LC_dir)
str_scatter_0 = '0.8'
str_fsat = '8.0'

def tabulate_wtheta_clustering_noW(RA, DEC, rand_RA , rand_DEC, out_file ):
	"""
	wtheta direct estimate
	path_2_data : path to the catalogue to correlate
	path_2_random
	"""
	#
	t0 = time.time()
	N = len(RA)
	rand_N = len(rand_RA)
	print(N, rand_N, out_file, time.time()-t0)
	bins = 10**n.arange(-3.0, 0.65, 0.2)
	nbins = len(bins)-1
	#print('bins', bins, bins.shape)
	x = (bins[1:]+bins[:-1])/2.
	cosmology = 2
	nthreads = 16
	# Auto pairs counts in DD
	autocorr=1
	DD_counts = DDtheta_mocks(autocorr, nthreads, bins, n.array(list(RA)).astype('float'), n.array(list(DEC)).astype('float'))
	#print(DD_counts)
	autocorr=0
	DR_counts = DDtheta_mocks(autocorr, nthreads, bins,
							n.array(list(RA)).astype('float'),
							n.array(list(DEC)).astype('float'),
							RA2=rand_RA.astype('float'),
							DEC2=rand_DEC.astype('float'))
	# Auto pairs counts in RR
	autocorr=1
	RR_counts = DDtheta_mocks(autocorr, nthreads, bins,
							rand_RA.astype('float'),
							rand_DEC.astype('float'))
	# All the pair counts are done, get the angular correlation function
	wtheta = convert_3d_counts_to_cf(N, N, rand_N, rand_N, DD_counts, DR_counts, DR_counts, RR_counts)
	t = Table()
	t.add_column(Column(data = bins[:-1], name='theta_min', unit='deg'  ) )
	t.add_column(Column(data = bins[1:], name='theta_max', unit='deg'  ) )
	t.add_column(Column(data = x, name='theta_mid', unit='deg'  ) )
	t.add_column(Column(data = wtheta, name='wtheta', unit=''  ) )
	t.add_column(Column(data = n.ones_like(x) * N, name='N_data', unit=''  ) )
	t.add_column(Column(data = n.ones_like(x) * rand_N, name='N_random', unit=''  ) )
	t.add_column(Column(data = DD_counts['npairs'], name='DD_counts', unit=''  ) )
	t.add_column(Column(data = DR_counts['npairs'], name='DR_counts', unit=''  ) )
	t.add_column(Column(data = RR_counts['npairs'], name='RR_counts', unit=''  ) )
	#print(t.info())
	#print(t)
	#repeat when removing random 10%
	wtheta_JK = n.zeros((50, len(wtheta)))
	for jj in n.arange(50):
		s1=(n.random.random(len(RA))<0.9)
		ra1, dec1 = n.array(list(RA)).astype('float')[s1], n.array(list(DEC)).astype('float')[s1]
		# Auto pairs counts in DD
		autocorr=1
		DD_counts = DDtheta_mocks(autocorr, nthreads, bins, ra1, dec1)
		# Auto pairs counts in DR
		autocorr=0
		DR_counts = DDtheta_mocks(autocorr, nthreads, bins, ra1, dec1, RA2=rand_RA.astype('float'), DEC2=rand_DEC.astype('float'))
		# All the pair counts are done, get the angular correlation function
		wtheta_JK[jj] = convert_3d_counts_to_cf(N, N, rand_N, rand_N, DD_counts, DR_counts, DR_counts, RR_counts)
	t.add_column(Column(data = wtheta_JK.mean(axis=0), name='wtheta_JK_mean', unit=''  ) )
	t.add_column(Column(data = wtheta_JK.std(axis=0), name='wtheta_JK_std', unit=''  ) )
	t.write(out_file, overwrite=True)
	print(out_file, 'written', time.time()-t0, 's')


for meta in C_AGN.LC_MetaData:
	#
	print(meta)
	# retrieve the resulting catalogues and meta data
	str_replic = 'replication_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz'])
	p_2_catalogue = os.path.join(os.environ['UCHUU'], LC_dir, z_dir, 'replication_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']), 'glist.fits')
	p_2_catal_MAG = os.path.join(os.environ['UCHUU'], LC_dir, z_dir, 'replication_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']), 'Kmatch_mags.fits')
	GAL = Table.read(p_2_catalogue)
	#MAG = Table.read(p_2_catal_MAG)

	#p_2_catal_AGNs = np.array( glob.glob( os.path.join( os.path.dirname(p_2_catalogue), 'AGN_list_sigma_*_fsat_*.fits' ) ) )
	p_2_catal_AGNs = np.array( glob.glob( os.path.join( os.path.dirname(p_2_catalogue), 'AGN_list_sigma_'+str_scatter_0+'_fsat_'+str_fsat+'.fits' ) ) )
	print(not os.path.isfile(p_2_catalogue))
	print(not os.path.isfile(p_2_catal_MAG))
	print(len(p_2_catal_AGNs)==0)
	if not os.path.isfile(p_2_catalogue) or not os.path.isfile(p_2_catal_MAG) or len(p_2_catal_AGNs)==0 :
		continue
	#AGNs = {}
	for p_2_catal_AGN in p_2_catal_AGNs:
		AGN = Table.read(p_2_catal_AGN)
		bn_agn = os.path.basename(p_2_catal_AGN)[:-5]
		#AGNs[AGN_cat_names[p_2_catal_AGN]] = t_i[(t_i['redshift_S']>=z_min)&(t_i['redshift_S']<z_max)]
		#z_min, z_max = np.min(GAL['redshift_S']), np.max(GAL['redshift_S'])
		#print('z_min, z_max=', z_min, z_max)
		#volume_mock = (( cosmo.comoving_volume(z_max) - cosmo.comoving_volume(z_min) ) * meta['mean_area'] * np.pi / 129600.).value
		#z_mean = np.mean(GAL['redshift_S'])
		#z_mins = n.array([0.05])#, 0.1, 0.2, 0.3, 0.4])
		#z_maxs = n.array([0.5 ])#, 0.2, 0.3, 0.4, 0.5])
		magr_maxs = n.array([22.5])
		r_max = 22.5
		fx_mins = n.array([-14.3, -14., -13.5, -13.3])
		#for r_max in magr_maxs:
		for fx_min in fx_mins:
			s_Z = ( AGN['SDSS_r_AB']<=r_max ) & ( AGN['FX_soft'] > fx_min )
			data = Table()
			data['RA'] = GAL['RA'][AGN['ID_glist'][s_Z]]
			data['DEC'] = GAL['DEC'][AGN['ID_glist'][s_Z]]
			data['Z'] = AGN['redshift_S'][s_Z]
			N_data = len(data)
			x_low, x_high = n.min(GAL['x'][AGN['ID_glist'][s_Z]]), n.max(GAL['x'][AGN['ID_glist'][s_Z]])
			y_low, y_high = n.min(GAL['y'][AGN['ID_glist'][s_Z]]), n.max(GAL['y'][AGN['ID_glist'][s_Z]])
			z_low, z_high = n.min(GAL['z'][AGN['ID_glist'][s_Z]]), n.max(GAL['z'][AGN['ID_glist'][s_Z]])

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
			##for jj, p2_rand in enumerate(p2_rands):
			#print(z_min, z_max, N_data,'data and ', N_RD,'randoms opened')
			#p_2_data = os.path.join(validation_dir_WPRP, 'LS_DR10_8R20_'+str(z_min) +'_zPHOT_' + str(z_max) +'-Nmax-R0-DATA.fits' )
			#data.write(p_2_data, format='fits', overwrite = True)
			#p_2_rand = os.path.join(validation_dir_WPRP, 'LS_DR10_8R20_'+str(z_min) +'_zPHOT_' + str(z_max) +'-Nmax-R0-RAND.fits' )
			#rand.write(p_2_rand, format='fits', overwrite = True)
			#
			p_2_2PCF = os.path.join(validation_dir_WPRP, z_dir+'_'+str_replic+'_'+bn_agn+'_rmag_' + str(r_max)+'_fx_' + str(fx_min) +'-wtheta-2pcf.fits' )
			#if os.path.isfile(p_2_2PCF)==False:
			tabulate_wtheta_clustering_noW(data['RA'], data['DEC'], rand['RA'] , rand['DEC'],  out_file=p_2_2PCF)
