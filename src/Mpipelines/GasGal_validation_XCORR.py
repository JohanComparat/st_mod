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

validation_dir       = os.path.join(os.environ['GIT_STMOD_DATA'], 'data', 'validation','validation_GasGal')
validation_dir_WPRP = os.path.join(validation_dir, 'XCORR')
os.system('mkdir -p ' + validation_dir_WPRP       )


def get_wtheta(RA, DEC, rand_RA , rand_DEC ):
	"""
	wtheta direct estimate
	path_2_data : path to the catalogue to correlate
	path_2_random
	"""
	#
	t0 = time.time()
	N = len(RA)
	rand_N = len(rand_RA)
	print(N, rand_N, time.time()-t0)
	bins = 10**np.arange(-3.0, 1.4, 0.2)
	nbins = len(bins)-1
	#print('bins', bins, bins.shape)
	x = (bins[1:]+bins[:-1])/2.
	cosmology = 2
	nthreads = 16
	# Auto pairs counts in DD
	autocorr=1
	DD_counts = DDtheta_mocks(autocorr, nthreads, bins, np.array(list(RA)).astype('float'), np.array(list(DEC)).astype('float'))
	#print(DD_counts)
	autocorr=0
	DR_counts = DDtheta_mocks(autocorr, nthreads, bins,
													np.array(list(RA)).astype('float'),
													np.array(list(DEC)).astype('float'),
													RA2=rand_RA.astype('float'),
													DEC2=rand_DEC.astype('float'))
	# Auto pairs counts in RR
	autocorr=1
	RR_counts = DDtheta_mocks(autocorr, nthreads, bins,
													rand_RA.astype('float'),
													rand_DEC.astype('float'))
	# All the pair counts are done, get the angular correlation function
	wtheta = convert_3d_counts_to_cf(N, N, rand_N, rand_N, DD_counts, DR_counts, DR_counts, RR_counts)
	return wtheta

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
	bins = 10**np.arange(-3.0, 1.4, 0.2)
	nbins = len(bins)-1
	#print('bins', bins, bins.shape)
	x = (bins[1:]+bins[:-1])/2.
	cosmology = 2
	nthreads = 16
	# Auto pairs counts in DD
	autocorr=1
	DD_counts = DDtheta_mocks(autocorr, nthreads, bins, np.array(list(RA)).astype('float'), np.array(list(DEC)).astype('float'))
	#print(DD_counts)
	autocorr=0
	DR_counts = DDtheta_mocks(autocorr, nthreads, bins,
							np.array(list(RA)).astype('float'),
							np.array(list(DEC)).astype('float'),
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
	t.add_column(Column(data = np.ones_like(x) * N, name='N_data', unit=''  ) )
	t.add_column(Column(data = np.ones_like(x) * rand_N, name='N_random', unit=''  ) )
	t.add_column(Column(data = DD_counts['npairs'], name='DD_counts', unit=''  ) )
	t.add_column(Column(data = DR_counts['npairs'], name='DR_counts', unit=''  ) )
	t.add_column(Column(data = RR_counts['npairs'], name='RR_counts', unit=''  ) )

	t.write(out_file, overwrite=True)
	print(out_file, 'written', time.time()-t0, 's')


def tabulate_XCORRwtheta_clustering_noW(RA, DEC, RA2, DEC2, rand_RA , rand_DEC, out_file ):
	"""
	wtheta direct estimate
	path_2_data : path to the catalogue to correlate
	path_2_random
	"""
	#
	t0 = time.time()
	N = len(RA)
	N2 = len(RA2)
	rand_N = len(rand_RA)
	print(N, N2, rand_N, out_file, time.time()-t0)
	bins = 10**np.arange(-3.0, 1.4, 0.2)
	nbins = len(bins)-1
	#print('bins', bins, bins.shape)
	x = (bins[1:]+bins[:-1])/2.
	cosmology = 2
	nthreads = 16
	# Auto pairs counts in DD
	autocorr=1
	DD_counts = DDtheta_mocks(autocorr, nthreads, bins, np.array(list(RA)).astype('float'), np.array(list(DEC)).astype('float'))
	#print(DD_counts)
	autocorr=0
	D1D2_counts = DDtheta_mocks(autocorr, nthreads, bins,
							np.array(list(RA)).astype('float'),
							np.array(list(DEC)).astype('float'),
							RA2=np.array(list(RA2)).astype('float'),
							DEC2=np.array(list(DEC2)).astype('float'))

	D1R_counts = DDtheta_mocks(autocorr, nthreads, bins,
							np.array(list(RA)).astype('float'),
							np.array(list(DEC)).astype('float'),
							RA2=rand_RA.astype('float'),
							DEC2=rand_DEC.astype('float'))

	D2R_counts = DDtheta_mocks(autocorr, nthreads, bins,
							np.array(list(RA2)).astype('float'),
							np.array(list(DEC2)).astype('float'),
							RA2=rand_RA.astype('float'),
							DEC2=rand_DEC.astype('float'))
	# Auto pairs counts in RR
	autocorr=1
	RR_counts = DDtheta_mocks(autocorr, nthreads, bins,
							rand_RA.astype('float'),
							rand_DEC.astype('float'))
	# All the pair counts are done, get the angular correlation function
	wtheta = convert_3d_counts_to_cf(N, N2, rand_N, rand_N, D1D2_counts, D1R_counts, D2R_counts, RR_counts)
	t = Table()
	t.add_column(Column(data = bins[:-1], name='theta_min', unit='deg'  ) )
	t.add_column(Column(data = bins[1:], name='theta_max', unit='deg'  ) )
	t.add_column(Column(data = x, name='theta_mid', unit='deg'  ) )
	t.add_column(Column(data = wtheta, name='wtheta', unit=''  ) )
	t.add_column(Column(data = np.ones_like(x) * N, name='N_data', unit=''  ) )
	t.add_column(Column(data = np.ones_like(x) * N2, name='N2_data', unit=''  ) )
	t.add_column(Column(data = np.ones_like(x) * rand_N, name='N_random', unit=''  ) )
	t.add_column(Column(data = D1D2_counts['npairs'], name='D1D2_counts', unit=''  ) )
	t.add_column(Column(data = D1R_counts['npairs'], name='D1R_counts', unit=''  ) )
	t.add_column(Column(data = D2R_counts['npairs'], name='D2R_counts', unit=''  ) )
	t.add_column(Column(data = RR_counts['npairs'], name='RR_counts', unit=''  ) )

	t.write(out_file, overwrite=True)
	print(out_file, 'written', time.time()-t0, 's')

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

agn_seed = '1' # sys.argv[1] # 1
clu_seed = '1' # sys.argv[2] # 1
LC_dir = 'LCerass'
top_dir = os.path.join(os.environ['UCHUU'], LC_dir)
nl = lambda sel : len(sel.nonzero()[0])

sky_map_hdu = Table.read(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS.fits') )

LC_dir = 'LCerass'
top_dir = os.path.join(os.environ['UCHUU'], LC_dir)

def get_srvmap(ra, dec):
    return sky_map_hdu['SRVMAP'].value[(sky_map_hdu['RA_MIN']<ra ) & ( sky_map_hdu['RA_MAX'] >= ra ) & ( sky_map_hdu['DE_MIN']<dec ) & ( sky_map_hdu['DE_MAX'] >= dec)]

def get_srvmap_rev(ra, dec, sky_tile):
    return (sky_tile['RA_MIN']<ra ) & ( sky_tile['RA_MAX'] >= ra ) & ( sky_tile['DE_MIN']<dec ) & ( sky_tile['DE_MAX'] >= dec)


size = int(50e6)
uu = np.random.uniform(size=size)
dec_fs = np.arccos(1 - 2 * uu) * 180 / np.pi - 90.
ra_fs = np.random.uniform(size=size) * 2 * np.pi * 180 / np.pi

for sky_tile in sky_map_hdu[(sky_map_hdu['OWNER']==2)|(sky_map_hdu['OWNER']==0)]:
	sky_tile_id = str(sky_tile['SRVMAP'])
	str_field = str(sky_tile['SRVMAP']).zfill(6)
	print(str_field)
	esass_dir = os.path.join(os.environ['UCHUU'], LC_dir, str_field,
							 'GE_e4_merge_AGNseed' + agn_seed.zfill(3) + '_SimBKG_CLUseed' + clu_seed.zfill(3))

	dir_2pcf = os.path.join(esass_dir, 'XCORR')
	os.system('mkdir -p ' + dir_2pcf)
	path_2_event_file = os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'glist.fits')
	path_2_simeventAGN_file = os.path.join(esass_dir, 'simAGNevt_'+str_field+'.fits')
	path_2_simeventCLU_file = os.path.join(esass_dir, 'simCLUevt_'+str_field+'.fits')
	path_2_simeventBKG_file = os.path.join(esass_dir, 'simBKGevt_'+str_field+'.fits')
	s_R = get_srvmap_rev(ra_fs, dec_fs, sky_tile)
	basename = 'GAL_m' + str(np.round(11.0, 1)) + '_evAGNevCLUevBKG'
	p_2_2PCF = os.path.join(dir_2pcf, basename + '_CROSSCORR_05E20.wtheta.2pcf.fits')

	if os.path.isfile(path_2_event_file) and os.path.isfile(path_2_simeventBKG_file) and os.path.isfile(path_2_simeventCLU_file) and os.path.isfile(path_2_simeventAGN_file) and not os.path.isfile(p_2_2PCF):
		GAL = Table.read(path_2_event_file)
		evBKG = Table.read(path_2_simeventBKG_file)
		evBKG = evBKG[(evBKG['SIGNAL']>=0.5)&(evBKG['SIGNAL']<=2.0)]
		evAGN = Table.read(path_2_simeventAGN_file)
		evAGN = evAGN[(evAGN['SIGNAL']>=0.5)&(evAGN['SIGNAL']<=2.0)]
		evCLU = Table.read(path_2_simeventCLU_file)
		evBKG = evBKG[(evBKG['SIGNAL']>=0.5)&(evBKG['SIGNAL']<=2.0)]
		m0 = 10.0
		z1 = 0.18
		s10 = (np.log10(GAL['obs_sm'])>=m0) & (GAL['redshift_S']>0.05) & (GAL['redshift_S']<z1)
		basename = 'GAL_m'+str(np.round(m0,1))+'_evGAS'
		p_2_2PCF = os.path.join(dir_2pcf, basename +'_CROSSCORR_05E20.wtheta.2pcf.fits' )
		if not os.path.isfile(p_2_2PCF):
			tabulate_XCORRwtheta_clustering_noW(evCLU['RA'], evCLU['DEC'],
											GAL['RA'][s10], GAL['DEC'][s10],
											ra_fs[s_R], dec_fs[s_R], p_2_2PCF)
		basename = 'GAL_m'+str(np.round(m0,1))+'_evBKG'
		p_2_2PCF = os.path.join(dir_2pcf, basename +'_CROSSCORR_05E20.wtheta.2pcf.fits' )
		if not os.path.isfile(p_2_2PCF):
			tabulate_XCORRwtheta_clustering_noW(evBKG['RA'], evBKG['DEC'],
											GAL['RA'][s10], GAL['DEC'][s10],
											ra_fs[s_R], dec_fs[s_R], p_2_2PCF)
		basename = 'GAL_m'+str(np.round(m0,1))+'_evAGN'
		p_2_2PCF = os.path.join(dir_2pcf, basename +'_CROSSCORR_05E20.wtheta.2pcf.fits' )
		if not os.path.isfile(p_2_2PCF):
			tabulate_XCORRwtheta_clustering_noW(evAGN['RA'], evAGN['DEC'],
											GAL['RA'][s10], GAL['DEC'][s10],
											ra_fs[s_R], dec_fs[s_R], p_2_2PCF)
		basename = 'GAL_m'+str(np.round(m0,1))+'_evAGNevCLU'
		p_2_2PCF = os.path.join(dir_2pcf, basename +'_CROSSCORR_05E20.wtheta.2pcf.fits' )
		if not os.path.isfile(p_2_2PCF):
			tabulate_XCORRwtheta_clustering_noW(np.hstack((evCLU['RA'], evAGN['RA'])), np.hstack((evCLU['DEC'],evAGN['DEC'])),
											GAL['RA'][s10], GAL['DEC'][s10],
											ra_fs[s_R], dec_fs[s_R], p_2_2PCF)
		basename = 'GAL_m'+str(np.round(m0,1))+'_evAGNevCLUevBKG'
		p_2_2PCF = os.path.join(dir_2pcf, basename +'_CROSSCORR_05E20.wtheta.2pcf.fits' )
		if not os.path.isfile(p_2_2PCF):
			tabulate_XCORRwtheta_clustering_noW(np.hstack((evCLU['RA'], evAGN['RA'], evBKG['RA'])), np.hstack((evCLU['DEC'], evAGN['DEC'], evBKG['DEC'])),
											GAL['RA'][s10], GAL['DEC'][s10],
											ra_fs[s_R], dec_fs[s_R], p_2_2PCF)

		m0 = 10.5
		z1 = 0.26
		s10 = (np.log10(GAL['obs_sm'])>=m0) & (GAL['redshift_S']>0.05) & (GAL['redshift_S']<z1)
		basename = 'GAL_m'+str(np.round(m0,1))+'_evGAS'
		p_2_2PCF = os.path.join(dir_2pcf, basename +'_CROSSCORR_05E20.wtheta.2pcf.fits' )
		if not os.path.isfile(p_2_2PCF):
			tabulate_XCORRwtheta_clustering_noW(evCLU['RA'], evCLU['DEC'],
											GAL['RA'][s10], GAL['DEC'][s10],
											ra_fs[s_R], dec_fs[s_R], p_2_2PCF)
		basename = 'GAL_m'+str(np.round(m0,1))+'_evBKG'
		p_2_2PCF = os.path.join(dir_2pcf, basename +'_CROSSCORR_05E20.wtheta.2pcf.fits' )
		if not os.path.isfile(p_2_2PCF):
			tabulate_XCORRwtheta_clustering_noW(evBKG['RA'], evBKG['DEC'],
											GAL['RA'][s10], GAL['DEC'][s10],
											ra_fs[s_R], dec_fs[s_R], p_2_2PCF)
		basename = 'GAL_m'+str(np.round(m0,1))+'_evAGN'
		p_2_2PCF = os.path.join(dir_2pcf, basename +'_CROSSCORR_05E20.wtheta.2pcf.fits' )
		if not os.path.isfile(p_2_2PCF):
			tabulate_XCORRwtheta_clustering_noW(evAGN['RA'], evAGN['DEC'],
											GAL['RA'][s10], GAL['DEC'][s10],
											ra_fs[s_R], dec_fs[s_R], p_2_2PCF)
		basename = 'GAL_m'+str(np.round(m0,1))+'_evAGNevCLU'
		p_2_2PCF = os.path.join(dir_2pcf, basename +'_CROSSCORR_05E20.wtheta.2pcf.fits' )
		if not os.path.isfile(p_2_2PCF):
			tabulate_XCORRwtheta_clustering_noW(np.hstack((evCLU['RA'], evAGN['RA'])), np.hstack((evCLU['DEC'],evAGN['DEC'])),
											GAL['RA'][s10], GAL['DEC'][s10],
											ra_fs[s_R], dec_fs[s_R], p_2_2PCF)
		basename = 'GAL_m'+str(np.round(m0,1))+'_evAGNevCLUevBKG'
		p_2_2PCF = os.path.join(dir_2pcf, basename +'_CROSSCORR_05E20.wtheta.2pcf.fits' )
		if not os.path.isfile(p_2_2PCF):
			tabulate_XCORRwtheta_clustering_noW(np.hstack((evCLU['RA'], evAGN['RA'], evBKG['RA'])), np.hstack((evCLU['DEC'], evAGN['DEC'], evBKG['DEC'])),
											GAL['RA'][s10], GAL['DEC'][s10],
											ra_fs[s_R], dec_fs[s_R], p_2_2PCF)



		m0 = 11.0
		z1 = 0.35
		s10 = (np.log10(GAL['obs_sm'])>=m0) & (GAL['redshift_S']>0.05) & (GAL['redshift_S']<z1)
		basename = 'GAL_m'+str(np.round(m0,1))+'_evGAS'
		p_2_2PCF = os.path.join(dir_2pcf, basename +'_CROSSCORR_05E20.wtheta.2pcf.fits' )
		if not os.path.isfile(p_2_2PCF):
			tabulate_XCORRwtheta_clustering_noW(evCLU['RA'], evCLU['DEC'],
											GAL['RA'][s10], GAL['DEC'][s10],
											ra_fs[s_R], dec_fs[s_R], p_2_2PCF)
		basename = 'GAL_m'+str(np.round(m0,1))+'_evBKG'
		p_2_2PCF = os.path.join(dir_2pcf, basename +'_CROSSCORR_05E20.wtheta.2pcf.fits' )
		if not os.path.isfile(p_2_2PCF):
			tabulate_XCORRwtheta_clustering_noW(evBKG['RA'], evBKG['DEC'],
											GAL['RA'][s10], GAL['DEC'][s10],
											ra_fs[s_R], dec_fs[s_R], p_2_2PCF)
		basename = 'GAL_m'+str(np.round(m0,1))+'_evAGN'
		p_2_2PCF = os.path.join(dir_2pcf, basename +'_CROSSCORR_05E20.wtheta.2pcf.fits' )
		if not os.path.isfile(p_2_2PCF):
			tabulate_XCORRwtheta_clustering_noW(evAGN['RA'], evAGN['DEC'],
											GAL['RA'][s10], GAL['DEC'][s10],
											ra_fs[s_R], dec_fs[s_R], p_2_2PCF)
		basename = 'GAL_m'+str(np.round(m0,1))+'_evAGNevCLU'
		p_2_2PCF = os.path.join(dir_2pcf, basename +'_CROSSCORR_05E20.wtheta.2pcf.fits' )
		if not os.path.isfile(p_2_2PCF):
			tabulate_XCORRwtheta_clustering_noW(np.hstack((evCLU['RA'], evAGN['RA'])), np.hstack((evCLU['DEC'],evAGN['DEC'])),
											GAL['RA'][s10], GAL['DEC'][s10],
											ra_fs[s_R], dec_fs[s_R], p_2_2PCF)
		basename = 'GAL_m'+str(np.round(m0,1))+'_evAGNevCLUevBKG'
		p_2_2PCF = os.path.join(dir_2pcf, basename +'_CROSSCORR_05E20.wtheta.2pcf.fits' )
		if not os.path.isfile(p_2_2PCF):
			tabulate_XCORRwtheta_clustering_noW(np.hstack((evCLU['RA'], evAGN['RA'], evBKG['RA'])), np.hstack((evCLU['DEC'], evAGN['DEC'], evBKG['DEC'])),
											GAL['RA'][s10], GAL['DEC'][s10],
											ra_fs[s_R], dec_fs[s_R], p_2_2PCF)




