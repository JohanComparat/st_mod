"""
Tabulates the AGN model from Comparat et al. 2019.
"""
import sys, os, glob, time
t0 = time.time()

from astropy.table import Table
from scipy.special import erf
##from scipy.interpolate import interp2d
from scipy.interpolate import LinearNDInterpolator # SmoothBivariateSpline

from scipy.interpolate import interp1d
import numpy as np
print('Creates AGN mock catalogue ')
print('------------------------------------------------')
print('------------------------------------------------')

# minimum luminosity kept (rejects sources less bright)
env = "UCHUU"
jj0 = int(sys.argv[1])

# link to X-ray K-correction and attenuation curves
path_2_hard_RF_obs_soft = os.path.join(
    os.environ['GIT_STMOD'],
    "data", "models", "model_AGN",
    "xray_k_correction",
    "v3_fraction_observed_A15_RF_hard_Obs_soft_fscat_002.txt")

path_2_RF_obs_hard = os.path.join(
    os.environ['GIT_STMOD'],
    "data", "models", "model_AGN",
    "xray_k_correction",
    "v3_fraction_observed_A15_RF_hard_Obs_hard_fscat_002.txt")

path_2_obs_hard_obs_soft = os.path.join(
    os.environ['GIT_STMOD'],
    "data", "models", "model_AGN",
    "xray_k_correction",
    "v3_fraction_observed_A15_Obs_hard_Obs_soft_fscat_002.txt")

path_2_hard_RF_obs_soft_3D = np.array( glob.glob( os.path.join(
    os.environ['GIT_STMOD'],
    "data", "models", "model_AGN",
    "xray_k_correction",
    "v3_fraction_observed_A15_RF_hard_Obs_soft_fscat_002_GALnH_*.txt") ) )
path_2_hard_RF_obs_soft_3D.sort()

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoUNIT = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
h = 0.6774
L_box = 1000.0 / h
cosmo = cosmoUNIT

z_array = np.arange(0, 7.5, 0.001)
dc_to_z = interp1d(cosmo.comoving_distance(z_array), z_array)

d_L = cosmo.luminosity_distance(z_array)
dl_cm = (d_L.to(u.cm)).value
dL_interpolation = interp1d(z_array, dl_cm)

path_2_agnOut_file = os.path.join( os.environ["UCHUU"], 'AGN_LX_tables')

# computes the cosmological volume
area = 129600/np.pi    # deg2
DZ = 0.01
z_bins = np.arange(0.0, 6., DZ)
##
log10_FX_lim = -18
LX_min_1 = np.log10( 4 * np.pi * dL_interpolation(z_bins[jj0])**2 * 10**log10_FX_lim )
LX_min_0 = 38.
LX_min = np.max([LX_min_0, LX_min_1])
print( 'z min = ', z_bins[jj0] )
print( LX_min, '=max(', LX_min_0, LX_min_1, ')' )

def tabulate_AGN(z_bins_i = z_bins[0]):
	print("="*100)
	print(z_bins_i)
	t0 = time.time()
	p_2_out = os.path.join( path_2_agnOut_file , 'LX_table_' +str(np.round(z_bins_i,2))+ '.fits' )
	print('starts ' , p_2_out, time.time() - t0)
	t_out = Table()
	zmin = z_bins_i
	zmax = z_bins_i + DZ
	z_mean = 0.5 * (zmin + zmax)
	if zmin==0:
		vol = (cosmo.comoving_volume(zmax).value) * area * np.pi / 129600.
	else:
		vol = (cosmo.comoving_volume(zmax).value - cosmo.comoving_volume(zmin).value) * area * np.pi / 129600.
	DL_mean_z = (cosmo.luminosity_distance(z_mean).to(u.cm)).value
	print(zmin, '<z<', zmax, ',volume', vol/1e9, 'Gpc3')
	# Hard LX Abundance Matching
	# Equations 2 and 3 of Comparat et al. 2019
	def kz_h(z): return 10**(-4.03 - 0.19 * (1 + z))

	def Ls_h(z): return 10**(44.84 - np.log10(((1 + 2.0) / (1 + z))** 3.87 + ((1 + 2.0) / (1 + z))**(-2.12)))

	def phi_h(L, z): return kz_h(z) / ((L / Ls_h(z))**0.48 + (L / Ls_h(z))**2.27)

	# LF in the mock, starting parameters
	dlogf = 0.05
	Lbin_min = 36
	fbins = np.arange(Lbin_min, 48, dlogf)
	xf = fbins[:-1] + dlogf / 2.

	# theoretical number counts and LF
	N_obs_th = phi_h(10**xf, z_mean * np.ones_like(xf)) * vol * dlogf

	t1 = time.time()
	# select bins with a number of AGN greater than 1 and smaller than 2x the
	# total number of agn, we want to simulate
	n_agn = 20_000_000
	bin_selection = (N_obs_th >= 0.5) & (N_obs_th < n_agn * 2.)
	print(np.transpose([xf[bin_selection], N_obs_th[bin_selection].astype('int') + 1]))
	# draw LX luminosities uniformly in each LX bin, the bins (dlogf = 0.05)
	# are small enough for a uniform sampling
	all_x = []
	for aa, bb, cc in zip(fbins[:-1][bin_selection], fbins[1:][bin_selection], N_obs_th[bin_selection].astype('int') + 1):
		all_x.append(np.random.uniform(low=aa, high=bb, size=int(cc)))

	X_luminosities = np.hstack((all_x))
	X_luminosities_sorted = X_luminosities[np.argsort(X_luminosities)]

	# select with the flux limit the right number of AGN :
	selected = X_luminosities_sorted > LX_min
	n_agn = len(X_luminosities_sorted[selected])
	print('N AGN=', n_agn, n_agn<10_000_000, 'dt=', time.time() - t0)
	z = np.random.uniform(zmin, zmax, n_agn)
	dl_cm = dL_interpolation(z)
	# output table
	#t_out['z'] = z
	t_out['LX_hard'] = X_luminosities_sorted[selected]
	lx = t_out['LX_hard']
	print('hard LX computed, dt=', time.time() - t0)

	#=============================
	# Obscured fractions
	# ===============================
	# model from equations 4-11, 12-15 of Comparat et al. 2019

	# too many CTK at high luminosity
	# Eq. 4
	#def f_thick(LXhard, z): return 0.30
	def thick_LL(z, lx0 = 41.5): return lx0 + np.arctan(z*5)*1.5
	def f_thick(LXhard, z): return 0.30 * (0.5 + 0.5 * erf((thick_LL(z) - LXhard) / 0.25))

	# too many absorbed ones
	# Eq. 7
	def f_2(LXhard, z): return 0.9 * (41 / LXhard)**0.5

	# fiducial
	# Eq. 8
	def f_1(LXhard, z): return f_thick(LXhard, z) + 0.01 + erf(z / 4.) * 0.3

	# Eq. 10
	def LL(z, lx0 = 43.2): return lx0 + erf(z) * 1.2

	# Eq. 5,6
	def fraction_ricci(LXhard, z, width = 0.6): return f_1(LXhard,z) + (f_2(LXhard, z) - f_1(LXhard,z)) * (0.5 + 0.5 * erf((LL(z) - LXhard) / width))

	# initializes logNH
	logNH = np.zeros(n_agn)

	# obscuration, after the equations above
	randomNH = np.random.rand(n_agn)

	# unobscured 20-22
	#frac_thin = fraction_ricci(lsar, z)
	frac_thin = fraction_ricci(lx, z)
	thinest = (randomNH >= frac_thin)

	# thick obscuration, 24-26
	thick = (randomNH < f_thick(lx, z))
	#thick = (randomNH < thick_fraction)

	# obscured 22-24
	obscured = (thinest == False) & (thick == False)

	# assigns logNH values randomly :
	logNH[thick] = np.random.uniform(24, 26, len(logNH[thick]))
	logNH[obscured] = np.random.uniform(22, 24, len(logNH[obscured]))
	logNH[thinest] = np.random.uniform(20, 22, len(logNH[thinest]))
	print('logNH computed, dt=', time.time() - t0)

	#print('=====================  AGN fractions and numbers vs NH values =================')
	#print(n_agn,
		#len(thick.nonzero()[0]) * 1. / n_agn,
		#len(obscured.nonzero()[0]) * 1. / n_agn,
		#len(thinest.nonzero()[0]) * 1. / n_agn)

	# ===============================
	# Assigns flux
	# ===============================

	NHS = np.arange(20, 26 + 0.05, 0.4)
	# hard X-ray 2-10 keV rest-frame ==>> 2-10 obs frame
	obscuration_z_grid, obscuration_nh_grid, obscuration_fraction_obs_erosita = np.loadtxt(
		path_2_RF_obs_hard, unpack=True)

	xy = np.c_[obscuration_z_grid, obscuration_nh_grid]
	obscuration_itp_H_H = LinearNDInterpolator(xy, obscuration_fraction_obs_erosita)
	##obscuration_itp_H_H_itp2d = interp2d(
		##obscuration_z_grid,
		##obscuration_nh_grid,
		##obscuration_fraction_obs_erosita)

	percent_observed_itp = interp1d(
		np.hstack((20 - 0.1, NHS, 26 + 0.1)),
		np.hstack((
			obscuration_itp_H_H(z_mean, 20.),
			np.array([obscuration_itp_H_H(z_i, logNH_i) for z_i, logNH_i in zip(z_mean * np.ones_like(NHS), NHS)]),
			obscuration_itp_H_H(z_mean, 26.))))
	percent_observed_H_H = percent_observed_itp(logNH)

	lx_obs_frame_2_10 = np.log10(10**lx * percent_observed_H_H)
	fx_2_10 = 10**(lx_obs_frame_2_10) / (4 * np.pi * (dl_cm)**2.) # / h**3
	print('flux 2-10 computed, dt=', time.time() - t0)


	def fraction_22p21_merloni(lx): return (
		0.5 + 0.5 * erf((-lx + 44.) / 0.9)) * 0.69 + 0.26


	def compute_agn_type(z, lx, logNH, fbins=fbins, n_agn=n_agn):
		"""
		Assigns a type to an AGN population

		parameters:
		- z: redshift
		- lx: hard X-ray luminosity (log10)
		- logNH: nH value (log10)

		return: array of AGN types
		"""
		# boundary between the 22 and the 21 populations
		limit = fraction_22p21_merloni((fbins[1:] + fbins[:-1]) * 0.5)
		# selection per obscuration intensity
		nh_21 = (logNH <= 22.)
		nh_23 = (logNH > 22.)  # &(logNH<=26.)
		# initiate columns to compute
		opt_type = np.zeros(n_agn).astype('int')
		rd = np.random.rand(n_agn)
		# compute histograms of LX for different obscurations
		nall = np.histogram(lx, fbins)[0]       # all
		nth = np.histogram(lx[nh_23], fbins)[0]  # thin
		nun = np.histogram(lx[nh_21], fbins)[0]  # unobscured
		fr_thk = nth * 1. / nall  # fraction of obscured
		fr_un = nun * 1. / nall  # fraction of unobscured
		# first get the type 12: NH absorption but optically unobscured
		# to be chosen in obscured population
		n_per_bin_12 = (fr_thk - limit) * nall
		sel_12 = (np.ones(len(z)) == 0)
		for bin_low, bin_high, num_needed, nn_un in zip(
				fbins[:-1], fbins[1:], n_per_bin_12.astype('int'), nth):
			if num_needed > 0 and nn_un > 0:
				frac_needed = num_needed * 1. / nn_un
				sel_12 = (sel_12) | (
					(lx > bin_low) & (
						lx < bin_high) & (nh_23) & (
						rd < frac_needed))
		t_12 = (nh_23) & (sel_12)
		# second the types 21
		# to be chosen in nun
		n_per_bin_21 = (-fr_thk + limit) * nall
		sel_21 = (np.ones(len(z)) == 0)
		for bin_low, bin_high, num_needed, nn_un in zip(
				fbins[:-1], fbins[1:], n_per_bin_21.astype('int'), nun):
			if num_needed > 0 and nn_un > 0:
				frac_needed = num_needed * 1. / nn_un
				sel_21 = (sel_21) | (
					(lx > bin_low) & (
						lx < bin_high) & (nh_21) & (
						rd < frac_needed))
		t_21 = (nh_21) & (sel_21)
		# finally the types 11 and 22
		t_11 = (nh_21) & (t_21 == False)
		t_22 = (nh_23) & (t_12 == False)
		opt_type[t_22] = 22
		opt_type[t_12] = 12
		opt_type[t_11] = 11
		opt_type[t_21] = 21
		return opt_type


	opt_type = compute_agn_type(z, lx, logNH)
	print('type added, dt=', time.time() - t0)

	# ===============================
	# Writing results
	# ===============================

	t_out['FX_hard'] = fx_2_10
	t_out['logNH'] = logNH
	t_out['agn_type'] = opt_type
	t_out.write( p_2_out, overwrite = True )
	print(p_2_out, 'written in ', time.time() - t0, 's')

tabulate_AGN(z_bins_i = z_bins[jj0])
