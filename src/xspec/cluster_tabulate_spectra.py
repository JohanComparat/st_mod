"""
What it does
------------
Tabulates APEC X-ray spectra

Command to run
--------------
# ds43
# conda activate eroconda
python cluster_tabulate_spectra.py

"""
from astropy_healpix import healpy
import sys
import os
import time
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import astropy.constants as cc
import astropy.io.fits as fits
from scipy.special import erf
from scipy.stats import norm
from scipy.interpolate import interp2d
from scipy.interpolate import interp1d
from astropy.table import Table, Column
from scipy.linalg import eigh, cholesky
from scipy.stats import norm

import h5py
import numpy as n
print('Tabulate K-correction')
print('------------------------------------------------')
print('------------------------------------------------')
t0 = time.time()

# import all pathes

# initializes pathes to files
dir_2_result = os.path.join(os.environ['GIT_STMOD'], 'data/models/model_GAS', 'xray_k_correction')
os.system('mkdir -p '+dir_2_result)

cosmoUCHUU = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
h = 0.6774
cosmo = cosmoUCHUU

import xspec
xspec.Xset.cosmo = "67.74 0. 0.6911"

nh_vals = 10**n.arange(-2,4+0.01,0.5)#0.05)
kT_vals = 10**n.arange(-2.09,1.8,0.01)
z_vals = n.hstack((n.array([0.]), n.arange(0.001, 0.01, 0.001), 10**n.arange(n.log10(0.01), n.log10(6.1), 0.01)))
#z_vals = n.hstack(( n.arange(0.,0.7,0.05), n.arange(0.8, 4.5, 0.1) ))
#kT_vals = n.hstack(([0.1, 0.2], n.arange(0.5,8,0.5), [10, 20, 30, 40, 50] ))

def get_fraction_obsF_RF(nh_val=1, kT=4, redshift=0, norm=1, metallicity=0.3):
	"""
	computes the fraction of flux observed in the eROSITA band w.r.t. RF 0.5-2 keV flux using a tbabs*apec model (cluster)
	model at a given redshift
	returns :
	  * flux(0.5-2,nH) / flux(0.5/(1+z)-2/(1+z),nH=0.01)
	nh_val=1 
	kT=4 
	redshift=0 
	norm=1
	"""
	print(nh_val, redshift)
	kev_min_erosita = 0.5 # 24.8
	kev_max_erosita = 2.0 # 6.19
	kev_min_erosita_RF = 0.5*(1+redshift)
	kev_max_erosita_RF = 2.0*(1+redshift)
	m1 = xspec.Model("tbabs*apec")
	m1.setPars(
		nh_val,      #   1    1   TBabs      nH         10^22    1.00000      +/-  0.0          
		kT,          #   2    2   apec       kT         keV      1.00000      +/-  0.0          
		metallicity, #   3    2   apec       Abundanc            1.00000      frozen  
		0.,          #   4    2   apec       Redshift            0.0          frozen
		norm,        #   5    2   apec       norm                1.00000      +/-  0.0          
		)
	# flux observed
	xspec.AllModels.calcFlux(str(kev_min_erosita)+" "+str(kev_max_erosita))
	flux_obs = m1.flux[0]#/(kev_max_erosita-kev_min_erosita)
	xspec.AllModels.show()
	# flux in the rest-frame without galactic absorption
	m1.TBabs.nH = 0.0
	xspec.AllModels.calcFlux(str(kev_min_erosita_RF)+" "+str(kev_max_erosita_RF))
	flux_intrinsic = m1.flux[0]#/(kev_max_erosita_RF-kev_min_erosita_RF)
	xspec.AllModels.show()
	fraction_observed = flux_obs / flux_intrinsic
	return fraction_observed


def get_fraction_soft_broadSoft(nh_val=1, kT=4, redshift=0, norm=1, metallicity=0.3):
	"""
	computes the fraction of flux observed in the eROSITA band w.r.t. RF 0.5-2 keV flux using a tbabs*apec model (cluster)
	model at a given redshift
	returns :
	  * flux(0.5-2,nH) / flux(0.5/(1+z)-2/(1+z),nH=0.01)
	nh_val=1 
	kT=4 
	redshift=0 
	norm=1
	"""
	print(nh_val, redshift)
	kev_min_erosita = 0.5*(1+redshift) # 24.8
	kev_max_erosita = 2.0*(1+redshift) # 6.19
	kev_min_erosita_RF = 0.1*(1+redshift)
	kev_max_erosita_RF = 2.4*(1+redshift)
	m1 = xspec.Model("tbabs*apec")
	m1.setPars(
		nh_val,      #   1    1   TBabs      nH         10^22    1.00000      +/-  0.0          
		kT,          #   2    2   apec       kT         keV      1.00000      +/-  0.0          
		metallicity, #   3    2   apec       Abundanc            1.00000      frozen  
		0.,          #   4    2   apec       Redshift            0.0          frozen
		norm,        #   5    2   apec       norm                1.00000      +/-  0.0          
		)
	# flux observed
	xspec.AllModels.calcFlux(str(kev_min_erosita)+" "+str(kev_max_erosita))
	flux_obs = m1.flux[0]#/(kev_max_erosita-kev_min_erosita)
	xspec.AllModels.show()
	# flux in the rest-frame without galactic absorption
	m1.TBabs.nH = 0.0
	xspec.AllModels.calcFlux(str(kev_min_erosita_RF)+" "+str(kev_max_erosita_RF))
	flux_intrinsic = m1.flux[0]#/(kev_max_erosita_RF-kev_min_erosita_RF)
	xspec.AllModels.show()
	fraction_observed = flux_obs / flux_intrinsic
	return fraction_observed


def get_attenuation(nh_val=1, kT = 2):
	"""
	computes the attenuation
	"""
	print(nh_val)
	kev_min_erosita = 0.5 # 24.8
	kev_max_erosita = 2.0 # 6.19
	kev_min_erosita_RF = 0.5
	kev_max_erosita_RF = 2.0
	m1 = xspec.Model("tbabs*apec")
	m1.setPars(nh_val, kT, 0.3, 0., 1.,)
	# flux observed
	xspec.AllModels.calcFlux(str(kev_min_erosita)+" "+str(kev_max_erosita))
	flux_obs = m1.flux[0]#/(kev_max_erosita-kev_min_erosita)
	xspec.AllModels.show()
	# flux in the rest-frame without galactic absorption
	m1.TBabs.nH = 0.0
	xspec.AllModels.calcFlux(str(kev_min_erosita_RF)+" "+str(kev_max_erosita_RF))
	flux_intrinsic = m1.flux[0]#/(kev_max_erosita_RF-kev_min_erosita_RF)
	xspec.AllModels.show()
	fraction_observed = flux_obs / flux_intrinsic
	return fraction_observed


if os.path.isfile(os.path.join( dir_2_result, "nh_attenuation_clusters_kt1p0.txt"))==False:
	nh_vals = 10**n.arange(-3,4+0.01,0.01)
	fraction_observed = n.zeros_like(nh_vals)

	for jj, nh_val in enumerate(nh_vals):
		fraction_observed[jj] = get_attenuation(nh_val, kT = 1)

	n.savetxt( os.path.join( dir_2_result, "nh_attenuation_clusters_kt1p0.txt"), n.transpose([22 + n.log10(nh_vals), fraction_observed]), header = 'log_nh fraction_observed' )


if os.path.isfile(os.path.join( dir_2_result, "nh_attenuation_clusters_kt2p0.txt"))==False:
	nh_vals = 10**n.arange(-3,4+0.01,0.01)
	fraction_observed = n.zeros_like(nh_vals)

	for jj, nh_val in enumerate(nh_vals):
		fraction_observed[jj] = get_attenuation(nh_val, kT = 2)

	n.savetxt( os.path.join( dir_2_result, "nh_attenuation_clusters_kt2p0.txt"), n.transpose([22 + n.log10(nh_vals), fraction_observed]), header = 'log_nh fraction_observed' )


if os.path.isfile(os.path.join( dir_2_result, "nh_attenuation_clusters_kt0p5.txt"))==False:
	nh_vals = 10**n.arange(-3,4+0.01,0.01)
	fraction_observed = n.zeros_like(nh_vals)

	for jj, nh_val in enumerate(nh_vals):
		fraction_observed[jj] = get_attenuation(nh_val, kT = 0.5)

	n.savetxt( os.path.join( dir_2_result, "nh_attenuation_clusters_kt0p5.txt"), n.transpose([22 + n.log10(nh_vals), fraction_observed]), header = 'log_nh fraction_observed' )


if os.path.isfile(os.path.join( dir_2_result, "fraction_observed_clusters.txt"))==False:
	# rest-frame 0.5-2 to observed frame 0.5-2
	XX_nh, YY_z, ZZ_kt = n.meshgrid(nh_vals, z_vals, kT_vals)

	XX_nh = n.hstack((n.hstack((XX_nh))  ))
	YY_z  = n.hstack((n.hstack((YY_z))   ))
	ZZ_kt = n.hstack((n.hstack((ZZ_kt))  ))
	fraction_observed = n.zeros_like(XX_nh)

	for jj, (nh_val, redshift_val, kt_val) in enumerate(zip(XX_nh, YY_z, ZZ_kt)):
		try :
			fraction_observed[jj] = get_fraction_obsF_RF(nh_val, kt_val, redshift_val)
		except(ZeroDivisionError):
			print('no data')

	n.savetxt( os.path.join( dir_2_result, "fraction_observed_clusters.txt"), n.transpose([22 + n.log10(XX_nh), YY_z, ZZ_kt, fraction_observed]), header = 'log_nh z kT fraction_observed' )

sys.exit()

if os.path.isfile(os.path.join( dir_2_result, "fraction_observed_clusters_no_nH.txt"))==False:
	# no NH

	# rest-frame 0.5-2 to observed frame 0.5-2
	YY_z, ZZ_kt = n.meshgrid(z_vals, kT_vals)

	YY_z  =n.hstack((YY_z))
	ZZ_kt =n.hstack((ZZ_kt))
	fraction_observed = n.zeros_like(YY_z)

	for jj, (redshift_val, kt_val) in enumerate(zip(YY_z, ZZ_kt)):
		try :
			fraction_observed[jj] = get_fraction_obsF_RF(0., kt_val, redshift_val)
		except(ZeroDivisionError):
			print('no data')

	n.savetxt( os.path.join( dir_2_result, "fraction_observed_clusters_no_nH.txt"), n.transpose([YY_z, ZZ_kt, fraction_observed]), header = 'z kT fraction_observed' )


"""
for jj, (nh_val, redshift_val, kt_val) in enumerate(zip(XX_nh, YY_z, ZZ_kt)):
	fraction_observed[jj] = get_fraction_soft_broadSoft(nh_val, kt_val, redshift_val) 

n.savetxt( os.path.join( dir_2_result, "fraction_05_20_01_24_no_nH.txt"), n.transpose([YY_z, ZZ_kt, fraction_observed]), header = 'z kT fraction_observed' )
"""
