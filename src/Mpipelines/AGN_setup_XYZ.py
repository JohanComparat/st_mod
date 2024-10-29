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

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoUNIT = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
h = 0.6774
L_box = 1000.0 / h
cosmo = cosmoUNIT

z_array = np.arange(0, 7.5, 0.001)
dc_to_z = interp1d(cosmo.comoving_distance(z_array), z_array)
z_to_dc = interp1d(z_array, cosmo.comoving_distance(z_array))

d_L = cosmo.luminosity_distance(z_array)
dl_cm = (d_L.to(u.cm)).value
dL_interpolation = interp1d(z_array, dl_cm)

path_2_agnOut_file = os.path.join( os.environ["UCHUU"], 'AGN_LX_tables')

# computes the cosmological volume
DZ = 0.01
z_bins = np.arange(0.0, 6., DZ)
##

def get_xyz(RARA, DECDEC, ZZZZ):
	phi   = ( RARA   - 180 ) * np.pi / 180.
	theta = ( DECDEC + 90 ) * np.pi / 180.
	rr    = z_to_dc( ZZZZ )
	xx = rr * np.cos( phi   ) * np.sin( theta )
	yy = rr * np.sin( phi   ) * np.sin( theta )
	zz = rr * np.cos( theta )
	return xx, yy, zz


def tabulate_AGN(z_bins_i = z_bins[0]):
	print("="*100)
	print(z_bins_i)
	t0 = time.time()
	p_2_AGN = os.path.join( path_2_agnOut_file , 'LX_table_' +str(np.round(z_bins_i,2))+ '.fits' )
	p_2_out = os.path.join( path_2_agnOut_file , 'XYZ_position_' +str(np.round(z_bins_i,2))+ '.fits' )
	print('starts ' , p_2_out, time.time() - t0)
	t_in = Table.read(p_2_AGN)
	t_in['z']
	# create random RA, DEC,
	t_out = Table()
	t_out['z'] = t_in['z']
	t_out['dC'] = z_to_dc(t_out['z'])
	size = len(t_in['z'])
	print(size)
	uu = np.random.uniform(size=size)
	t_out['DEC'] = np.arccos(1 - 2 * uu) * 180 / np.pi - 90.
	t_out['RA'] = np.random.uniform(size=size) * 2 * np.pi * 180 / np.pi
	# deduce X,Y,Z
	X, Y, Z = get_xyz(t_out['RA'], t_out['DEC'], t_out['z'])
	t_out['X'] = X
	t_out['Y'] = Y
	t_out['Z'] = Z
	print(t_out)
	t_out.write( p_2_out, overwrite = True )
	print(p_2_out, 'written in ', time.time() - t0, 's')

tabulate_AGN(z_bins_i = z_bins[jj0])
