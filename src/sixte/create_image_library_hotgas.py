
import time
t0 = time.time()
import astropy.io.fits as fits
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import os, glob
from scipy.interpolate import interp1d
from scipy.interpolate import RegularGridInterpolator
from scipy.stats import norm
from astropy.table import Table
import numpy as np


def write_img(matrix, out='spherical_cc', n_pixel = 120, angularSize_per_pixel = 0.033):
    prihdr = fits.Header()
    prihdr['HDUCLASS'] = 'HEASARC/SIMPUT'
    prihdr['HDUCLAS1'] = 'IMAGE'
    prihdr['HDUVERS'] = '1.1.0'
    prihdr['EXTNAME'] = 'IMAGE'
    prihdr['CTYPE1'] = ('RA---TAN', 'first axis (column) is Right Ascension')
    prihdr['CRPIX1'] = ((n_pixel+1) / 2., 'middle pixel of array in col direction')
    prihdr['CRVAL1'] = (0, 'Dec of this middle pixel, in degrees')
    prihdr['CDELT1'] = (-angularSize_per_pixel/60., 'move 1column forward,decrease RA by CDELT1/deg')
    prihdr['CROTA1'] = 0
    prihdr['CUNIT1'] = 'deg'
    prihdr['CTYPE2'] = ('DEC--TAN', 'first axis (column) is Declination')
    prihdr['CRPIX2'] = ((n_pixel+1) / 2., 'middle pixel of array in row direction')
    prihdr['CRVAL2'] = (0, 'RA of this middle pixel, in degrees')
    prihdr['CDELT2'] = (angularSize_per_pixel/60., 'move 1column forward,increase Dec by CDELT1/deg')
    prihdr['CROTA2'] = 0
    prihdr['CUNIT2'] = 'deg'
    prihdr['EQUINOX'] = 2000
    prihdr['RADECSYS'] = 'FK5'
    prihdr['aMinpPix'] = angularSize_per_pixel
    prihdu = fits.PrimaryHDU(matrix, header=prihdr)
    #if os.path.isfile(out):
        #os.remove(out)
    prihdu.writeto(out, overwrite=True)


# b_a=0.71
def create_matrix(profile, n_pixel = 121, b_a = 0.7, truncation_radius = 12.):
	x_max = truncation_radius
	#sel = (profile.y < profile.y.max()/20)
	#x_max = np.min([np.min(profile.x[sel]), profile.x[-2]])
	angularSize_per_pixel = x_max/(n_pixel/2.)
	matrix = np.zeros((n_pixel, n_pixel))
	xxx = (np.arange(n_pixel) - (n_pixel-1) / 2.) * angularSize_per_pixel
	x_matrix, y_matrix = np.meshgrid(xxx, xxx)
	r_matrix = ((x_matrix / b_a)**2 + y_matrix**2)**0.5
	matrix = profile(r_matrix)
	matrix = matrix / np.sum(matrix)
	return matrix, angularSize_per_pixel


cosmoUCHUU = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
h = 0.6774
cosmo = cosmoUCHUU
zs = np.arange(0.0000001, 7.1, 0.001)
dm_itp = interp1d(zs, cosmo.distmod(zs).value)

from colossus.cosmology import cosmology
cosmology.setCosmology('planck18')
from colossus.halo import mass_so
from colossus.halo import mass_defs
from colossus.halo import concentration

xgrid_ext = np.loadtxt(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/model_GAS', 'radial_binning.txt'))

p_2_profiles = os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/model_GAS', 'profiles_010z015_1e14M2e14.fits')
#
snapshot, redshift, scale_factor = np.loadtxt( os.path.join(os.environ['UCHUU'], 'snap_list.txt'), unpack = True)
all_zs = redshift[(redshift>0)&(redshift<1.5)]

pixel_size_image = 20 # arc second

# loop over redshift
for j_z in np.arange(len(all_zs)):
    print('='*100)
    z_bar = all_zs[j_z]
    z_str = 'z'+str(np.round(z_bar,5))
    img_dir = os.path.join(os.environ['UCHUU'], 'cluster_images', z_str )
    os.system('mkdir -p '+img_dir)
    print(img_dir)

    # loop over ellipticity
    j_e = 0
    b_to_a_500c = [0.75]
    b_a = b_to_a_500c[j_e]
    e_str = str( np.round( b_a, 2) )
    conversion_arcmin = cosmo.kpc_proper_per_arcmin(z_bar).value
    print('creates SIMPUT images',z_bar,conversion_arcmin,b_a)
    path_2_images = []
    t_prof = Table.read( p_2_profiles )
    p_2_profiles_images = os.path.join(img_dir, 'profiles_010z015_1e14M2e14_p2images.fits')
    angularSize_per_pixel = np.zeros(len(t_prof))

    # loop over profile
    for j_p in np.arange(len(t_prof)):
        prf = t_prof[j_p]
        file_name = 'profileLineID_'+str(j_p).zfill(5)+'_ba_'+e_str+'_'+z_str
        image_file = os.path.join(img_dir, file_name+'.fits')
        #print(image_file)
        path_2_images.append( image_file )
        r500c_i = prf['R500c']
        profile_i = prf['profiles']/prf['profiles'].max()
        x_coord = np.hstack(( 0., xgrid_ext*r500c_i/conversion_arcmin, 10*xgrid_ext[-1]*r500c_i/conversion_arcmin, 10000000 ))
        y_coord = np.hstack(( profile_i[0], profile_i, 0., 0. ))
        profile = interp1d(x_coord, y_coord)
        truncation_radius = 2 * r500c_i/ conversion_arcmin
        n_pixel = 2*int(truncation_radius*60/pixel_size_image)+1
        matrix, angularSize_per_pixel_j = create_matrix(profile, n_pixel = n_pixel, b_a = b_a, truncation_radius = truncation_radius)
        angularSize_per_pixel[j_p] = angularSize_per_pixel_j
        #print(r500c_i, conversion_arcmin, r500c_i/ conversion_arcmin, 60*angularSize_per_pixel_j, n_pixel, image_file)
        #if os.path.isfile(image_file)==False:
        write_img(matrix, image_file, n_pixel = n_pixel, angularSize_per_pixel=angularSize_per_pixel_j)
        print(image_file, 'written')

    t_prof['path_2_images'] = path_2_images
    t_prof['angularSize_per_pixel'] = angularSize_per_pixel
    t_prof.write(p_2_profiles_images, overwrite = True)
    print(p_2_profiles_images, 'written')
