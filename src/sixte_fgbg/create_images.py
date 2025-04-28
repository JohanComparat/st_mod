"""
On he1srv 


cd software/linux
sxrbg-data.tar.gz
sxrbg.tar.gz

tar -xzf sxrbg-data.tar.gz
tar -xzf sxrbg.tar.gz

cd sxrbg 

hmake 

cp sxrbg.par /home/comparat/pfiles/

tail sxrbg.hlp
sxrbg lii=115.26 bii=-30.08 disio=1.0 disin=0.5 searchmode=annulus
sxrbg lii=115.26 bii=-30.08 disio=1.0 searchmode=cone


count rate x 1.38e-11 => FLUX erg/cm2/s


"""
import numpy as n
from astropy_healpix import healpy
from astropy.table import Table, Column
import sys, os, time
import astropy.units as u
from astropy.coordinates import SkyCoord
import subprocess
import astropy.io.fits as fits
from astropy import wcs

path_2_SMPT_images = '/data40s/erosim/eRASS/bkg_erosita_simput_full_sky/images'

NSIDE = 128 
N_pixels_all = healpy.nside2npix(NSIDE)
area_per_pixel =  healpy.nside2pixarea(NSIDE, degrees=True)
radius = (area_per_pixel/n.pi)**0.5

pix_ids = n.arange(N_pixels_all)
ra_cen, dec_cen = healpy.pix2ang(NSIDE, pix_ids, lonlat=True, nest=True)
#healpy.ang2pix(NSIDE, ra_cen, dec_cen, lonlat=True, nest=True) 

coords = SkyCoord(ra_cen, dec_cen, unit='deg', frame='icrs')
g_lat = coords.galactic.b.value
g_lon = coords.galactic.l.value

def write_image_000(pixel_id = 0, out = 'test.fits', NSIDE = 128 ):
	"""
	Sets the WCS and the matrix is filled with zeros
	"""
	print(pixel_id)
	area_per_pixel =  healpy.nside2pixarea(NSIDE, degrees=True)
	radius = (area_per_pixel/n.pi)**0.5
	CDELT = 60 # arcseconds
	n_pixel = int( int(1.5*radius*3600/CDELT) * 2 ) + 1
	matrix = n.zeros((n_pixel, n_pixel)).astype('int')
	xxx = (n.arange(n_pixel) - (n_pixel-1) / 2.) 
	x_matrix, y_matrix = n.meshgrid(xxx, xxx)
	ra_mat, dec_mat = healpy.pix2ang(NSIDE, pixel_id, lonlat=True, nest=True)
	matrix = n.zeros((n_pixel, n_pixel))
	prihdr = fits.Header()
	prihdr['HDUCLASS'] = 'HEASARC/SIMPUT'
	prihdr['HDUCLAS1'] = 'IMAGE'
	prihdr['HDUVERS'] = '1.1.0'
	prihdr['EXTNAME'] = 'IMAGE'
	prihdr['CTYPE1'] = ('RA---TAN', 'first axis (column) is Right Ascension')
	prihdr['CRPIX1'] = ((n_pixel+1) / 2., 'middle pixel of array in col direction')
	prihdr['CRVAL1'] = (ra_mat, 'RA of this middle pixel, in degrees')
	prihdr['CDELT1'] = (-CDELT/3600., 'move 1column forward,decrease RA by CDELT1/deg')
	prihdr['CROTA1'] = 0
	prihdr['CUNIT1'] = 'deg'
	prihdr['CTYPE2'] = ('DEC--TAN', 'first axis (column) is Declination')
	prihdr['CRPIX2'] = ((n_pixel+1) / 2., 'middle pixel of array in row direction')
	prihdr['CRVAL2'] = (dec_mat, 'Dec of this middle pixel, in degrees')
	prihdr['CDELT2'] = (CDELT/3600., 'move 1column forward,increase Dec by CDELT1/deg')
	prihdr['CROTA2'] = 0
	prihdr['CUNIT2'] = 'deg'
	prihdr['EQUINOX'] = 2000
	prihdr['RADECSYS'] = 'FK5'
	prihdu = fits.PrimaryHDU(matrix, header=prihdr)
	if os.path.isfile(out):
		os.remove(out)
	prihdu.writeto(out, overwrite=True)

def write_matrix(pixel_id = 0, out = 'test.fits', NSIDE = 128 ):
	"""
	computes the values of the matrix
	"""
	print(pixel_id)
	hdu = fits.open( out )
	w = wcs.WCS(hdu[0].header)
	matrix = hdu[0].data # n.zeros((n_pixel, n_pixel)).astype('int')
	n_pixel = len(matrix[0])
	xxx = n.arange(n_pixel)  
	x_matrix, y_matrix = n.meshgrid(xxx, xxx)
	for ii, (x_element, y_element) in enumerate( zip(x_matrix, y_matrix) ):
		world_ra, world_dec = w.all_pix2world(n.transpose([x_element, y_element]), 0, ra_dec_order=True).T
		pixel_values = healpy.ang2pix(NSIDE, world_ra, world_dec, lonlat=True, nest=True) 
		#print(set(pixel_values))
		matrix_values = n.zeros_like(pixel_values)
		matrix_values[pixel_values==pixel_id] = 1
		matrix[ii] = matrix_values

	prihdu = fits.PrimaryHDU(matrix, header=hdu[0].header)
	if os.path.isfile(out):
		os.remove(out)
	prihdu.writeto(out, overwrite=True)

print(NSIDE)
print("writes empty images")
out8 = n.array([write_image_000(iii, os.path.join(path_2_SMPT_images, 'image_'+str(iii).zfill(6)+'.fits'), NSIDE = NSIDE ) for iii in n.arange(healpy.nside2npix(NSIDE)) ])
print("writes image values")
n.array([write_matrix(iii, os.path.join(path_2_SMPT_images, 'image_'+str(iii).zfill(6)+'.fits'), NSIDE = NSIDE ) for iii in n.arange(healpy.nside2npix(NSIDE)) ])
