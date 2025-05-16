import time
t0 = time.time()
import astropy.io.fits as fits
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import os, glob, sys
from scipy.interpolate import interp1d
from scipy.interpolate import RegularGridInterpolator
from scipy.stats import norm
from astropy.table import Table, vstack
import numpy as np
import astropy.io.fits as fits


def write_img(matrix, out='spherical_cc', n_pixel = 120, angularSize_per_pixel = 0.033, RA_CEN=0.0, DEC_CEN=0.0):
    prihdr = fits.Header()
    prihdr['HDUCLASS'] = 'HEASARC/SIMPUT'
    prihdr['HDUCLAS1'] = 'IMAGE'
    prihdr['HDUVERS'] = '1.1.0'
    prihdr['EXTNAME'] = 'IMAGE'
    prihdr['CTYPE1'] = ('RA---TAN', 'first axis (column) is Right Ascension')
    prihdr['CRPIX1'] = ((n_pixel+1) / 2., 'middle pixel of array in col direction')
    prihdr['CRVAL1'] = (0, 'RA of this middle pixel, in degrees')
    prihdr['CDELT1'] = (-angularSize_per_pixel/60., 'move 1column forward,decrease RA by CDELT1/deg')
    prihdr['CROTA1'] = 0
    prihdr['CUNIT1'] = 'deg'
    prihdr['CTYPE2'] = ('DEC--TAN', 'first axis (column) is Declination')
    prihdr['CRPIX2'] = ((n_pixel+1) / 2., 'middle pixel of array in row direction')
    prihdr['CRVAL2'] = (0, 'DEC of this middle pixel, in degrees')
    prihdr['CDELT2'] = (angularSize_per_pixel/60., 'move 1column forward,increase Dec by CDELT1/deg')
    prihdr['CROTA2'] = 0
    prihdr['CUNIT2'] = 'deg'
    prihdr['EQUINOX'] = 2000
    prihdr['RADECSYS'] = 'FK5'
    prihdr['aMinpPix'] = angularSize_per_pixel
    prihdu = fits.PrimaryHDU(matrix, header=prihdr)
    if os.path.isfile(out):
        os.remove(out)
    prihdu.writeto(out, overwrite=True)


nl = lambda sel : len(sel.nonzero()[0])

LC_dir = 'LCerass'

sky_map_hdu = Table.read(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS.fits') )

# merge catalog
# for srv_val in sky_map_hdu['SRVMAP'][(sky_map_hdu['OWNER']==1)]:
for srv_val in sky_map_hdu['SRVMAP'][(sky_map_hdu['OWNER']==2)|(sky_map_hdu['OWNER']==0)][:1]:
    t0 = time.time()
    str_field = str(srv_val).zfill(6)
    t_in = Table.read( os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'Xgas_bHS0.8_simput.fits') )
    # re-writes the images
    for clu in t_in:
        p2_im_out = os.path.join(os.environ['UCHUU'], LC_dir, str_field, clu["IMAGE"][:-7])
        out_dir = os.path.dirname(p2_im_out)
        os.system('mkdir -p '+out_dir)
        p2_im_in = os.path.join(os.environ['UCHUU'], clu["IMAGE"][:-7])
        im = fits.open(p2_im_in)[0]
        n_pixel = len(im.data)
        matrix = im.data
        write_img(matrix, out=p2_im_out, n_pixel =n_pixel, angularSize_per_pixel = 0.033, RA_CEN=clu['RA'], DEC_CEN=clu['DEC'])
        print(p2_im_out)

