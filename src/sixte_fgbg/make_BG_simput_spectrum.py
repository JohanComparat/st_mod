import os, glob, sys
from scipy.stats import scoreatpercentile
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as n
from scipy.interpolate import interp1d
from astropy.table import Table, Column


BG_out_file = os.path.join('/data26s/mpecl/eRASS1', 'E_hist_eRASS1_BG_evt.fits')
print('histogram', BG_out_file)
bg = fits.open(BG_out_file)[1].data
energy = (bg['E_MIN']+bg['E_max'])/2./1000.
fluxes = bg['H_e_t']

path_2_SMPT_dir = os.path.join('/data40s/erosim/eRASS', 'bkg_erosita_simput_full_sky', 'spectra')
if os.path.isdir(path_2_SMPT_dir) == False:
	os.system('mkdir -p ' + path_2_SMPT_dir)

p_2_spectrum = os.path.join(path_2_SMPT_dir, 'mean_bg_spectrum.fits' )
print('output', p_2_spectrum)

Col_E = fits.column.Column(array=n.array([energy, ], dtype=n.object), name='ENERGY', format='PE()', unit='keV')
Col_F = fits.column.Column(	array=n.array([	fluxes, ], dtype=n.object), name='FLUXDENSITY',	format='PE()',	unit='photon/s/cm**2/keV')
HDU1 = fits.BinTableHDU.from_columns([Col_E, Col_F])
HDU1.header['EXTNAME'] = 'SPECTRUM'
HDU1.header['HDUCLAS1'] = 'SPECTRUM'
HDU1.writeto(p_2_spectrum, overwrite=True)
print('done')
