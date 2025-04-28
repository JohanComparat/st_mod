"""
on ds54	

python healpixelize_catalogue.py

create 768 catalogues: one per NSIDE 8  
read the catalogue at NSIDE = 64 = 2**6

4.7e-14 erg/cm2/s on eFEDs / GAMA field

mean over efeds : 
2radii 3.75e-11, divide by 797.87
1radii 9.53e-12, divide by 202.76
conversion factor 0.2-2.3 : 1.092e+12

"""
import sys, os, time
t0=time.time()
import numpy as n
from astropy_healpix import healpy
from astropy.table import Table, Column
import astropy.units as u
from astropy.coordinates import SkyCoord
import subprocess
import astropy.io.fits as fits
#from astropy import wcs

NSIDE = 128 

p_2_map = os.path.join('/data26s/mpecl/eRASS1', 'maps', 'BG_map_eRASS1_june2021', str(NSIDE), 'map_C1.fits')
# 0.3-2.3 keV 

path_2_SMPT_dir = os.path.join('/data40s/erosim/eRASS', 'bkg_erosita_simput_full_sky')
#path_2_SMPT_dir = '/data40s/erosim/eRASS/bkg_erosita_simput_full_sky'

N_pixels_all = healpy.nside2npix(NSIDE)
area_per_pixel =  healpy.nside2pixarea(NSIDE, degrees=True)
radius = (area_per_pixel/n.pi)**0.5
pix_ids = n.arange(N_pixels_all)

theta, phi = healpy.pix2ang(NSIDE, pix_ids, nest=True)
dec_pix = (n.pi/2.-theta)*180/n.pi
ra_pix = phi*180/n.pi

# conversion to galactic and ecliptic coordinates

coords = SkyCoord(ra_pix, dec_pix, unit='deg', frame='icrs')
g_lat_pix = coords.galactic.b.value
g_lon_pix = coords.galactic.l.value

ecl_lat_pix = coords.barycentrictrueecliptic.lat.value
ecl_lon_pix = coords.barycentrictrueecliptic.lon.value

g_lon_sym = 180 - (g_lon_pix - 180)
g_lat_sym = g_lat_pix
coords_sym = SkyCoord(g_lon_sym, g_lat_sym, unit='deg', frame='galactic')

pix_ids_sym = healpy.ang2pix(NSIDE, n.pi/2. - coords_sym.icrs.dec.value*n.pi/180., coords_sym.icrs.ra.value*n.pi/180., nest=True)

#print( len(n.unique(g_lon_pix[g_lon_pix<180])), len(g_lon_pix[g_lon_pix<180]) )

halfSky = fits.open(p_2_map)[1].data

sym_rates = halfSky['count_rates']
sym_rates[g_lon_pix<180] = halfSky['count_rates'][pix_ids_sym[g_lon_pix<180]]

#t = Table()
#t.add_column(Column(name='count_rates', data=sym_rates, unit='s**-1',dtype=n.float32 ) )
#t.write( os.path.join(directory, 'SYM_map_C1.fits'), overwrite=True)

# ECF: in the band 0.3-2.3 1.092e+12
#energy_conversion = 11.2 / 1.092e+12 
energy_conversion = 1.0 / ( 1.092e+12 )

images = n.array(['images/image_'+str(pix).zfill(6)+'.fits[IMAGE]' for pix in pix_ids ])
spectra = n.array(['spectra/mean_bg_spectrum.fits'+ """[SPECTRUM][#row==1]""" for pix in pix_ids ])

# makes the complete table
allsky = Table()
allsky.add_column(Column(name='SRC_ID'   ,  data = pix_ids ) )
allsky.add_column(Column(name='RA'       ,  data=ra_pix, unit='deg' ) )
allsky.add_column(Column(name='DEC'      ,  data=dec_pix, unit='deg' ) )
allsky.add_column(Column(name='E_MIN'    ,  data=n.ones_like(ra_pix)*0.3, unit='keV' ) )
allsky.add_column(Column(name='E_MAX'    ,  data=n.ones_like(ra_pix)*2.3, unit='keV' ) )
allsky.add_column(Column(name='FLUX'     ,  data=sym_rates*energy_conversion, unit='erg * s**(-1) * cm**(-2)' ) )
allsky.add_column(Column(name='SPECTRUM' ,  data=spectra, unit='' ) )
allsky.add_column(Column(name='IMAGE'    ,  data=images, unit='' ) )
allsky.write(os.path.join(path_2_SMPT_dir, 'catalogue.fits'), overwrite = True )


def get_pp1(p):
	return n.array([4*p, 4*p + 1, 4*p + 2, 4*p + 3])

def write_simput_file(pix_id_8 = 0):
	#pix_id_8 = 0
	path_2_SMPT_pixel_catalog =os.path.join(path_2_SMPT_dir, 'SIMPUT_'+str(pix_id_8).zfill(6)+'.fits')
	# get the set of pixels
	# pix_id_array_64 = get_pp1( get_pp1( get_pp1( pix_id_8 ) ) ).ravel()
	pix_id_array_128 = get_pp1( get_pp1(get_pp1(get_pp1(pix_id_8)))).ravel()
	# selection = n.in1d(allsky['SRC_ID'], pix_id_array_64, assume_unique=True)
	selection = n.in1d(allsky['SRC_ID'], pix_id_array_128, assume_unique=True)
	hdu_cols = fits.ColDefs([
		fits.Column(name="SRC_ID"  , format='K'   , unit='',            array=allsky["SRC_ID"  ][selection] ) ,
		fits.Column(name="RA"      , format='D'   , unit='deg',         array=allsky["RA"      ][selection] ) ,
		fits.Column(name="DEC"     , format='D'   , unit='deg',         array=allsky["DEC"     ][selection] ) ,
		fits.Column(name="E_MIN"   , format='D'   , unit='keV',         array=allsky["E_MIN"   ][selection] ) ,
		fits.Column(name="E_MAX"   , format='D'   , unit='keV',         array=allsky["E_MAX"   ][selection] ) , 
		fits.Column(name="FLUX"    , format='D'   , unit='erg/s/cm**2', array=allsky["FLUX"    ][selection] ) ,
		fits.Column(name="SPECTRUM", format='100A', unit='',            array=allsky["SPECTRUM"][selection] ) ,
		fits.Column(name="IMAGE"   , format='100A', unit='',            array=allsky["IMAGE"   ][selection] ) 
		])
	hdu = fits.BinTableHDU.from_columns(hdu_cols)
	hdu.name = 'SRC_CAT'
	hdu.header['HDUCLASS'] = 'HEASARC/SIMPUT'
	hdu.header['HDUCLAS1'] = 'SRC_CAT'
	hdu.header['HDUVERS'] = '1.1.0'
	hdu.header['RADESYS'] = 'FK5'
	hdu.header['EQUINOX'] = 2000.0

	outf = fits.HDUList([fits.PrimaryHDU(), hdu])  # ,  ])
	if os.path.isfile(path_2_SMPT_pixel_catalog):
		os.system("rm " + path_2_SMPT_pixel_catalog)
	outf.writeto(path_2_SMPT_pixel_catalog, overwrite=True)
	print(path_2_SMPT_pixel_catalog, 'written', time.time() - t0)

for pix_id_8 in n.arange(healpy.nside2npix(8)):
	write_simput_file(pix_id_8)
