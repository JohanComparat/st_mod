"""
cd $GIT_ERASS_SIM/sixte
python BG_map.py 64

"""
import numpy as n
import os, sys, glob
import astropy.io.fits as fits
from astropy.table import Table, Column
import healpy
from astropy_healpix import healpy
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as p
from astropy.wcs import WCS
from sklearn.neighbors import BallTree
import time
t0 = time.time()
import pandas

# parameters : 
# masking radius for point sources
NSIDE = 64 # int(sys.argv[1])
N_pixels_all = healpy.nside2npix(NSIDE)

energy_bins = n.hstack(( n.arange(0.05, 1.05, 0.1), n.arange(1,11,1) ))
energy_bins_min = energy_bins[:-1]
energy_bins_max = energy_bins[1:]

directory = '/data/data/erosim/'
p2_sourceExcl_cat = os.path.join( directory, 'all_t0erass_ccd1_evt.fits' )

dir_2_BGmaps = os.path.join(directory, str(NSIDE) )
fig_dir = os.path.join(os.environ['GIT_SXRBG'], 'figures', 'BG_map', 'simulation')
data_dir = os.path.join(os.environ['GIT_SXRBG'], 'data', 'simulated_BG_map')
if os.path.isdir(dir_2_BGmaps) == False:
	os.system('mkdir -p ' + dir_2_BGmaps)
	os.system('chmod ogu+rX '+dir_2_BGmaps)
	os.system('mkdir -p ' + fig_dir)
	os.system('chmod ogu+rX '+fig_dir)
	os.system('mkdir -p ' + data_dir)
	os.system('chmod ogu+rX '+data_dir)


new_bg = fits.open(p2_sourceExcl_cat)[1].data
t = Table()
t.add_column(Column(data=new_bg['RA'], name='RA', unit='deg'  ) )
t.add_column(Column(data=new_bg['DEC'], name='DEC', unit='deg'  ) )
t.add_column(Column(data=new_bg['SIGNAL']*1000, name='PHA', unit='keV'  ) )
evts = t.to_pandas()
E_ra     = evts.RA 
E_dec    = evts.DEC 
E_energy = evts.PHA

evts['HEALPIX_VAL'] = healpy.ang2pix(NSIDE, n.pi/2. - evts.DEC *n.pi/180. , evts.RA*n.pi/180. , nest=True)

out = evts.groupby(by='HEALPIX_VAL')
out2 = out.agg({'exposure_time' : ['mean', 'median']})
exposure_time_in_pixel = out2['exposure_time']['median'] 
U_pix, U_counts = n.unique(evts['HEALPIX_VAL'], return_counts=True)
CountRate = U_counts / exposure_time_in_pixel

counts      = n.zeros(N_pixels_all)
count_rates = n.zeros(N_pixels_all)
counts      [U_pix] = U_counts
count_rates [U_pix] = CountRate
path_2_out = os.path.join(dir_2_BGmaps, "map_E_all.fits" )

t = Table()
t.add_column(Column(name='counts'     , data=counts.astype('int'), unit='',dtype=n.int32 ) )
t.add_column(Column(name='count_rates', data=count_rates, unit='s**-1',dtype=n.float32 ) )
t.write(path_2_out, overwrite=True)


	
# Band definition from M. Freyberg email
"""
bands = n.array([  
	["C1" ,    300 , 2300],
	["C2" ,    600 , 2300],
	["B1" ,    200 ,  600],  
	["B2" ,    600 , 1000],
	["B3" ,   1000 , 2300],
	["B4" ,   2300 , 5000],
	["E00",    100 ,  150], # *\ could be merged, different joint lower threshold   
	["E01",    150 ,  200], # */
	["E02",    200 ,  335], #    C v @ 308 eV
	["E03",    335 ,  470], #    c vi @ 367 eV, N vi @ 431 eV
	["E04",    470 ,  600], #    N vii @ 500 eV, O vii @ 570 eV
	["E05",    600 ,  730], #    O viii @ 654 eV, O vii @ 666 and 698 eV
	["E06",    730 ,  870], #    Fe xvii @ 826 eV
	["E07",    870 , 1000], #    Ne ix @ 915 and 922 eV, Fe xx at 996 eV !!
	["E08",   1000 , 1200], #    Ne x @ 1022 eV
	["E09",   1200 , 1410], #    Mg xi @ 1340 and 1352 eV
	["E10",   1410 , 1600], #    Al-K @ 1489 eV, Mg xii @ 1471 eV
	["E11",   1600 , 1800], # 
	["E12",   1800 , 2000], #    Si xiii @ 1863 eV
	["E13",   2000 , 2300], # 
	["E14",   2300 , 3100], #    S lines ?
	["E15",   3100 , 3800], #    Ca-K @ 3.6 keV, Ar lines? S lines?
	["E16",   3800 , 5000], #    Ti-Ka @ 4.5 keV
	["E17",   5000 , 5700], #    Cr-Ka @ 5.4 keV
	["E18",   5700 , 7300], #    Fe-Ka @ 6.4 keV, Fe xxv @ 6.64, 6.68 6.70 keV Fe xxvi @ 6.93 keV, Fe-Kb @ 7.06 keV, Co-Ka @ 6.93 keV
	["E19",   7300 , 9000], #    Ni-Ka @ 7.48 keV Cu-Ka @ 8.04 keV, Zn-Ka @ 8.64 keV
	["E20",   9000 , 200000]   #    only useful for pileup and background
	])


for element in bands:
	BAND_NAME, E_min, E_max = element[0], float(element[1]), float(element[2])
	print(BAND_NAME, E_min, E_max)
	path_2_out = os.path.join(dir_2_BGmaps, 'map_'+BAND_NAME+ ".fits" )
	selection = (evts['PHA']>=E_min) & (evts['PHA']<E_max)
	U_pix, U_counts = n.unique(evts[selection]['HEALPIX_VAL'], return_counts=True)
	print(U_pix, U_counts)
	if len(U_pix)>5 :
		out = evts[selection].groupby(by='HEALPIX_VAL')
		out2 = out.agg({'exposure_time' : ['mean', 'median']})
		exposure_time_in_pixel = out2['exposure_time']['median'] 
		CountRate = U_counts / exposure_time_in_pixel
		counts      = n.zeros(N_pixels_all)
		count_rates = n.zeros(N_pixels_all)
		counts      [U_pix] = U_counts
		count_rates [U_pix] = CountRate
		t = Table()
		t.add_column(Column(name='counts'     , data=counts.astype('int'), unit='',dtype=n.int32 ) )
		t.add_column(Column(name='count_rates', data=count_rates, unit='s**-1',dtype=n.float32 ) )
		t.write(path_2_out, overwrite=True)
"""
emin = 10**n.arange(2,4,0.05)
emax = 10**n.arange(2+0.05,4+0.05,0.05)
band_names = n.array(['N'+str(n.round(el,2)) for el in n.arange(2,4,0.05) ])


for BAND_NAME, E_min, E_max in zip(band_names, emin, emax) :
	print(BAND_NAME, E_min, E_max)
	path_2_out = os.path.join(dir_2_BGmaps, 'map_'+BAND_NAME+ ".fits" )
	selection = (evts['PHA']>=E_min) & (evts['PHA']<E_max)
	U_pix, U_counts = n.unique(evts[selection]['HEALPIX_VAL'], return_counts=True)
	print(U_pix, U_counts)
	if len(U_pix)>3 :
		out = evts[selection].groupby(by='HEALPIX_VAL')
		out2 = out.agg({'exposure_time' : ['mean', 'median']})
		exposure_time_in_pixel = out2['exposure_time']['median'] 
		CountRate = U_counts / exposure_time_in_pixel
		counts      = n.zeros(N_pixels_all)
		count_rates = n.zeros(N_pixels_all)
		counts      [U_pix] = U_counts
		count_rates [U_pix] = CountRate
		t = Table()
		t.add_column(Column(name='counts'     , data=counts.astype('int'), unit='',dtype=n.int32 ) )
		t.add_column(Column(name='count_rates', data=count_rates, unit='s**-1',dtype=n.float32 ) )
		t.write(path_2_out, overwrite=True)
