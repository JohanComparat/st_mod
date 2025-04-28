"""

"""
import os, glob, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm 
from matplotlib import rc
from scipy.stats import scoreatpercentile
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import healpy as hp

def mk_radec_plot(hpx_map, nside=32, filename=' ', min_hpx=0.00001, max_hpx=10, clabel='count', title=" ",  titleTop = " "):
	plt.figure(0, (14,10))
	rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
	cmap = plt.cm.inferno
	cmap.set_bad('w')
	cmap.set_under('w')
	# make the map with meridians overplotted
	hp.mollview(hpx_map, 
				coord='C', rot=(180,0), nest=True,
				unit=r'Count rate [s$^{-1}$]', xsize = 1000, cmap=cmap, min=min_hpx, max=max_hpx, norm='log', cbar=None)
	hp.graticule(dpar=15, dmer=20, verbose= True, color='white')
	# HA labels
	params = {'mathtext.default': 'regular' }          
	plt.rcParams.update(params)
	#hp.projtext(20,  32,  '$\mathrm{\mathsf{0^h}}$',     color = 'black', fontsize = 16, lonlat=True)
	#hp.projtext(88,  32,  '$\mathrm{\mathsf{6^h}}$',     color = 'black', fontsize = 16, lonlat=True)
	#hp.projtext(177, 32,  '$\mathrm{\mathsf{12^h}}$',    color = 'black', fontsize = 16, lonlat=True) #rotation
	#hp.projtext(264, 32,  '$\mathrm{\mathsf{18^h}}$',    color = 'black', fontsize = 16, lonlat=True)
	#hp.projtext(351, 32,  '$\mathrm{\mathsf{24^h}}$',    color = 'black', fontsize = 16, lonlat=True)
	hp.projtext(20,  32,  '$\mathrm{\mathsf{0^\circ}}$',     color = 'white', fontsize = 14, lonlat=True)
	hp.projtext(88,  32,  '$\mathrm{\mathsf{90^\circ}}$',     color = 'white', fontsize = 14, lonlat=True)
	hp.projtext(177, 32,  '$\mathrm{\mathsf{180^\circ}}$',    color = 'white', fontsize = 14, lonlat=True) #rotation
	hp.projtext(264, 32,  '$\mathrm{\mathsf{270^\circ}}$',    color = 'white', fontsize = 14, lonlat=True)
	hp.projtext(351, 32,  '$\mathrm{\mathsf{360^\circ}}$',    color = 'white', fontsize = 14, lonlat=True)
	plt.title(titleTop, size=10)
	ax = plt.gca()
	im = ax.get_images()[0]
	# DEC labels
	fig= plt.gcf()
	for ax in fig.get_axes():
		if type(ax) is hp.projaxes.HpxMollweideAxes:
			ax.set_ylim(-1, 0.51)
			ax.set_position([0.02, 0.03, 0.94, 0.95])
			ax.annotate(r'$\, \mathrm{\mathsf{ 0^\circ}}$', xy=(2.04,-0.02),    size=14, annotation_clip=False) 
			ax.annotate(r'$\, \mathrm{\mathsf{ 30^\circ}}$', xy=(1.90,0.38),     size=14) 
			ax.annotate(r'$\, \mathrm{\mathsf{-30^\circ}}$', xy=(1.86,-0.42),    size=14) 
			ax.annotate(r'$\, \mathrm{\mathsf{-60^\circ}}$', xy=(1.38, -0.78),   size=14) 
			ax.annotate(title, xy=(-2.0, -0.92),  size=14) # optional, can be removed

	# create colour bar
	cbaxes = fig.add_axes([0.1, 0.15, 0.8, 0.04]) # [left, bottom, width, height]
	tks=10**np.arange(-4,4,0.1) #[1e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000, 10000]
	cb = plt.colorbar(im,  orientation='horizontal', cax = cbaxes, ticks=tks)             
	cb.set_label(clabel, fontsize=8)
	cb.ax.tick_params(labelsize=8)
	#plt.xlim([0,1])
	# save plot to PDF file
	#plt.savefig(filename + '.pdf')
	plt.savefig(filename + '.png')


if source == 'BG':
	directory = '/home/comparat/data/erosita/eRASS/BGmaps/c946/unique_events_source_excluded'
	dir_2_maps = os.path.join(directory, NSIDE )
	fig_dir = os.path.join(os.environ['GIT_SXRBG'], 'figures', 'BG_map', str(NSIDE))

if source == 'SC':
	directory = '/home/comparat/data/erosita/eRASS/BGmaps/c946/unique_events_source_included'
	dir_2_maps = os.path.join(directory, NSIDE )
	fig_dir = os.path.join(os.environ['GIT_SXRBG'], 'figures', 'SC_map', str(NSIDE))

if os.path.isdir(fig_dir) == False:
	os.system('mkdir -p ' + fig_dir)

#path_2_maps = np.array( glob.glob( os.path.join(dir_2_maps, 'map_E_??.fits' ) ) )
#path_2_maps.sort()

#path_2_outs = np.array([os.path.join( fig_dir, os.path.basename( el )[:-5] ) for el in path_2_maps ])
#path_2_outs_CR = np.array([os.path.join( fig_dir, 'CR_' + os.path.basename(  el )[:-5] ) for el in path_2_maps ])

# count_rates
## counts
#for filename, p2map in zip(path_2_outs, path_2_maps): 
	#out_hpx = fits.open(p2map)[1].data['counts'] 
	#min_hpx = np.min(out_hpx[out_hpx>0])
	#out_hpx[out_hpx<min_hpx] = min_hpx/100
	#mk_radec_plot(out_hpx, nside=int(NSIDE), filename=filename, min_hpx= 2, clabel='counts')

#for filename, p2map in zip(path_2_outs_CR, path_2_maps): 
	#out_hpx = fits.open(p2map)[1].data['count_rates'] 
	#no_zero_map = out_hpx[out_hpx>0]
	#min_hpx = np.min(no_zero_map)
	#Val01, Val99 = scoreatpercentile(, no_zero_map, [1, 99])
	#out_hpx[out_hpx<min_hpx] = min_hpx/100
	#mk_radec_plot(out_hpx, nside=int(NSIDE), filename=filename, min_hpx= 0.001, clabel='count rates')


# merging energy bins 
#NSIDEs = np.array(['16', '32', '64', '128', '256'])
#energy_bins = np.hstack(( np.arange(0.05, 1.05, 0.1), np.arange(1,11,1) ))
#E_min = data=energy_bins[:-1]
#E_max = data=energy_bins[1:] 
#for ii, (e1, e2) in enumerate(zip(E_min, E_max)):
	#print(ii, np.round(e1,2), np.round(e2,2))
	
# Band definition from M. Freyberg email

bands = np.array([  
	["C1" ,    300 , 2300  , " "                                                                                                      ],
	["C2" ,    600 , 2300  , " "                                                                                                      ],
	["B1" ,    200 ,  600  , " "                                                                                                      ],
	["B2" ,    600 , 1000  , " "                                                                                                      ],
	["B3" ,   1000 , 2300  , " "                                                                                                      ],
	["B4" ,   2300 , 5000  , " "                                                                                                      ],
	["E00",    100 ,  150  , " "                                                                                                      ],
	["E01",    150 ,  200  , " "                                                                                                      ],
	["E02",    200 ,  335  , "C v @ 308 eV"                                                                                           ],
	["E03",    335 ,  470  , "c vi @ 367 eV, N vi @ 431 eV"                                                                           ],
	["E04",    470 ,  600  , "N vii @ 500 eV, O vii @ 570 eV"                                                                         ],
	["E05",    600 ,  730  , "O viii @ 654 eV, O vii @ 666 and 698 eV"                                                                ],
	["E06",    730 ,  870  , "Fe xvii @ 826 eV"                                                                                       ],
	["E07",    870 , 1000  , "Ne ix @ 915 and 922 eV, Fe xx at 996 eV"                                                                ],
	["E08",   1000 , 1200  , "Ne x @ 1022 eV"                                                                                         ],
	["E09",   1200 , 1410  , "Mg xi @ 1340 and 1352 eV"                                                                               ],
	["E10",   1410 , 1600  , "Al-K @ 1489 eV, Mg xii @ 1471 eV"                                                                       ],
	["E11",   1600 , 1800  , " "                                                                                                      ],
	["E12",   1800 , 2000  , "Si xiii @ 1863 eV"                                                                                      ],
	["E13",   2000 , 2300  , " "                                                                                                      ],
	["E14",   2300 , 3100  , "S lines ?"                                                                                              ],
	["E15",   3100 , 3800  , "Ca-K @ 3.6 keV, Ar lines? S lines?"                                                                     ],
	["E16",   3800 , 5000  , "Ti-Ka @ 4.5 keV"                                                                                        ],
	["E17",   5000 , 5700  , "Cr-Ka @ 5.4 keV"                                                                                        ],
	["E18",   5700 , 7300  , "Fe-Ka @ 6.4 keV, Fe xxv @ 6.64, 6.68 6.70 keV Fe xxvi @ 6.93 keV, Fe-Kb @ 7.06 keV, Co-Ka @ 6.93 keV"   ],
	["E19",   7300 , 9000  , "Ni-Ka @ 7.48 keV Cu-Ka @ 8.04 keV, Zn-Ka @ 8.64 keV"                                                    ],
	["E20",   9000 , 200000, "pileup and background"                                                                                  ]
	])                                                                                                                                

def get_hpx(name_str, col_name = 'counts' ):
	p2map = os.path.join(dir_2_maps, 'map_'+name_str+'.fits' )
	out_hpx = fits.open(p2map)[1].data[col_name]
	min_hpx = np.min(out_hpx[out_hpx>0])
	#out_hpx[out_hpx<min_hpx] = min_hpx/100
	no_zero_map = out_hpx[out_hpx>0]
	Val01, Val99 = scoreatpercentile(no_zero_map, [1, 99])
	print(Val01, Val99)
	area = hp.nside2pixarea(int(NSIDE), degrees = True)*3600
	return out_hpx/area, Val01/area, Val99/area

if source == 'SC':
	sourceDef = 'photons within <20arcSeconds or 2xEXT or sources \n'
if source == 'BG':
	sourceDef = 'photons away from >30arcSeconds or 3xEXT of any source \n'

for el in bands:
	name_str, emin, emax, comment = el
	filename = os.path.join( fig_dir, 'E_countRate_band_'+name_str )
	print(filename)
	hpx, min_hpx, max_hpx = get_hpx(name_str, col_name = 'count_rates')
	mk_radec_plot(hpx, nside=int(NSIDE), filename=filename, min_hpx=min_hpx, max_hpx=max_hpx, clabel=r'count rate [$s^{-1}$ arcmin$^{-2}$]', title=name_str+' '+emin+'<E[eV]<'+emax, titleTop= sourceDef+emin+'<E[eV]<'+emax+'\n'+ comment)

EMIN = 10**np.arange(2,4,0.05)
EMAX = 10**np.arange(2+0.05,4+0.05,0.05)
band_names = np.array(['N'+str(np.round(el,2)) for el in np.arange(2,4,0.05) ])

for name_str, eminF, emaxF in zip(band_names, EMIN, EMAX):
	comment = " "
	emin, emax = str(np.round(eminF,1)), str(np.round(emaxF,1))
	filename = os.path.join( fig_dir, 'E_countRate_band_'+name_str )
	print(filename)
	hpx, min_hpx, max_hpx = get_hpx(name_str, col_name = 'count_rates')
	mk_radec_plot(hpx, nside=int(NSIDE), filename=filename, min_hpx=min_hpx, max_hpx=max_hpx, clabel=r'count rate [$s^{-1}$ arcmin$^{-2}$]', title=name_str+' '+emin+'<E[eV]<'+emax, titleTop= sourceDef+emin+'<E[eV]<'+emax+'\n'+ comment)
