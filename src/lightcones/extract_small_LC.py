"""
 * 2 square degrees e.g. COSMOS selected by 149.2<ra<151.1 && 1.3<Dec<3.1 up to redshift 6
 * 60 square degrees e.g. GAMA 9h selected by 129<ra<141 && -2<Dec<3      up to redshift 6
 * 1800 square degrees e.g. fat equatorial stripe (HSC Wide) selected by abs(dec)<2.5 & All R.A., likely needs to be re-cut in redshift to be load-able in memory.

"""
import time
t0 = time.time()

import os, glob, sys
from astropy.table import Table
import numpy as np

nl = lambda sel : len(sel.nonzero()[0])

LC_dir = sys.argv[1]

all_z_dirs = np.array([  'z0p00',
                    'z0p02',
                    'z0p05',
                    'z0p09',
                    'z0p14',
                    'z0p19',
                    'z0p25',
                    'z0p30',
                    'z0p36',
                    'z0p43',
                    'z0p49',
                    'z0p56',
                    'z0p63',
                    'z0p70',
                    'z0p78',
                    'z0p86',
                    'z0p94',
                    'z1p03',
                    'z1p12',
                    'z1p22',
                    'z1p32',
                    'z1p43',
                    'z1p54',
                    'z1p65',
                    'z1p77',
                    'z1p90',
                    'z2p03',
                    'z2p17',
                    'z2p31',
                    'z2p46',
                    'z2p62',
                    'z2p78',
                    'z2p95',
                    'z3p13',
                    'z3p32',
                    'z3p61',
                    'z3p93',
                    'z4p27',
                    'z4p63',
                    'z5p15',
                    'z5p73' ])

if LC_dir == 'LC0002' :
    ra0 = 149.2
    ra1 = 151.1
    de0 = 1.3
    de1 = 3.1
    z_dirs = all_z_dirs
if LC_dir == 'LC0060' :
    ra0 = 129
    ra1 = 141
    de0 = -2
    de1 = 3
    z_dirs = all_z_dirs[:37]
if LC_dir == 'LC1800' :
    ra0 = 0.0
    ra1 = 360.
    de0 = -2.5
    de1 = 2.5
    z_dirs = all_z_dirs[:27]

for z_dir in z_dirs:
    print('='*100)
    print(z_dir)
    print('='*100)
    p_2_catalogues = np.array( glob.glob( os.path.join(os.environ['UCHUU'], 'GPX8', z_dir, 'replication_*_*_*', 'glist.fits') ) )
    p_2_catalogues.sort()
    for p_2_catalogue in p_2_catalogues:
        from astropy.io import fits
        hdu1 = fits.open(p_2_catalogue)[1]
        selection = (hdu1.data['RA']>=ra0) & (hdu1.data['RA']<=ra1)&(hdu1.data['DEC']>=de0) & (hdu1.data['DEC']<=de1)
        N_selected = nl(selection)
        if N_selected>0:
            t_in = Table.read(p_2_catalogue)
            dir_4_out = os.path.join(os.environ['UCHUU'], LC_dir, z_dir , p_2_catalogue.split('/')[-2])
            os.system('mkdir -p ' + dir_4_out)
            p_2_catalogue_out = os.path.join( dir_4_out, 'glist.fits')
            t_in[selection].write(p_2_catalogue_out, overwrite = True)
            print(N_selected, p_2_catalogue_out, 'written')

