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
#from models import AGN as GG

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

for z_dir in z_dirs[::-1]:
    print('='*100)
    print(z_dir)
    print('='*100)
    all_reps = np.array( glob.glob( os.path.join(os.environ['UCHUU'], LC_dir, z_dir, 'replication_*_*_*') ))
    rps_ids = np.array([el.split('/')[-1] for el in all_reps])
    print(rps_ids)
    p_2_catalogues = np.array([ os.path.join(os.environ['UCHUU'], 'GPX8', z_dir, el, 'glist.fits')  for el in rps_ids ])
    p_2_catalogues.sort()
    #C_AGN = GG.AGN(z_dir, LC_dir=LC_dir, scatter_0 = 0.8, f_sat = 0 )
    p_2_meta = os.path.join(os.environ['UCHUU'], 'area_per_replica_'+z_dir+'.fits')
    p_2_meta_out = os.path.join(os.environ['UCHUU'], LC_dir, 'area_per_replica_'+z_dir+'.fits')
    t_meta = Table.read( p_2_meta )
    t_meta['N_frac_'+LC_dir] = 0.
    #t_meta['N_frac_FullSky'] = 1.
    for p_2_catalogue in p_2_catalogues:
        print(p_2_catalogue)
        from astropy.io import fits
        hdu1 = fits.open(p_2_catalogue)[1]
        selection = (hdu1.data['RA']>=ra0) & (hdu1.data['RA']<=ra1)&(hdu1.data['DEC']>=de0) & (hdu1.data['DEC']<=de1)
        N_selected = nl(selection)
        if N_selected>0:
            N_total = len(hdu1.data['RA'])
            replication_dir = p_2_catalogue.split('/')[-2]
            print(replication_dir)
            ix = float(replication_dir.split('_')[1])
            iy = float(replication_dir.split('_')[2])
            iz = float(replication_dir.split('_')[3])
            #print(ix, iy, iz)
            selection_area = ( t_meta['jx'] == ix ) & ( t_meta['jy'] == iy ) & ( t_meta['jz'] == iz )
            t_meta['N_frac_'+LC_dir][selection_area] = N_selected*1./N_total
            print(t_meta['N_frac_'+LC_dir][selection_area])
    print(t_meta['N_frac_'+LC_dir])
    t_meta.write(p_2_meta_out, overwrite = True)

