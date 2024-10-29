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

LC_dir = 'FullSky'

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

for z_dir in all_z_dirs:
    print('='*100)
    print(z_dir)
    print('='*100)
    #C_AGN = GG.AGN(z_dir, LC_dir=LC_dir, scatter_0 = 0.8, f_sat = 0 )
    p_2_meta = os.path.join(os.environ['UCHUU'], 'area_per_replica_'+z_dir+'.fits')
    p_2_meta_out = os.path.join(os.environ['UCHUU'], LC_dir, 'area_per_replica_'+z_dir+'.fits')
    t_meta = Table.read( p_2_meta )
    t_meta['N_frac_'+LC_dir] = 1.
    print(t_meta['N_frac_'+LC_dir])
    t_meta.write(p_2_meta_out, overwrite = True)

