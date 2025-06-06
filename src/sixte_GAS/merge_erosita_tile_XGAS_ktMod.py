import time
#t0 = time.time()
import os, glob, sys
from astropy.table import Table, vstack
import numpy as np

nl = lambda sel : len(sel.nonzero()[0])

LC_dir = 'LCerass'
basename = sys.argv[1]

sky_map_hdu = Table.read(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS.fits') )

# merge catalog
#for srv_val in sky_map_hdu['SRVMAP'][(sky_map_hdu['OWNER']==1)]:
for srv_val in sky_map_hdu['SRVMAP'][(sky_map_hdu['OWNER']==2)|(sky_map_hdu['OWNER']==0)]:
    t0 = time.time()
    str_field = str(srv_val).zfill(6)
    all_tile_catalogues = np.array( glob.glob( os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'z?p??', 'replication_*', basename+'.fits') ) )
    print(str_field, 'merging',len(all_tile_catalogues), 'catalogs')
    if len(all_tile_catalogues)>=1:
        full_cat = []
        for el in all_tile_catalogues:
            full_cat.append(Table.read(el))

        merge_cat = vstack(( full_cat ))
        merge_cat.write(os.path.join(os.environ['UCHUU'], LC_dir, str_field, basename+'_XGAS.fits'), overwrite = True )
        print(os.path.join(os.environ['UCHUU'], LC_dir, str_field, basename+'_XGAS.fits'), 'written, time spent=', time.time()-t0)
