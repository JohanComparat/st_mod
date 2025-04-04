import time
t0 = time.time()
import os, glob, sys
from astropy.table import Table, vstack
import numpy as np

nl = lambda sel : len(sel.nonzero()[0])

LC_dir = 'LCerass'

sky_map_hdu = Table.read(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS.fits') )

# merge catalog
for srv_val in sky_map_hdu['SRVMAP']:
    print(srv_val)
    str_field = str(srv_val).zfill(6)
    all_tile_catalogues = np.array( glob.glob( os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'z?p??', 'replication_*', 'Xgas_bHS0.8_simput.fits') ) )
    print(len(all_tile_catalogues), all_tile_catalogues)
    command = 'rm -rf ' + os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'z?p??')
    print(command)
    os.system( command )
    print('time spent=', time.time()-t0)
