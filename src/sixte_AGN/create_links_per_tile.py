
import time
t0 = time.time()
import os, glob, sys
from astropy.table import Table, vstack
import numpy as np

LC_dir = 'LCerass'

sky_map_hdu = Table.read(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS.fits') )

# merge catalog
for srv_val in sky_map_hdu['SRVMAP']:
    print(srv_val)
    str_field = str(srv_val).zfill(6)
    print(os.path.join(os.environ['UCHUU'], LC_dir, str_field))
    os.chdir(os.path.join(os.environ['UCHUU'], LC_dir, str_field))
    command2 = 'ln -s /home/idies/workspace/erosim/simput/AGNspectra_V2 .'
    print(command2)
    os.system(command2)
