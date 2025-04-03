import time
#t0 = time.time()
import os, glob, sys
from astropy.table import Table, vstack
import numpy as np
import astropy.io.fits as fits

nl = lambda sel : len(sel.nonzero()[0])

LC_dir = 'LCerass'

sky_map_hdu = Table.read(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS.fits') )

# merge catalog
for srv_val in sky_map_hdu['SRVMAP'][(sky_map_hdu['OWNER']==2)|(sky_map_hdu['OWNER']==0)]:
    t0 = time.time()
    str_field = str(srv_val).zfill(6)
    folder_erosim = os.path.join(os.environ['UCHUU'], LC_dir, str_field)
    folder_erosim_s4 = os.path.join(os.environ['UCHUU'], LC_dir, str_field, 's4_c030')
    folder_erosim_s5 = os.path.join(os.environ['UCHUU'], LC_dir, str_field, 's5_c030')
    folder_data_s4 = os.path.join(os.environ['HOME'], 'workspace/Storage/comparat/persistent/data/data_s4_c030', str(srv_val), 'c030', '*')
    folder_data_s5 = os.path.join(os.environ['HOME'], 'workspace/Storage/comparat/persistent/data/data_s5_c030', str(srv_val), 'c030', '*')
    mk_s4 = 'mkdir -p '+folder_erosim_s4
    mk_s5 = 'mkdir -p '+folder_erosim_s5
    cp_s4 = 'rsync -avz '+folder_data_s4+' '+folder_erosim_s4+'/'
    cp_s5 = 'rsync -avz '+folder_data_s5+' '+folder_erosim_s5+'/'
    print(mk_s4)
    print(mk_s5)
    print(cp_s4)
    print(cp_s5)
    os.system(mk_s4)
    os.system(mk_s5)
    os.system(cp_s4)
    os.system(cp_s5)
    folder_erosim_s4 = os.path.join(os.environ['UCHUU'], LC_dir, str_field, 's4_eSASS')
    folder_erosim_s5 = os.path.join(os.environ['UCHUU'], LC_dir, str_field, 's5_eSASS')
    folder_data_s4 = os.path.join(os.environ['HOME'], 'workspace/Storage/comparat/persistent/data/data_s4_c030', str(srv_val), 'eSASS', '*')
    folder_data_s5 = os.path.join(os.environ['HOME'], 'workspace/Storage/comparat/persistent/data/data_s5_c030', str(srv_val), 'eSASS', '*')
    mk_s4 = 'mkdir -p '+folder_erosim_s4
    mk_s5 = 'mkdir -p '+folder_erosim_s5
    cp_s4 = 'rsync -avz '+folder_data_s4+' '+folder_erosim_s4+'/'
    cp_s5 = 'rsync -avz '+folder_data_s5+' '+folder_erosim_s5+'/'
    print(mk_s4)
    print(mk_s5)
    print(cp_s4)
    print(cp_s5)
    os.system(mk_s4)
    os.system(mk_s5)
    os.system(cp_s4)
    os.system(cp_s5)
