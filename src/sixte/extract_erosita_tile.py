import time
t0 = time.time()
import os, glob, sys
from astropy.table import Table, vstack
import numpy as np

nl = lambda sel : len(sel.nonzero()[0])

tile_ii = int(sys.argv[1])
#tile_ii= 0

log10FXmin = -15.7

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

tile = Table.read(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS.fits') )[tile_ii]
str_field = str(tile['SRVMAP']).zfill(6)
LC_dir = 'LCerass'

ra0 = tile['RA_MIN']
ra1 = tile['RA_MAX']
de0 = tile['DE_MIN']
de1 = tile['DE_MAX']
z_dirs = all_z_dirs[:23]

N_tot=[]
for z_dir in z_dirs:
    print('='*100)
    print(z_dir)
    print('='*100)
    p_2_catalogues = np.array( glob.glob( os.path.join(os.environ['UCHUU'], 'FullSky', z_dir, 'replication_*_*_*', 'Xgas_bHS0.8_simput.fits') ) )
    p_2_catalogues.sort()
    for p_2_catalogue in p_2_catalogues:
        t_in = Table.read(p_2_catalogue)
        selection = (t_in['FLUX']>=log10FXmin) & (t_in['RA']>=ra0) & (t_in['RA']<=ra1)&(t_in['DEC']>=de0) & (t_in['DEC']<=de1)
        N_selected = nl(selection)
        print(N_selected, 'in', p_2_catalogue)
        if N_selected>0:
            dir_4_out = os.path.join(os.environ['UCHUU'], LC_dir, str_field, z_dir , p_2_catalogue.split('/')[-2])
            os.system('mkdir -p ' + dir_4_out)
            p_2_catalogue_out = os.path.join( dir_4_out, 'Xgas_bHS0.8_simput.fits')
            t_in[selection].write(p_2_catalogue_out, overwrite = True)
            print( p_2_catalogue_out, 'written')
            N_tot.append(N_selected)

N_haloes =  np.sum(N_tot)
print(N_haloes, 'to be simulated')
# merge catalog
all_tile_catalogues = np.array( glob.glob( os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'z?p??', 'replication_*', 'Xgas_bHS0.8_simput.fits') ) )
print('merging',len(all_tile_catalogues), 'catalogs')
full_cat = []
for el in all_tile_catalogues:
    full_cat.append(Table.read(el))

merge_cat = vstack(( full_cat ))
merge_cat.write(os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'Xgas_bHS0.8_simput.fits'), overwrite = True )
print(os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'Xgas_bHS0.8_simput.fits'), 'written, time spent=', time.time()-t0)

os.system( 'rm -rf ' + os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'z?p??') )
print('time spent=', time.time()-t0)
