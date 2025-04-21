import time
t0 = time.time()
import os, glob, sys
from astropy.table import Table, vstack
import numpy as np

nl = lambda sel : len(sel.nonzero()[0])

log10FXmin = -15.7
jj_zdir = int(sys.argv[1])

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

sky_map_hdu = Table.read(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS.fits') )
def get_srvmap(ra, dec):
    return sky_map_hdu['SRVMAP'].value[(sky_map_hdu['RA_MIN']<ra ) & ( sky_map_hdu['RA_MAX'] >= ra ) & ( sky_map_hdu['DE_MIN']<dec ) & ( sky_map_hdu['DE_MAX'] >= dec)]

LC_dir = 'LCerass'
z_dirs = all_z_dirs#[:40]

N_tot=[]
z_dir = z_dirs[jj_zdir]
print('='*100)
print(z_dir)
print('='*100)
p_2_catalogues = np.array( glob.glob( os.path.join(os.environ['UCHUU'], 'FullSky', z_dir, 'replication_*_*_*', 'glist.fits') ) )
p_2_catalogues.sort()
for p_2_catalogue in p_2_catalogues:
    t_in = Table.read(p_2_catalogue)
    #selection = (t_in['FLUX']>=log10FXmin)
    #t_in = t_in[selection]
    if len(t_in)>0:
        t_in['RA'][t_in['RA']==360.]=359.9999
        t_in['SRVMAP'] = np.array([get_srvmap(ra, dec)[0] for (ra, dec) in zip(t_in['RA'], t_in['DEC']) ])
        U_srvmap_val = np.unique(t_in['SRVMAP'])
        for srv_val in U_srvmap_val:
            dir_4_out = os.path.join(os.environ['UCHUU'], LC_dir, str_field, z_dir , p_2_catalogue.split('/')[-2])
            os.system('mkdir -p ' + dir_4_out)
            p_2_catalogue_out = os.path.join( dir_4_out, 'glist.fits')
            if os.path.isfile(p_2_catalogue_out)==False:
                s2 = t_in['SRVMAP']==srv_val
                str_field = str(srv_val).zfill(6)
                N_selected = nl(s2)
                t_in[s2].write(p_2_catalogue_out, overwrite = True)
                print(N_selected, p_2_catalogue_out, 'written')
                N_tot.append(N_selected)

N_haloes =  np.sum(N_tot)
print(N_haloes, 'to be simulated')
