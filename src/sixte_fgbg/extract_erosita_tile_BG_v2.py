import time
t0 = time.time()
import os, glob, sys
from astropy.table import Table, vstack
import numpy as np
import healpy
nl = lambda sel : len(sel.nonzero()[0])

sky_map_hdu = Table.read(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS.fits') )
def get_srvmap(ra, dec):
    return sky_map_hdu['SRVMAP'].value[(sky_map_hdu['RA_MIN']<ra ) & ( sky_map_hdu['RA_MAX'] >= ra ) & ( sky_map_hdu['DE_MIN']<dec ) & ( sky_map_hdu['DE_MAX'] >= dec)]

#/home/idies/workspace/erosim/eRASS8_events_BG_FG/???/t0erass_ccd1_evt.fits | wc -l

LC_dir = 'LCerass'
pbg_dir = 'pBG'
pbg_dir = 'pBG2'

N_tot=[]
#z_dir = z_dirs[jj_zdir]
#print('='*100)
#print(z_dir)
#print('='*100)
p_2_catalogues = np.array( glob.glob( os.path.join(os.environ['HOME'], 'workspace/erosim/eRASS8_events_BG_FG', '???', 't0erass_ccd?_evt.fits') ) )
p_2_catalogues.sort()
for p_2_catalogue in p_2_catalogues[4510:]:
    t0 = time.time()
    print('Reads', p_2_catalogue)
    hpx_val = int(p_2_catalogue.split('/')[-2])
    t_in = Table.read(p_2_catalogue)
    t_in['HEALPIX_VAL'] = healpy.ang2pix(8, np.pi/2. - t_in['DEC'] *np.pi/180. , t_in['RA']*np.pi/180. , nest=True)
    t_in = t_in[t_in['HEALPIX_VAL']==hpx_val]
    t_in['RA'][t_in['RA']==360.]=359.9999
    t_in['RA'][t_in['RA']==0.]=0.000001
    t_in['DEC'][t_in['DEC']==90.]=89.99999
    t_in['DEC'][t_in['DEC']==-90.]=-89.99999
    t_in['SRVMAP'] = np.array([get_srvmap(ra, dec)[0] for (ra, dec) in zip(t_in['RA'], t_in['DEC']) ])
    U_srvmap_val = np.unique(t_in['SRVMAP'])
    print(len(U_srvmap_val), 'fields')
    for srv_val in U_srvmap_val:
        str_field = str(srv_val).zfill(6)
        dir_4_out = os.path.join(os.environ['UCHUU'], LC_dir, str_field, pbg_dir )
        os.system('mkdir -p ' + dir_4_out)
        p_2_catalogue_out = os.path.join( dir_4_out,  'pix_'+str(hpx_val).zfill(3)+'_'+os.path.basename(p_2_catalogue))
        if os.path.isfile(p_2_catalogue_out)==False:
            s2 = t_in['SRVMAP']==srv_val
            N_selected = nl(s2)
            t_in[s2].write(p_2_catalogue_out, overwrite = True)
            print(N_selected, p_2_catalogue_out, 'written', time.time()-t0)
        else:
            print(os.path.isfile(p_2_catalogue_out), 'done', time.time()-t0)
