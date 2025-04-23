import subprocess
import os
import errno
import sys, glob
import healpy
import numpy as np
from astropy.io import fits
from astropy.table import Table, Column
sky_map_hdu = Table.read(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS.fits'))
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as plt

sky_tile_num = 164087 # int(sys.argv[1])
seed = 1
LC_dir = 'LCerass'
erass_option = "eRASS8"
print(seed, LC_dir, erass_option, sky_tile_num)#, env, erass_option)

sky_tile = sky_map_hdu[sky_map_hdu['SRVMAP']==sky_tile_num]
sky_tile_id = str(sky_tile_num)
str_field = sky_tile_id.zfill(6)
ra_cen = sky_tile['RA_CEN']
dec_cen = sky_tile['DE_CEN']
data_dir = os.path.join(os.environ['UCHUU'], LC_dir, str_field, erass_option + "_SEED_"+str(seed).zfill(3) +"_events_cluster_2025_04" )
print('events here',data_dir)
figure_dir = os.path.join(os.environ['UCHUU'], LC_dir, str_field, erass_option + "_SEED_"+str(seed).zfill(3) +"_events_cluster_2025_04", 'figures' )
print('figures here',figure_dir)
os.system('mkdir -p '+figure_dir )

ccd1 = fits.open(os.path.join( data_dir, 't0erass_ccd1_evt.fits' ))[1].data
ccd2 = fits.open(os.path.join( data_dir, 't0erass_ccd2_evt.fits' ))[1].data
ccd3 = fits.open(os.path.join( data_dir, 't0erass_ccd3_evt.fits' ))[1].data
ccd4 = fits.open(os.path.join( data_dir, 't0erass_ccd4_evt.fits' ))[1].data
ccd5 = fits.open(os.path.join( data_dir, 't0erass_ccd5_evt.fits' ))[1].data
ccd6 = fits.open(os.path.join( data_dir, 't0erass_ccd6_evt.fits' ))[1].data
ccd7 = fits.open(os.path.join( data_dir, 't0erass_ccd7_evt.fits' ))[1].data
simput = Table.read( os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'Xgas_bHS0.8_simput.fits') )


fig_out = os.path.join(os.environ['UCHUU'], 'ra-dec-events-simput.png' )
plt.figure(2, (5.5, 5.5))
N_evt = len(ccd1['DEC']) + len(ccd2['DEC']) + len(ccd3['DEC']) +  len(ccd4['DEC']) + len(ccd5['DEC']) + len(ccd6['DEC']) + len(ccd7['DEC'])
plt.plot(ccd1['RA'], ccd1['DEC'], 'ko', markersize=3, alpha=0.5)
plt.plot(ccd2['RA'], ccd2['DEC'], 'ko', markersize=3, alpha=0.5)
plt.plot(ccd3['RA'], ccd3['DEC'], 'ko', markersize=3, alpha=0.5)
plt.plot(ccd4['RA'], ccd4['DEC'], 'ko', markersize=3, alpha=0.5)
plt.plot(ccd5['RA'], ccd5['DEC'], 'ko', markersize=3, alpha=0.5)
plt.plot(ccd6['RA'], ccd6['DEC'], 'ko', markersize=3, alpha=0.5)
plt.plot(ccd7['RA'], ccd7['DEC'], 'ko', markersize=3, alpha=0.5, label=str(N_evt)+' events')
N_max = int(np.log10(np.max(simput['FLUX'])))
N_min = int(np.log10(np.min(simput['FLUX'])))
FX_bds = np.arange(N_min-1, N_max+1,1)
mks_all = (1+np.arange(len(FX_bds)))*4
for fx_lo, mks in zip(FX_bds, mks_all):
    f_cut = (simput['FLUX']>10**(float(fx_lo))) & (simput['FLUX']<10**(float(fx_lo+1)))
    if len(simput['RA'][f_cut])>0:
        plt.plot(simput['RA'][f_cut], simput['DEC'][f_cut], marker='o', fillstyle='none', markersize=mks, ls='none',
             label='DM haloes '+str(int(fx_lo)) +'<FX<'+str(int(fx_lo+1)) )

plt.xlabel(r'R.A.')
plt.ylabel(r'Dec.')
plt.legend(loc=2, bbox_to_anchor=(-0.1, 1.15), fontsize=10, ncol=2)
plt.tight_layout()
plt.savefig(fig_out)
plt.clf()
print(fig_out, 'written')
