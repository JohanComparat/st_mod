#import subprocess
import os
#import errno
import sys, glob
#import healpy
import numpy as np
from astropy.io import fits
from astropy.table import Table, Column
sky_map_hdu = Table.read(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS.fits'))
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as plt
LC_dir = 'LCerass'

for srv_val in sky_map_hdu['SRVMAP'][(sky_map_hdu['OWNER']==2)|(sky_map_hdu['OWNER']==0)]:
    str_field = str(srv_val).zfill(6)
    all_glist = np.array(
        glob.glob(os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'z?p??', 'replication_*', "glist.fits")))
    print(str_field, len(all_glist))

sky_tile_num = 121048 # int(sys.argv[1])
seed = 1
LC_dir = 'LCerass'
erass_option = "eRASS8"
print(seed, LC_dir, erass_option, sky_tile_num)#, env, erass_option)

sky_tile = sky_map_hdu[sky_map_hdu['SRVMAP']==sky_tile_num]
sky_tile_id = str(sky_tile_num)
str_field = sky_tile_id.zfill(6)
ra_cen = sky_tile['RA_CEN']
dec_cen = sky_tile['DE_CEN']
gal = Table.read(os.path.join(os.environ['UCHUU'], LC_dir, str_field, "glist.fits" ))
xagn = Table.read(os.path.join(os.environ['UCHUU'], LC_dir, str_field, "AGN_list_sigma_0.8_fsat_8.0.fits" ))

all_glist = np.array( glob.glob( os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'z?p??', 'replication_*', "glist.fits" ) ) )
print('N_files', len(all_glist))


figure_dir = os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'figures')
os.system('mkdir -p '+figure_dir )
print('figures here', figure_dir)

fig_out = os.path.join(figure_dir, 'glist-ra-dec.png' )
plt.figure(2, (5.5, 5.5))
plt.plot(gal['RA'], gal['DEC'], 'k+', markersize=0.2, alpha=0.5, rasterized=True)
plt.xlabel(r'R.A.')
plt.ylabel(r'Dec.')
plt.tight_layout()
plt.savefig(fig_out)
plt.clf()
print(fig_out, 'written')


fig_out = os.path.join(figure_dir, 'glist-z-hist.png' )
plt.figure(2, (5.5, 5.5))
plt.hist(gal['redshift_S'], bins=np.arange(0, 4.1, 0.05))
plt.hist(xagn['redshift_S'], bins=np.arange(0, 4.1, 0.05), histtype='step')
plt.xlabel(r'R.A.')
plt.ylabel(r'Dec.')
plt.yscale('log')
plt.tight_layout()
plt.savefig(fig_out)
plt.clf()
print(fig_out, 'written')

