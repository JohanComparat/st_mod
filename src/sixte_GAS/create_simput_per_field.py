import numpy as np
import numpy
import h5py
import healpy
import os, sys, glob

from sixte_simput_lib import make_agn_simput

from scipy.interpolate import interp1d

import astropy.io.fits as fits
from astropy.table import Table, Column, vstack, hstack
import astropy.units as u
from astropy.coordinates import SkyCoord
import time
#t0 = time.time()

nl = lambda selection : len(selection.nonzero()[0])
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoUNIT = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
cosmo2 = FlatLambdaCDM(H0=68 * u.km / u.s / u.Mpc, Om0=0.31)
h = 0.6774
L_box = 1000.0 / h
cosmo = cosmoUNIT

data_s5 = '/data56s/comparat/erosim/data_s5'
data_s5 = os.environ['DATA_S5']

sky_tile_id = '015114'

# output directories
dir_sky_tile_uchuu = os.path.join(data_s5, sky_tile_id,'UCHUU')
dir_simput_gas     = os.path.join(dir_sky_tile_uchuu, 'simput_gas'          )
dir_simput_gasI    = os.path.join(dir_sky_tile_uchuu, 'simput_gas', 'images')
dir_simput_gasS    = os.path.join(dir_sky_tile_uchuu, 'simput_gas', 'gas_Xspectra')
dir_simput_agn     = os.path.join(dir_sky_tile_uchuu, 'simput_agn'          )
# make soft link for AGN spectra
dir_simput_agnS     = os.path.join(dir_sky_tile_uchuu, 'simput_agn', 'agn_Xspectra' )
os.system('ln -s /afs/mpe/www/people/comparat/eROSITA_AGN_mock/spectra/Xray/spectra '+dir_simput_agnS)
dir_simput_xrb     = os.path.join(dir_sky_tile_uchuu, 'simput_xrb'          )
dir_simput_xrbS     = os.path.join(dir_sky_tile_uchuu, 'simput_xrb', 'xrb_Xspectra' )
dir_simput_bkg     = os.path.join(dir_sky_tile_uchuu, 'simput_bkg')
dir_simput_bkgI    = os.path.join(dir_sky_tile_uchuu, 'simput_bkg', 'images')
dir_simput_bkgI    = os.path.join(dir_sky_tile_uchuu, 'simput_bkg', 'bkg_Xspectra')

os.system('mkdir -p ' + dir_sky_tile_uchuu )
os.system('mkdir -p ' + dir_simput_gas     )
os.system('mkdir -p ' + dir_simput_gasI    )
os.system('mkdir -p ' + dir_simput_agn     )
os.system('mkdir -p ' + dir_simput_xrb     )
os.system('mkdir -p ' + dir_simput_bkg     )
os.system('mkdir -p ' + dir_simput_bkgI    )

skymap = Table.read( os.path.join( os.environ['GIT_ER4W'], 'data' , 'SKYMAPS.fits') )
field_selection = ( skymap['SRVMAP']==int(sky_tile_id) )
sti = skymap[field_selection]

# boundary condition from real data
obs_events = Table.read( os.path.join(data_s5, sky_tile_id, 'c030', 'sb05_015114_020_Image_c030.fits') )
evt_ra_min = np.min(obs_events['RA'])
evt_ra_max = np.max(obs_events['RA'])
evt_de_min = np.min(obs_events['DEC'])
evt_de_max = np.max(obs_events['DEC'])

# take simulated file set
uchuu_z_str = 'z0p09'
uchuu_metadata = Table.read( os.path.join( os.environ['UCHUU'], 'area_per_replica_'+uchuu_z_str+'.fits' ) )
# choose the replicas
top_dir = os.path.join( os.environ['UCHUU'], 'FullSky', uchuu_z_str,'replication_*' )
all_GAS_files = np.array( glob.glob( os.path.join( top_dir, 'Xgas.fits' ) ) )
all_GAL_files = np.array( glob.glob( os.path.join( top_dir, 'glist.fits' ) ) )
all_AGN_files = np.array( glob.glob( os.path.join( top_dir, 'AGN_list_sigma_0.8_fsat_8.0.fits' ) ) )

t_Ref = time.time()
N_in_replica = np.zeros(len(all_GAL_files))
for jj, p_2_GAL in enumerate(all_GAL_files):
    t0 = Table.read(p_2_GAL)
    area_selection = (t0['RA']>=evt_ra_min)&(t0['RA']<=evt_ra_max)&(t0['DEC']>=evt_de_min)&(t0['DEC']<=evt_de_max)
    N_in_replica[jj] = nl(area_selection)
    print(jj, N_in_replica[jj], time.time()-t_Ref)

replicas_selected = (N_in_replica>1)

# retrieve all files and simulated objects in that field (for that redshift slice)
# tabulated pathes to outputs
GAL, AGN, GAS = [], [], []
p_2_simput_GAL, p_2_simput_AGN, p_2_simput_GAS = [], [], []
for jj, (p_2_GAL, p_2_GAS, p_2_AGN) in enumerate(zip(all_GAL_files[replicas_selected], all_GAS_files[replicas_selected], all_AGN_files[replicas_selected])):
    GAL_j = Table.read(p_2_GAL)
    GAL_j['line_ID'] = np.arange(len(GAL_j))
    area_selection = (GAL_j['RA']>=evt_ra_min)&(GAL_j['RA']<=evt_ra_max)&(GAL_j['DEC']>=evt_de_min)&(GAL_j['DEC']<=evt_de_max)
    GAS_j = Table.read(p_2_GAS)
    AGN_j = Table.read(p_2_AGN)
    GAL_AGN_j = GAL_j[AGN_j['ID_glist']]
    GAL_GAS_j = GAL_j[GAS_j['ID_glist']]
    GAL_AGN_j.keep_columns(['line_ID', 'RA', 'DEC', 'nH'])
    GAL_GAS_j.keep_columns(['line_ID', 'RA', 'DEC', 'nH' ])
    GAL_j = GAL_j[area_selection]
    sel_GAS_j = np.isin(GAS_j['ID_glist'], GAL_j['line_ID'])
    sel_AGN_j = np.isin(AGN_j['ID_glist'], GAL_j['line_ID'])
    GAS_j = GAS_j[ sel_GAS_j ]
    AGN_j = AGN_j[ sel_AGN_j ]
    GAL_GAS_j = GAL_GAS_j[ sel_GAS_j ]
    GAL_AGN_j = GAL_AGN_j[ sel_AGN_j ]
    GAL.append(GAL_j)
    AGN.append( hstack((AGN_j, GAL_AGN_j)) )
    GAS.append( hstack((GAS_j, GAL_GAS_j)) )
    p_2_simput_GAL.append( os.path.join(dir_simput_xrb, p_2_GAL.split('/')[-3]+'-'+p_2_GAL.split('/')[-2]+'-'+p_2_GAL.split('/')[-1] ) )
    p_2_simput_AGN.append( os.path.join(dir_simput_agn, p_2_AGN.split('/')[-3]+'-'+p_2_AGN.split('/')[-2]+'-'+p_2_AGN.split('/')[-1] ) )
    p_2_simput_GAS.append( os.path.join(dir_simput_gas, p_2_GAS.split('/')[-3]+'-'+p_2_GAS.split('/')[-2]+'-'+p_2_GAS.split('/')[-1] ) )

# turn AGN into a simput file
make_agn_simput(AGN[0], path_2_SMPT_catalog=p_2_simput_AGN[0], FX_LIM_value_cen = -16)

# turn the GAL into an XRB simput file

# turn GAS into a simput file and link to images
