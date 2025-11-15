"""
This script resambles https://gitlab.mpcdf.mpg.de/joco/erosita_sxrbg/-/blob/main/esass/erassX_write_scripts.py
It creates 3 scripts:
1) create images
2) start detection: 3 loops of erbox and erbackmap
3) finish detection: ermldet, apetool, srctool
/data26s/mpecl/eRASS1/??????/c946
/data26s/mpecl/eRASS1/358144/c946
/data26s/mpecl/eRASS1/358144/c946/*events*

nohup bash /data26s/comparat/simulations/erosim/eRASS1_UNIT_fA1i_2021_10_12_SKYMAP/086066/086066_pipeline_img.sh > /data26s/comparat/simulations/erosim/eRASS1_UNIT_fA1i_2021_10_12_SKYMAP/086066/086066_pipeline_img_RS.log &
nohup bash /data26s/comparat/simulations/erosim/eRASS1_UNIT_fA1i_2021_10_12_SKYMAP/086066/086066_pipeline_det1.sh > /data26s/comparat/simulations/erosim/eRASS1_UNIT_fA1i_2021_10_12_SKYMAP/086066/086066_pipeline_det1_RS.log &
nohup bash /data26s/comparat/simulations/erosim/eRASS1_UNIT_fA1i_2021_10_12_SKYMAP/086066/086066_pipeline_det2.sh > /data26s/comparat/simulations/erosim/eRASS1_UNIT_fA1i_2021_10_12_SKYMAP/086066/086066_pipeline_det2_RS.log &

nohup bash /data26s/comparat/simulations/erosim/eRASS1_UNIT_fA1i_2021_10_12_SKYMAP/092084/092084_pipeline_img.sh > /data26s/comparat/simulations/erosim/eRASS1_UNIT_fA1i_2021_10_12_SKYMAP/092084/092084_pipeline_img_RS.log &
nohup bash /data26s/comparat/simulations/erosim/eRASS1_UNIT_fA1i_2021_10_12_SKYMAP/092084/092084_pipeline_det1.sh > /data26s/comparat/simulations/erosim/eRASS1_UNIT_fA1i_2021_10_12_SKYMAP/092084/092084_pipeline_det1_RS.log &
nohup bash /data26s/comparat/simulations/erosim/eRASS1_UNIT_fA1i_2021_10_12_SKYMAP/092084/092084_pipeline_det2.sh > /data26s/comparat/simulations/erosim/eRASS1_UNIT_fA1i_2021_10_12_SKYMAP/092084/092084_pipeline_det2_RS.log &
"""
# !/usr/bin/env python
import sys, os, glob
from astropy.table import Table, vstack
import astropy.io.fits as fits
sky_map_hdu = Table.read(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS.fits') )
import numpy as np

sky_map_hdu['N_files_BG']  = np.zeros(len(sky_map_hdu))
sky_map_hdu['N_files_Att_eRASS8_sixte_v27_SEED_001_events_cluster']  = np.zeros(len(sky_map_hdu))
sky_map_hdu['N_files_Att_eRASS8_sixte_v27_SEED_002_events_cluster']  = np.zeros(len(sky_map_hdu))
sky_map_hdu['N_files_Att_eRASS8_sixte_v27_SEED_003_events_cluster']  = np.zeros(len(sky_map_hdu))
sky_map_hdu['N_files_Att_eRASS8_sixte_v27_SEED_004_events_cluster']  = np.zeros(len(sky_map_hdu))
sky_map_hdu['N_files_Att_eRASS8_sixte_v27_SEED_005_events_cluster']  = np.zeros(len(sky_map_hdu))
sky_map_hdu['N_files_Att_eRASS8_sixte_v27_SEED_006_events_cluster']  = np.zeros(len(sky_map_hdu))
sky_map_hdu['N_files_Att_eRASS8_sixte_v27_SEED_007_events_cluster']  = np.zeros(len(sky_map_hdu))
sky_map_hdu['N_files_Att_eRASS8_sixte_v27_SEED_008_events_cluster']  = np.zeros(len(sky_map_hdu))
sky_map_hdu['N_files_Att_eRASS8_sixte_v27_SEED_009_events_cluster']  = np.zeros(len(sky_map_hdu))
sky_map_hdu['N_files_eRASS8_SEED_001_events_AGN_2025_04'] = np.zeros(len(sky_map_hdu))
sky_map_hdu['N_files_eRASS8_SEED_002_events_AGN_2025_04'] = np.zeros(len(sky_map_hdu))
sky_map_hdu['N_files_eRASS8_SEED_003_events_AGN_2025_04'] = np.zeros(len(sky_map_hdu))
sky_map_hdu['N_files_eRASS8_SEED_004_events_AGN_2025_04'] = np.zeros(len(sky_map_hdu))
sky_map_hdu['N_files_eRASS8_SEED_005_events_AGN_2025_04'] = np.zeros(len(sky_map_hdu))
sky_map_hdu['N_files_eRASS8_SEED_006_events_AGN_2025_04'] = np.zeros(len(sky_map_hdu))
sky_map_hdu['N_files_eRASS8_SEED_007_events_AGN_2025_04'] = np.zeros(len(sky_map_hdu))
sky_map_hdu['N_files_eRASS8_SEED_008_events_AGN_2025_04'] = np.zeros(len(sky_map_hdu))
sky_map_hdu['N_files_eRASS8_SEED_009_events_AGN_2025_04'] = np.zeros(len(sky_map_hdu))
sky_map_hdu['N_files_stars'] = np.zeros(len(sky_map_hdu))
LC_dir = "LCerass"

for jj, sky_tile_value in enumerate(sky_map_hdu['SRVMAP']):
    sky_tile_id = str(sky_tile_value)
    str_field = sky_tile_id.zfill(6)
    sky_map_hdu['N_files_BG'][jj] = len(n.array( glob.glob( os.path.join( os.environ['UCHUU'], LC_dir, str_field, 'pBG2', '*.fits' ) ) ))
    sky_map_hdu['N_files_Att_eRASS8_sixte_v27_SEED_001_events_cluster'] = len(np.array( glob.glob( os.path.join( "/home/idies/workspace/erosim/Uchuu/", LC_dir, str_field, 'Att_eRASS8_sixte_v27_SEED_001_events_cluster_Xgas_bHS0.8', '*.fits' ) ) ))
    sky_map_hdu['N_files_Att_eRASS8_sixte_v27_SEED_002_events_cluster'] = len(np.array( glob.glob( os.path.join( "/home/idies/workspace/erosim/Uchuu/", LC_dir, str_field, 'Att_eRASS8_sixte_v27_SEED_002_events_cluster_Xgas_bHS0.8', '*.fits' ) ) ))
    sky_map_hdu['N_files_Att_eRASS8_sixte_v27_SEED_003_events_cluster'] = len(np.array( glob.glob( os.path.join( "/home/idies/workspace/erosim/Uchuu/", LC_dir, str_field, 'Att_eRASS8_sixte_v27_SEED_003_events_cluster_Xgas_bHS0.8', '*.fits' ) ) ))
    sky_map_hdu['N_files_Att_eRASS8_sixte_v27_SEED_004_events_cluster'] = len(np.array( glob.glob( os.path.join( "/home/idies/workspace/erosim/Uchuu/", LC_dir, str_field, 'Att_eRASS8_sixte_v27_SEED_004_events_cluster_Xgas_bHS0.8', '*.fits' ) ) ))
    sky_map_hdu['N_files_Att_eRASS8_sixte_v27_SEED_005_events_cluster'] = len(np.array( glob.glob( os.path.join( "/home/idies/workspace/erosim/Uchuu/", LC_dir, str_field, 'Att_eRASS8_sixte_v27_SEED_005_events_cluster_Xgas_bHS0.8', '*.fits' ) ) ))
    sky_map_hdu['N_files_Att_eRASS8_sixte_v27_SEED_006_events_cluster'] = len(np.array( glob.glob( os.path.join( "/home/idies/workspace/erosim/Uchuu/", LC_dir, str_field, 'Att_eRASS8_sixte_v27_SEED_006_events_cluster_Xgas_bHS0.8', '*.fits' ) ) ))
    sky_map_hdu['N_files_Att_eRASS8_sixte_v27_SEED_007_events_cluster'] = len(np.array( glob.glob( os.path.join( "/home/idies/workspace/erosim/Uchuu/", LC_dir, str_field, 'Att_eRASS8_sixte_v27_SEED_007_events_cluster_Xgas_bHS0.8', '*.fits' ) ) ))
    sky_map_hdu['N_files_Att_eRASS8_sixte_v27_SEED_008_events_cluster'] = len(np.array( glob.glob( os.path.join( "/home/idies/workspace/erosim/Uchuu/", LC_dir, str_field, 'Att_eRASS8_sixte_v27_SEED_008_events_cluster_Xgas_bHS0.8', '*.fits' ) ) ))
    sky_map_hdu['N_files_Att_eRASS8_sixte_v27_SEED_009_events_cluster'] = len(np.array( glob.glob( os.path.join( "/home/idies/workspace/erosim/Uchuu/", LC_dir, str_field, 'Att_eRASS8_sixte_v27_SEED_009_events_cluster_Xgas_bHS0.8', '*.fits' ) ) ))
    sky_map_hdu['N_files_eRASS8_SEED_001_events_AGN_2025_04'] = len(np.array( glob.glob( os.path.join( "/home/idies/workspace/erosim/Uchuu/", LC_dir, str_field, 'eRASS8_SEED_001_events_AGN_2025_04', '*.fits' ) ) ))
    sky_map_hdu['N_files_eRASS8_SEED_002_events_AGN_2025_04'] = len(np.array( glob.glob( os.path.join( "/home/idies/workspace/erosim/Uchuu/", LC_dir, str_field, 'eRASS8_SEED_002_events_AGN_2025_04', '*.fits' ) ) ))
    sky_map_hdu['N_files_eRASS8_SEED_003_events_AGN_2025_04'] = len(np.array( glob.glob( os.path.join( "/home/idies/workspace/erosim/Uchuu/", LC_dir, str_field, 'eRASS8_SEED_003_events_AGN_2025_04', '*.fits' ) ) ))
    sky_map_hdu['N_files_eRASS8_SEED_004_events_AGN_2025_04'] = len(np.array( glob.glob( os.path.join( "/home/idies/workspace/erosim/Uchuu/", LC_dir, str_field, 'eRASS8_SEED_004_events_AGN_2025_04', '*.fits' ) ) ))
    sky_map_hdu['N_files_eRASS8_SEED_005_events_AGN_2025_04'] = len(np.array( glob.glob( os.path.join( "/home/idies/workspace/erosim/Uchuu/", LC_dir, str_field, 'eRASS8_SEED_005_events_AGN_2025_04', '*.fits' ) ) ))
    sky_map_hdu['N_files_eRASS8_SEED_006_events_AGN_2025_04'] = len(np.array( glob.glob( os.path.join( "/home/idies/workspace/erosim/Uchuu/", LC_dir, str_field, 'eRASS8_SEED_006_events_AGN_2025_04', '*.fits' ) ) ))
    sky_map_hdu['N_files_eRASS8_SEED_007_events_AGN_2025_04'] = len(np.array( glob.glob( os.path.join( "/home/idies/workspace/erosim/Uchuu/", LC_dir, str_field, 'eRASS8_SEED_007_events_AGN_2025_04', '*.fits' ) ) ))
    sky_map_hdu['N_files_eRASS8_SEED_008_events_AGN_2025_04'] = len(np.array( glob.glob( os.path.join( "/home/idies/workspace/erosim/Uchuu/", LC_dir, str_field, 'eRASS8_SEED_008_events_AGN_2025_04', '*.fits' ) ) ))
    sky_map_hdu['N_files_eRASS8_SEED_009_events_AGN_2025_04'] = len(np.array( glob.glob( os.path.join( "/home/idies/workspace/erosim/Uchuu/", LC_dir, str_field, 'eRASS8_SEED_009_events_AGN_2025_04', '*.fits' ) ) ))
    sky_map_hdu['N_files_stars'] = len(n.array( glob.glob( os.path.join( "/home/idies/workspace/erosim/Uchuu/", LC_dir, str_field, 'stars', '*.fits' ) ) ))
    print(sky_tile_value)
p_2_out = os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS_event_simulated.fits')
sky_map_hdu.write(p_2_out, overwrite = True)
print(p_2_out)
