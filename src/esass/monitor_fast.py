# !/usr/bin/env python
import os
import numpy as np
import astropy.table as tbl
from multiprocessing import Pool
from functools import partial

#Function for processing single tile
def make_summary_skymap_fast(sky_tile_el, other_elements):

    #Separate function elements
    jj, sky_tile_value = sky_tile_el
    GE_name = other_elements[0]

    #Get skytile id
    sky_tile_id = str(sky_tile_value)
    str_field = sky_tile_id.zfill(6)

    #Set up paths
    LC_dir = "LCerass"
    esass_dir = os.path.join("/home/idies/workspace/erosim/Uchuu/", LC_dir, str_field, GE_name, 'eSASS')
    GE_dir = os.path.join("/home/idies/workspace/erosim/Uchuu/", LC_dir, str_field, GE_name)
    VerBand = str(4)
    outprefix = str_field + "_"

    #Create array of tile status
    tile_status = [jj,
        os.path.isfile(os.path.join("/home/idies/workspace/erosim/Uchuu/LCerass/", str_field, GE_name, 'evt_'+str_field+'.fits')),
        os.path.isfile(os.path.join(esass_dir, f"{outprefix}02{VerBand}_EvtImg.fits")),
        os.path.isfile(os.path.join(esass_dir, f"{outprefix}02{VerBand}_ExpMap.fits")),
        os.path.isfile(os.path.join(esass_dir, f"{outprefix}02{VerBand}_DetMsk.fits")),
        os.path.isfile(os.path.join(esass_dir, f"{outprefix}02{VerBand}_Bo1Cat.fits")),
        os.path.isfile(os.path.join(esass_dir, f"{outprefix}02{VerBand}_Bg1Map.fits")),
        os.path.isfile(os.path.join(esass_dir, f"{outprefix}02{VerBand}_Bo2Cat.fits")),
        os.path.isfile(os.path.join(esass_dir, f"{outprefix}02{VerBand}_Bg2Map.fits")),
        os.path.isfile(os.path.join(esass_dir, f"{outprefix}02{VerBand}_Bo3Cat.fits")),
        os.path.isfile(os.path.join(esass_dir, f"{outprefix}02{VerBand}_Bg3Map.fits")),
        os.path.isfile(os.path.join(esass_dir, f"{outprefix}02{VerBand}_CheMsk.fits")),
        os.path.isfile(os.path.join(esass_dir, f"{outprefix}02{VerBand}_ML1Cat.fits")),
        os.path.isfile(os.path.join(esass_dir, f"{outprefix}02{VerBand}_Src1Map.fits")),
        os.path.isfile(os.path.join(esass_dir, f"{outprefix}02{VerBand}_Sc1Cat.fits")),
        os.path.isfile(os.path.join(esass_dir, f"{outprefix}02{VerBand}_PSFMap_eef0.75.fits")),
        os.path.isfile(os.path.join(esass_dir, f"{outprefix}02{VerBand}_ApeSen_eef0.75.fits")),
        os.path.isfile(os.path.join(esass_dir, f"{outprefix}02{VerBand}_ApeCat_eef0.75.fits")),
        os.path.isfile(os.path.join(esass_dir, f"{outprefix}_02_23_srcAUTO.reg")),
        os.path.isfile(os.path.join(esass_dir, f"{outprefix}_02_23_bkgAUTO.reg")),
        os.path.isfile(os.path.join(GE_dir, f'Sc1_' + str_field + '_IDMatch_Uniq_Tot' + str(MinTotalCts) + '.fits'))
        ]

    return tile_status

# =============================================================================

#Define environmental variables
os.environ['UCHUU']='/home/idies/workspace/erosim/Uchuu'
os.environ['GIT_STMOD']='/home/idies/workspace/erosim/software/st_mod'
os.environ['GIT_STMOD_DATA']='/home/idies/workspace/erosim/software/st_mod_data'
debug == True

#Define seed list
clu_seed = list(range(19,118))
agn_seed = list(range(1,10))*11
t_tot_seed = list(zip(agn_seed, clu_seed))
tot_seed = [(str(seed[0]).zfill(3), str(seed[1]).zfill(3)) for seed in t_tot_seed]
exp_name = ['GE_e4_merge_AGNseed{0}_SimBKG_CLUseed{1}'.format(seed[0], seed[1]) for seed in tot_seed] + ['GE_e5_merge_AGNseed{0}_SimBKG_CLUseed{1}'.format(seed[0], seed[1]) for seed in tot_seed]

#Read skymap
sky_map_hdu_general = tbl.Table.read(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS.fits'))
MinTotalCts = 2
colname_list = ['has_merged_events', 'has_EvtImg', 'has_ExpMap', 'has_DetMsk', 'has_Bo1Cat', 'has_Bg1Map', 'has_Bo2Cat', 'has_Bg2Map', 'has_Bo3Cat', 'has_Bg3Map', 'has_CheMsk', 'has_ML1Cat', 'has_Src1Map', 'has_Sc1Cat', 'has_PSFMap_eef0.75', 'has_ApeSen_eef0.75', 'has_ApeCat_eef0.75', 'has_srcAUTOreg', 'has_bkgAUTOreg', 'has_IDMatch_Uniq_Tot'+ str(MinTotalCts)]

#Cycle experiment names
for GE_name in exp_name:

    if debug:
        print(GE_name)

    #Create copy of skymap
    sky_map_hdu = tbl.Table(sky_map_hdu_general, copy = True)

    #Create needed columns
    for cn in colname_list:
        sky_map_hdu[cn] = np.ones(len(sky_map_hdu))==0
        
    #Define function for multiprocessing
    onepool_func = partial(make_summary_skymap_fast, other_elements = [GE_name])

    #Map to cores    
    with Pool(20) as p:
        all_tiles_array = p.map(onepool_func, list(enumerate(sky_map_hdu['SRVMAP'])))

    for row in all_tiles_array:
        tile_idx = row[0]
        one_tile_status = row[1:]
        for vali, val in enumerate(one_tile_status):
            sky_map_hdu[colname_list[vali]][tile_idx] = val

    #Save to file
    sky_map_hdu.write(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS_'+GE_name+'.fits'), overwrite = True)

    #Save ram
    sky_map_hdu = 0
    all_tiles_array = 0