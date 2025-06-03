# !/usr/bin/env python
import sys, os, glob
import numpy as n
from astropy.table import Table, vstack
import astropy.io.fits as fits
sky_map_hdu = Table.read(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS.fits') )
import numpy as np
GE_name = sys.argv[1]

#sky_map_hdu['N_BG_files']  = np.zeros(len(sky_map_hdu))
#sky_map_hdu['N_BG_files2']  = np.zeros(len(sky_map_hdu))
sky_map_hdu['has_merged_events']  = np.ones(len(sky_map_hdu))==0
sky_map_hdu['has_EvtImg']         = np.ones(len(sky_map_hdu))==0
sky_map_hdu['has_ExpMap']         = np.ones(len(sky_map_hdu))==0
sky_map_hdu['has_DetMsk']         = np.ones(len(sky_map_hdu))==0
sky_map_hdu['has_Bo1Cat']         = np.ones(len(sky_map_hdu))==0
sky_map_hdu['has_Bg1Map']         = np.ones(len(sky_map_hdu))==0
sky_map_hdu['has_Bo2Cat']         = np.ones(len(sky_map_hdu))==0
sky_map_hdu['has_Bg2Map']         = np.ones(len(sky_map_hdu))==0
sky_map_hdu['has_Bo3Cat']         = np.ones(len(sky_map_hdu))==0
sky_map_hdu['has_Bg3Map']         = np.ones(len(sky_map_hdu))==0
sky_map_hdu['has_CheMsk']         = np.ones(len(sky_map_hdu))==0
sky_map_hdu['has_ML1Cat']         = np.ones(len(sky_map_hdu))==0
sky_map_hdu['has_Src1Map']        = np.ones(len(sky_map_hdu))==0
sky_map_hdu['has_Sc1Cat']         = np.ones(len(sky_map_hdu))==0
sky_map_hdu['has_PSFMap_eef0.75'] = np.ones(len(sky_map_hdu))==0
sky_map_hdu['has_ApeSen_eef0.75'] = np.ones(len(sky_map_hdu))==0
sky_map_hdu['has_ApeCat_eef0.75'] = np.ones(len(sky_map_hdu))==0
sky_map_hdu['has_srcAUTOreg']     = np.ones(len(sky_map_hdu))==0
sky_map_hdu['has_bkgAUTOreg']     = np.ones(len(sky_map_hdu))==0
MinTotalCts = 2
sky_map_hdu['has_IDMatch_Uniq_Tot'+ str(MinTotalCts)] = np.ones(len(sky_map_hdu))==0

LC_dir = "LCerass"

for jj, sky_tile_value in enumerate(sky_map_hdu['SRVMAP']):

    sky_tile_id = str(sky_tile_value)
    str_field = sky_tile_id.zfill(6)

    esass_dir = os.path.join("/home/idies/workspace/erosim/Uchuu/", LC_dir, str_field, GE_name, 'eSASS')
    GE_dir = os.path.join("/home/idies/workspace/erosim/Uchuu/", LC_dir, str_field, GE_name)
    VerBand = str(4)
    outprefix = str_field + "_"  # ""
    #bg_dir      = os.path.join( os.environ['UCHUU'], LC_dir, str_field, 'pBG' ) # 'evt_particle_???.fits' )
    #BG_evt_files = n.array( glob.glob( os.path.join( bg_dir, 'evt_particle_???.fits' ) ) )
    #sky_map_hdu['N_BG_files'][jj] = len(BG_evt_files)
    #bg_dir2      = os.path.join( os.environ['UCHUU'], LC_dir, str_field, 'pBG2' ) # 'evt_particle_???.fits' )
    #BG_evt_files = n.array( glob.glob( os.path.join( bg_dir2, '*.fits' ) ) )
    #sky_map_hdu['N_BG_files'][jj] = len(BG_evt_files)
    sky_map_hdu['has_merged_events']  [jj] = os.path.isfile(os.path.join("/home/idies/workspace/erosim/Uchuu/LCerass/", str_field, GE_name, 'evt_'+str_field+'.fits'))
    sky_map_hdu['has_EvtImg']         [jj] = os.path.isfile(os.path.join(esass_dir, f"{outprefix}02{VerBand}_EvtImg.fits"))
    sky_map_hdu['has_ExpMap']         [jj] = os.path.isfile(os.path.join(esass_dir, f"{outprefix}02{VerBand}_ExpMap.fits"))
    sky_map_hdu['has_DetMsk']         [jj] = os.path.isfile(os.path.join(esass_dir, f"{outprefix}02{VerBand}_DetMsk.fits"))
    sky_map_hdu['has_Bo1Cat']         [jj] = os.path.isfile(os.path.join(esass_dir, f"{outprefix}02{VerBand}_Bo1Cat.fits"))
    sky_map_hdu['has_Bg1Map']         [jj] = os.path.isfile(os.path.join(esass_dir, f"{outprefix}02{VerBand}_Bg1Map.fits"))
    sky_map_hdu['has_Bo2Cat']         [jj] = os.path.isfile(os.path.join(esass_dir, f"{outprefix}02{VerBand}_Bo2Cat.fits"))
    sky_map_hdu['has_Bg2Map']         [jj] = os.path.isfile(os.path.join(esass_dir, f"{outprefix}02{VerBand}_Bg2Map.fits"))
    sky_map_hdu['has_Bo3Cat']         [jj] = os.path.isfile(os.path.join(esass_dir, f"{outprefix}02{VerBand}_Bo3Cat.fits"))
    sky_map_hdu['has_Bg3Map']         [jj] = os.path.isfile(os.path.join(esass_dir, f"{outprefix}02{VerBand}_Bg3Map.fits"))
    sky_map_hdu['has_CheMsk']         [jj] = os.path.isfile(os.path.join(esass_dir, f"{outprefix}02{VerBand}_CheMsk.fits"))
    sky_map_hdu['has_ML1Cat']         [jj] = os.path.isfile(os.path.join(esass_dir, f"{outprefix}02{VerBand}_ML1Cat.fits"))
    sky_map_hdu['has_Src1Map']        [jj] = os.path.isfile(os.path.join(esass_dir, f"{outprefix}02{VerBand}_Src1Map.fits"))
    sky_map_hdu['has_Sc1Cat']         [jj] = os.path.isfile(os.path.join(esass_dir, f"{outprefix}02{VerBand}_Sc1Cat.fits"))
    sky_map_hdu['has_PSFMap_eef0.75'] [jj] = os.path.isfile(os.path.join(esass_dir, f"{outprefix}02{VerBand}_PSFMap_eef0.75.fits"))
    sky_map_hdu['has_ApeSen_eef0.75'] [jj] = os.path.isfile(os.path.join(esass_dir, f"{outprefix}02{VerBand}_ApeSen_eef0.75.fits"))
    sky_map_hdu['has_ApeCat_eef0.75'] [jj] = os.path.isfile(os.path.join(esass_dir, f"{outprefix}02{VerBand}_ApeCat_eef0.75.fits"))
    sky_map_hdu['has_srcAUTOreg']     [jj] = os.path.isfile(os.path.join(esass_dir, f"{outprefix}_02_23_srcAUTO.reg"))
    sky_map_hdu['has_bkgAUTOreg']     [jj] = os.path.isfile(os.path.join(esass_dir, f"{outprefix}_02_23_bkgAUTO.reg"))
    sky_map_hdu['has_IDMatch_Uniq_Tot' + str(MinTotalCts)][jj] = os.path.isfile(os.path.join(GE_dir, f'Sc1_' + str_field + '_IDMatch_Uniq_Tot' + str(MinTotalCts) + '.fits'))

sky_map_hdu.write(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS_'+GE_name+'.fits'), overwrite = True)
