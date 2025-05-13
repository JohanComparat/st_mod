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
import numpy as n
from astropy.table import Table, vstack
import astropy.io.fits as fits
sky_map_hdu = Table.read(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS.fits') )

GE_names = ['GE_e4_merge_AGNseed001_SimBKG', 'GE_e4_merge_AGNseed001_SimBKG_CLUseed001', 'GE_e4_merge_SimBKG_CLUseed001', 'GE_e4_merge_SimBKG']

for sky_tile in sky_map_hdu[(sky_map_hdu['OWNER']==2)|(sky_map_hdu['OWNER']==0)][2:48]:
    sky_tile_id = str(sky_tile['SRVMAP'])
    str_field = str(sky_tile['SRVMAP']).zfill(6)

    for GE_name in GE_names:
        indir = os.path.join("/home/idies/workspace/erosim/Uchuu/LCerass/", str_field, GE_name)
        esass_dir = os.path.join(indir, 'eSASS')
        git_dir = os.path.join(os.environ['GIT_STMOD'], 'src/esass' )
        print ("cd "+esass_dir )
        print ("sh "+str_field+"_pipeline_img1.sh")
        print ("sh "+str_field+"_pipeline_det1.sh")
        print ("sh "+str_field+"_pipeline_Src1.sh")
        print ("cd "+git_dir)
        if not GE_name=='GE_e4_merge_SimBKG':
            print ("python photon_matching_RS.py "+GE_name+" "+str_field )
        print('# =====')



