import sys, os, glob
import numpy as n
import numpy as np
from astropy.table import Table, vstack
import astropy.io.fits as fits

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as plt

sky_map_hdu = Table.read(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS.fits') )

GE_names = ['GE_e4_merge_AGNseed001_SimBKG', 'GE_e4_merge_AGNseed001_SimBKG_CLUseed001', 'GE_e4_merge_SimBKG_CLUseed001',
            'GE_e4_merge_AGNseed002_SimBKG', 'GE_e4_merge_AGNseed002_SimBKG_CLUseed002', 'GE_e4_merge_SimBKG_CLUseed002',
            'GE_e4_merge_AGNseed003_SimBKG', 'GE_e4_merge_AGNseed003_SimBKG_CLUseed003', 'GE_e4_merge_SimBKG_CLUseed003',
            'GE_e4_merge_AGNseed004_SimBKG', 'GE_e4_merge_AGNseed004_SimBKG_CLUseed004', 'GE_e4_merge_SimBKG_CLUseed004',
            'GE_e4_merge_AGNseed005_SimBKG', 'GE_e4_merge_AGNseed005_SimBKG_CLUseed005', 'GE_e4_merge_SimBKG_CLUseed005',
            'GE_e4_merge_AGNseed006_SimBKG', 'GE_e4_merge_AGNseed006_SimBKG_CLUseed006', 'GE_e4_merge_SimBKG_CLUseed006',
            'GE_e4_merge_AGNseed007_SimBKG', 'GE_e4_merge_AGNseed007_SimBKG_CLUseed007', 'GE_e4_merge_SimBKG_CLUseed007',
            'GE_e4_merge_AGNseed008_SimBKG', 'GE_e4_merge_AGNseed008_SimBKG_CLUseed008', 'GE_e4_merge_SimBKG_CLUseed008',
            'GE_e4_merge_SimBKG',
            ]
SKYMAP = {}
for GE_name in GE_names:
    SKYMAP[GE_name] = Table.read(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS_'+GE_name+'.fits'))

datata = []
#for sky_tile in sky_map_hdu[(sky_map_hdu['OWNER']==2)|(sky_map_hdu['OWNER']==0)][2:48]:
for GE_name in GE_names:
    # print('='*100)
    sky_map_hdu = SKYMAP[GE_name]
    to_process = ((sky_map_hdu['OWNER']==2)|(sky_map_hdu['OWNER']==0))&(sky_map_hdu['has_merged_events'])&(sky_map_hdu['has_Sc1Cat']==False)
    already_done = ((sky_map_hdu['OWNER']==2)|(sky_map_hdu['OWNER']==0))&(sky_map_hdu['has_merged_events'])&(sky_map_hdu['has_Sc1Cat'])
    print(len(sky_map_hdu[to_process]), len(sky_map_hdu[already_done]), GE_name)
    print(GE_name, len(sky_map_hdu[to_process]), 'tiles to process', len(sky_map_hdu[already_done]), 'tiles already done')
    datata.append([len(sky_map_hdu[to_process]), len(sky_map_hdu[already_done])])
    p2fig = os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'ra-dec-SKYMAPS_' + GE_name + '.png')
    plt.plot(sky_map_hdu['RA_CEN'][((sky_map_hdu['OWNER']==2)|(sky_map_hdu['OWNER']==0))], sky_map_hdu['DE_CEN'][((sky_map_hdu['OWNER']==2)|(sky_map_hdu['OWNER']==0))], 'k,', label='eRO DE')
    plt.plot(sky_map_hdu['RA_CEN'][already_done], sky_map_hdu['DE_CEN'][already_done], 'k+', label=str(len(sky_map_hdu['DE_CEN'][already_done]))+' done')
    plt.plot(sky_map_hdu['RA_CEN'][to_process], sky_map_hdu['DE_CEN'][to_process], 'rx', label=str(len(sky_map_hdu['DE_CEN'][to_process]))+' todo')
    plt.legend(loc=0)
    plt.title(GE_name)
    plt.savefig(p2fig)
    plt.clf()
    if len(sky_map_hdu[to_process])>0:
        for kk in np.arange(0, len(sky_map_hdu[to_process]), 50):
            out_im1 = os.path.join(os.environ['GIT_STMOD'], 'src/esass', 'runs', GE_name + '_processing_'+str(kk).zfill(4)+'.sh')
            f_out = open(out_im1, 'w')
            f_out.write("""#!/bin/bash/ \n""")
            f_out.write(" \n ")

            for sky_tile in sky_map_hdu[to_process][kk: kk+50]:
                sky_tile_id = str(sky_tile['SRVMAP'])
                str_field = str(sky_tile['SRVMAP']).zfill(6)
                indir = os.path.join("/home/idies/workspace/erosim/Uchuu/LCerass/", str_field, GE_name)
                esass_dir = os.path.join(indir, 'eSASS')
                git_dir = os.path.join(os.environ['GIT_STMOD'], 'src/esass' )
                path_2_event_file = os.path.join(indir, 'evt_'+str_field+'.fits')
                if os.path.isfile(path_2_event_file) and os.path.isfile(os.path.join(esass_dir,str_field+"_pipeline_img1.sh")) :
                    f_out.write ("cd "+esass_dir +" \n")
                    f_out.write ("sh "+str_field+"_pipeline_img1.sh"+" \n")
                    f_out.write ("sh "+str_field+"_pipeline_det1.sh"+" \n")
                    f_out.write ("sh "+str_field+"_pipeline_Src1.sh"+" \n")
                    f_out.write ("cd "+git_dir+" \n")
                    if not GE_name=='GE_e4_merge_SimBKG':
                        f_out.write ("python photon_matching_RS.py "+GE_name+" "+str_field +" \n")
                    f_out.write('# ====='+" \n")
            f_out.close()
            print(out_im1, 'written')

print('='*100)
print('fields')
print(np.transpose(datata)[0].sum(), 'todo')
print(np.transpose(datata)[1].sum(), 'done')
