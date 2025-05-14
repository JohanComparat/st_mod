#!/bin/bash

# only need to do the detection chain as the flares have been removed from the event files
# script that generates teh scripts for the real data processing

python eRASSX_write_scripts.py 4 030108 # DONE

cd /data56s/comparat/erosim/data_s4_c030

cd /home/comparat/software/erass5_c30_processing/src/esass


# DONE
# execute the scripts to make bash scripts
import os
import numpy as n
from astropy.table import Table, Column
sky_map_hdu = Table.read(os.path.join(os.environ['GIT_ERASS_SIM'], 'data', 'SKYMAPS.fits'))
sky_tile_id = sky_map_hdu['SRVMAP'][(sky_map_hdu['OWNER']==2)|(sky_map_hdu['OWNER']==0)]
for tile_id in sky_tile_id:
    cmd = 'python eRASSX_write_scripts.py 4 '+ str(tile_id).zfill(6)
    os.system(cmd)
    cmd = 'python eRASSX_write_scripts.py 5 '+ str(tile_id).zfill(6)
    os.system(cmd)
# DONE
mkdir /home/comparat/software/erass5_c30_processing/src/esass/data_s4_proc

cd /home/comparat/software/erass5_c30_processing/src/esass/data_s4_proc
# DONE
import os, glob
import numpy as n
from astropy.table import Table, Column
sky_map_hdu = Table.read(os.path.join(os.environ['GIT_ERASS_SIM'], 'data', 'SKYMAPS.fits'))
sky_tile_id = sky_map_hdu['SRVMAP'][ ((sky_map_hdu['OWNER']==2)|(sky_map_hdu['OWNER']==0)) &( abs( sky_map_hdu['GLAT_CEN'] )>=20 )]
# for seed in n.arange(1,100,1):
dir_2_simPh2 = "/data56s/comparat/erosim/data_s4_c030"
for jj_0 in n.arange(0, len(sky_tile_id), 100) :
    print("nohup sh data_s4_c030_processing_on_SKYMAPS_EXGAL_"+str(jj_0).zfill(5)+".sh > logs/data_s4_c030_processing_on_SKYMAPS_EXGAL_"+str(jj_0).zfill(5)+".log &")
    f = open("data_s4_c030_processing_on_SKYMAPS_EXGAL_"+str(jj_0).zfill(5)+".sh", 'w')
    f.write("#!/bin/bash")
    f.write('\n')
    for tile_id in sky_tile_id[jj_0:jj_0+100]:
        dir_Data = dir_2_simPh2+"/" + str(tile_id).zfill(6) + '/eSASS'
        indir = dir_2_simPh2+"/" + str(tile_id).zfill(6) + '/c030'
        # input files
        EvtFiles_list = n.array(glob.glob(os.path.join(indir, '*Image_c030.fits.gz')))
        #print(EvtFiles_list)
        if len(EvtFiles_list) == 1:
            f.write("cd "+dir_Data )
            f.write('\n')
            f.write("sh " + str(tile_id).zfill(6) +"_pipeline_img1_RS.sh")
            f.write('\n')
            f.write("sh " + str(tile_id).zfill(6) +"_pipeline_det1_RS.sh")
            f.write('\n')
            f.write("sh " + str(tile_id).zfill(6) +"_pipeline_Src1_RS.sh")
            f.write('\n')
            f.write('\n')
    f.close()
# DONE
sky_tile_id = sky_map_hdu['SRVMAP'][((sky_map_hdu['OWNER']==2)|(sky_map_hdu['OWNER']==0))&(abs(sky_map_hdu['GLAT_CEN'])<20)]
# for seed in n.arange(1,100,1):
dir_2_simPh2 = "/data56s/comparat/erosim/data_s4_c030"
for jj_0 in n.arange(0, len(sky_tile_id), 100) :
    print("nohup sh data_s4_c030_processing_on_SKYMAPS_INGAL_"+str(jj_0).zfill(5)+".sh > logs/data_s4_c030_processing_on_SKYMAPS_INGAL_"+str(jj_0).zfill(5)+".log &")
    f = open("data_s4_c030_processing_on_SKYMAPS_INGAL_"+str(jj_0).zfill(5)+".sh", 'w')
    f.write("#!/bin/bash")
    f.write('\n')
    for tile_id in sky_tile_id[jj_0:jj_0+100]:
        dir_Data = dir_2_simPh2+"/" + str(tile_id).zfill(6) + '/eSASS'
        indir = dir_2_simPh2+"/" + str(tile_id).zfill(6) + '/c030'
        # input files
        EvtFiles_list = n.array(glob.glob(os.path.join(indir, '*Image_c030.fits.gz')))
        #print(EvtFiles_list)
        if len(EvtFiles_list) == 1:
            f.write("cd "+dir_Data )
            f.write('\n')
            f.write("sh " + str(tile_id).zfill(6) +"_pipeline_img1_RS.sh")
            f.write('\n')
            f.write("sh " + str(tile_id).zfill(6) +"_pipeline_det1_RS.sh")
            f.write('\n')
            f.write("sh " + str(tile_id).zfill(6) +"_pipeline_Src1_RS.sh")
            f.write('\n')
            f.write('\n')
    f.close()
# DONE
mkdir /home/comparat/software/erass5_c30_processing/src/esass/data_s4_proc/logs
# DONE

source /home/erosita/sw/eSASSusers_240410/bin/esass-init.sh
cd /home/comparat/software/erass5_c30_processing/src/esass/data_s4_proc
# nohup sh data_s4_c030_processing_on_SKYMAPS_EXGAL_00000.sh > logs/data_s4_c030_processing_on_SKYMAPS_EXGAL_00000.log & # DONE
# nohup sh data_s4_c030_processing_on_SKYMAPS_EXGAL_00100.sh > logs/data_s4_c030_processing_on_SKYMAPS_EXGAL_00100.log & # DONE
# nohup sh data_s4_c030_processing_on_SKYMAPS_EXGAL_00200.sh > logs/data_s4_c030_processing_on_SKYMAPS_EXGAL_00200.log & # DONE
# nohup sh data_s4_c030_processing_on_SKYMAPS_EXGAL_00300.sh > logs/data_s4_c030_processing_on_SKYMAPS_EXGAL_00300.log & # DONE
# nohup sh data_s4_c030_processing_on_SKYMAPS_EXGAL_00400.sh > logs/data_s4_c030_processing_on_SKYMAPS_EXGAL_00400.log & # DONE
# nohup sh data_s4_c030_processing_on_SKYMAPS_EXGAL_00500.sh > logs/data_s4_c030_processing_on_SKYMAPS_EXGAL_00500.log & # DONE
# nohup sh data_s4_c030_processing_on_SKYMAPS_EXGAL_00600.sh > logs/data_s4_c030_processing_on_SKYMAPS_EXGAL_00600.log & # DONE
# nohup sh data_s4_c030_processing_on_SKYMAPS_EXGAL_00700.sh > logs/data_s4_c030_processing_on_SKYMAPS_EXGAL_00700.log & # DONE
# nohup sh data_s4_c030_processing_on_SKYMAPS_EXGAL_00800.sh > logs/data_s4_c030_processing_on_SKYMAPS_EXGAL_00800.log & # DONE
# nohup sh data_s4_c030_processing_on_SKYMAPS_EXGAL_00900.sh > logs/data_s4_c030_processing_on_SKYMAPS_EXGAL_00900.log & # DONE
# nohup sh data_s4_c030_processing_on_SKYMAPS_EXGAL_01000.sh > logs/data_s4_c030_processing_on_SKYMAPS_EXGAL_01000.log & # DONE
# nohup sh data_s4_c030_processing_on_SKYMAPS_EXGAL_01100.sh > logs/data_s4_c030_processing_on_SKYMAPS_EXGAL_01100.log & # DONE
# nohup sh data_s4_c030_processing_on_SKYMAPS_EXGAL_01200.sh > logs/data_s4_c030_processing_on_SKYMAPS_EXGAL_01200.log & # DONE
# nohup sh data_s4_c030_processing_on_SKYMAPS_EXGAL_01300.sh > logs/data_s4_c030_processing_on_SKYMAPS_EXGAL_01300.log & # DONE
# nohup sh data_s4_c030_processing_on_SKYMAPS_EXGAL_01400.sh > logs/data_s4_c030_processing_on_SKYMAPS_EXGAL_01400.log & # DONE
# nohup sh data_s4_c030_processing_on_SKYMAPS_EXGAL_01500.sh > logs/data_s4_c030_processing_on_SKYMAPS_EXGAL_01500.log & # DONE
# nohup sh data_s4_c030_processing_on_SKYMAPS_EXGAL_01600.sh > logs/data_s4_c030_processing_on_SKYMAPS_EXGAL_01600.log & # DONE
# nohup sh data_s4_c030_processing_on_SKYMAPS_INGAL_00000.sh > logs/data_s4_c030_processing_on_SKYMAPS_INGAL_00000.log & # DONE
# nohup sh data_s4_c030_processing_on_SKYMAPS_INGAL_00100.sh > logs/data_s4_c030_processing_on_SKYMAPS_INGAL_00100.log & # DONE
# nohup sh data_s4_c030_processing_on_SKYMAPS_INGAL_00200.sh > logs/data_s4_c030_processing_on_SKYMAPS_INGAL_00200.log & # DONE
# nohup sh data_s4_c030_processing_on_SKYMAPS_INGAL_00300.sh > logs/data_s4_c030_processing_on_SKYMAPS_INGAL_00300.log & # DONE
# nohup sh data_s4_c030_processing_on_SKYMAPS_INGAL_00400.sh > logs/data_s4_c030_processing_on_SKYMAPS_INGAL_00400.log & # DONE
# nohup sh data_s4_c030_processing_on_SKYMAPS_INGAL_00500.sh > logs/data_s4_c030_processing_on_SKYMAPS_INGAL_00500.log & # DONE
# nohup sh data_s4_c030_processing_on_SKYMAPS_INGAL_00600.sh > logs/data_s4_c030_processing_on_SKYMAPS_INGAL_00600.log & # DONE
# nohup sh data_s4_c030_processing_on_SKYMAPS_INGAL_00700.sh > logs/data_s4_c030_processing_on_SKYMAPS_INGAL_00700.log & # DONE
# nohup sh data_s4_c030_processing_on_SKYMAPS_INGAL_00800.sh > logs/data_s4_c030_processing_on_SKYMAPS_INGAL_00800.log & # DONE

#
# verify which are finished, which are not.
# logging all files and if they are present in a SKYMAP file.
#
# DONE
import os, glob, sys
import numpy as n
from astropy.table import Table, Column
sky_map_hdu = Table.read(os.path.join(os.environ['GIT_ERASS_SIM'], 'data', 'SKYMAPS.fits'))
sky_map_hdu['has_events']      = (1.==0.)
sky_map_hdu['has_ExpMapFiles'] = (1.==0.)
sky_map_hdu['has_UnVigExpMap'] = (1.==0.)
sky_map_hdu['has_DetMask']     = (1.==0.)
sky_map_hdu['has_SrcMapFiles'] = (1.==0.)
sky_map_hdu['has_BkgMapFiles'] = (1.==0.)
sky_map_hdu['has_CheMskFiles'] = (1.==0.)
sky_map_hdu['has_BoxCats']     = (1.==0.)
sky_map_hdu['has_MLCats']      = (1.==0.)
sky_map_hdu['has_ExtCat']      = (1.==0.)
sky_map_hdu['has_MLinCat']     = (1.==0.)
sky_map_hdu['has_ApeCat']      = (1.==0.)
sky_map_hdu['has_ApeSenFiles'] = (1.==0.)
sky_map_hdu['has_PSFMapFiles'] = (1.==0.)
sky_map_hdu['has_SrcCats']     = (1.==0.)
sky_map_hdu['has_SrcReg']      = (1.==0.)
sky_map_hdu['has_BkgReg']      = (1.==0.)
eRASSn = '4'
VerBand4 = str(4)
VerBand = VerBand4
Nmlin = 3
for jj, tile_id_0 in enumerate(sky_map_hdu['SRVMAP']):
    tile_id = str(tile_id_0).zfill(6)
    indir = os.path.join("/data56s/comparat/erosim/data_s"+eRASSn+"_c030", tile_id, 'c030')
    outdir = os.path.join("/data56s/comparat/erosim/data_s"+eRASSn+"_c030", tile_id, 'eSASS')
    EvtFiles_list = n.array(glob.glob(os.path.join(indir, '*Image_c030.fits.gz')))
    sky_map_hdu['has_events'][jj] = len(EvtFiles_list) == 1
    # print(EvtFiles_list)
    if len(EvtFiles_list) == 1:
        outprefix = tile_id + "_"  # ""
        ExpMapFiles = os.path.join(outdir, f"{outprefix}02{VerBand}_ExpMap.fits")
        UnVigExpMap = os.path.join(outdir, f"{outprefix}02{VerBand}_UnvExpMap.fits")
        DetMask     = os.path.join(outdir, f"{outprefix}02{VerBand}_DetMsk.fits")
        SrcMapFiles = {n: os.path.join(outdir, f"{outprefix}02{VerBand}_Src{n:d}Map.fits") for n in (1, 2, 3)}
        BkgMapFiles = {n: os.path.join(outdir, f"{outprefix}02{VerBand}_Bg{n:d}Map.fits") for n in (1, 2, 3)}
        CheMskFiles = os.path.join(outdir, f"{outprefix}02{VerBand}_CheMsk.fits")
        BoxCats     = {n: os.path.join(outdir, f"{outprefix}02{VerBand}_Bo{n:d}Cat.fits") for n in (1, 2, 3, 4, 5)}
        MLCats      = {n: os.path.join(outdir, f"{outprefix}02{VerBand}_ML{n}Cat.fits") for n in (1, 2, 'A')} #uses only 1
        ExtCat      = os.path.join(outdir, f"{outprefix}02{VerBand}_ExtCat.fits")
        MLinCat     = {1: BoxCats[Nmlin], 2: os.path.join(outdir, f"{outprefix}02{VerBand}_ML1Cat.fits"), 'A': os.path.join(outdir, f"{outprefix}02{VerBand}_ML2Cat.fits")}
        ApeCat      = lambda eef: os.path.join(outdir, f"{outprefix}02{VerBand}_ApeCat_eef" + eef + ".fits")
        ApeSenFiles = lambda eef: os.path.join(outdir, f"{outprefix}02{VerBand}_ApeSen_eef" + eef + ".fits")
        PSFMapFiles = lambda eef: os.path.join(outdir, f"{outprefix}02{VerBand}_PSFMap_eef" + eef + ".fits")
        SrcCats     = {n: os.path.join(outdir, f"{outprefix}02{VerBand}_Sc{n:d}Cat.fits") for n in (1, 2, 3)} #only use 1
        SrcReg      = os.path.join(outdir, f"{outprefix}_02_23_srcAUTO.reg")
        BkgReg      = os.path.join(outdir, f"{outprefix}_02_23_bkgAUTO.reg")
        sky_map_hdu['has_ExpMapFiles'][jj] = os.path.isfile( ExpMapFiles )
        sky_map_hdu['has_UnVigExpMap'][jj] = os.path.isfile( UnVigExpMap )
        sky_map_hdu['has_DetMask']    [jj] = os.path.isfile( DetMask     )
        sky_map_hdu['has_SrcMapFiles'][jj] = os.path.isfile( SrcMapFiles[3] )
        sky_map_hdu['has_BkgMapFiles'][jj] = os.path.isfile( BkgMapFiles[3] )
        sky_map_hdu['has_CheMskFiles'][jj] = os.path.isfile( CheMskFiles )
        sky_map_hdu['has_BoxCats']    [jj] = os.path.isfile( BoxCats[3]     )
        sky_map_hdu['has_MLCats']     [jj] = os.path.isfile( MLCats[1]      )
        sky_map_hdu['has_ExtCat']     [jj] = os.path.isfile( ExtCat      )
        sky_map_hdu['has_MLinCat']    [jj] = os.path.isfile( MLinCat[1]     )
        sky_map_hdu['has_ApeCat']     [jj] = os.path.isfile( ApeCat("0.75")      )
        sky_map_hdu['has_ApeSenFiles'][jj] = os.path.isfile( ApeSenFiles("0.75") )
        sky_map_hdu['has_PSFMapFiles'][jj] = os.path.isfile( PSFMapFiles("0.75") )
        sky_map_hdu['has_SrcCats']    [jj] = os.path.isfile( SrcCats[1]     )
        sky_map_hdu['has_SrcReg']     [jj] = os.path.isfile( SrcReg      )
        sky_map_hdu['has_BkgReg']     [jj] = os.path.isfile( BkgReg      )
    else:
        print(indir)

        file_out = "/data56s/comparat/erosim/data_s"+eRASSn+"_c030/SKYMAPS_eSASS_files.fits"
sky_map_hdu.write(file_out, overwrite = True)
# DONE
# cd $DATA_S4
# rsync comparat@ds43:/data56s/comparat/erosim/data_s4_c030/SKYMAPS_eSASS_files.fits .
# topcat SKYMAPS_eSASS_files.fits


# DONE
# create commands for the remaining ones:
cd /home/comparat/software/erass5_c30_processing/src/esass/data_s4_proc

import os, glob
import numpy as n
from astropy.table import Table, Column
file_out = "/data56s/comparat/erosim/data_s4_c030/SKYMAPS_eSASS_files.fits"
sky_map_hdu = Table.read(file_out)
sky_tile_id = sky_map_hdu['SRVMAP'][ ((sky_map_hdu['OWNER']==2)|(sky_map_hdu['OWNER']==0)) &( sky_map_hdu['has_events'] ) & ( sky_map_hdu['has_BkgReg']==False )]
sky_map_todo = sky_map_hdu[ ((sky_map_hdu['OWNER']==2)|(sky_map_hdu['OWNER']==0)) &( sky_map_hdu['has_events'] ) & ( sky_map_hdu['has_BkgReg']==False )]
# for seed in n.arange(1,100,1):
dir_2_simPh2 = "/data56s/comparat/erosim/data_s4_c030"
for jj_0 in n.arange(0, len(sky_tile_id), 50) :
    print("nohup sh remainingFinal_data_s4_c030_processing_on_SKYMAPS_"+str(jj_0).zfill(5)+".sh > logs/remainingFinal_data_s4_c030_processing_on_SKYMAPS_"+str(jj_0).zfill(5)+".log &")
    f = open("remainingFinal_data_s4_c030_processing_on_SKYMAPS_"+str(jj_0).zfill(5)+".sh", 'w')
    f.write("#!/bin/bash")
    f.write('\n')
    for n_tile in n.arange(len(sky_map_todo[jj_0:jj_0+50])):
        sky_tile = sky_map_todo[jj_0:jj_0+50][n_tile]
        tile_id = sky_tile['SRVMAP']
        dir_Data = dir_2_simPh2+"/" + str(tile_id).zfill(6) + '/eSASS'
        indir = dir_2_simPh2+"/" + str(tile_id).zfill(6) + '/c030'
        # input files
        # print(dir_Data)
        EvtFiles_list = n.array(glob.glob(os.path.join(indir, '*Image_c030.fits.gz')))
        #print(EvtFiles_list)
        if len(EvtFiles_list) == 1:
            f.write("cd "+dir_Data )
            f.write('\n')
            if sky_tile['has_ExpMapFiles']:
                f.write('rm *.fits \n')
                f.write('rm *.reg \n')
            f.write('\n')
            f.write("sh " + str(tile_id).zfill(6) +"_pipeline_img1_RS.sh")
            f.write('\n')
            f.write("sh " + str(tile_id).zfill(6) +"_pipeline_det1_RS.sh")
            f.write('\n')
            f.write("sh " + str(tile_id).zfill(6) +"_pipeline_Src1_RS.sh")
            f.write('\n')
            f.write('\n')
    f.close()
# 0 remaining
# DONE

# data s4 processing finished.


