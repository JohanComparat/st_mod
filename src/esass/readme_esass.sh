#!/bin/bash/
export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/esass/

szr16jdjsgpd:q5=v0jJ

 out=Table.read(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS_'+'GE_e4_merge_AGNseed001_SimBKG_CLUseed001'+'.fits'))
 ero_de = (out['OWNER']==2)|(out['OWNER']==0)
 to_process = ((out['OWNER'] == 2) | (out['OWNER'] == 0)) & (out['has_merged_events']) & (out['has_Sc1Cat'] == False)
 already_done = ((out['OWNER'] == 2) | (out['OWNER'] == 0)) & (out['has_merged_events']) & (out['has_Sc1Cat'])
 todo = (ero_de)&(to_process==False)&(already_done==False)
 out['fail'] = 0
 out['fail'][todo] = np.array(fails)
 out.write(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS_'+'GE_e4_merge_AGNseed001_SimBKG_CLUseed001'+'_withFailReason.fits'), overwrite = True)
 out = Table.read(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS_'+'GE_e4_merge_AGNseed001_SimBKG_CLUseed001'+'_withFailReason.fits'))
 tosync_srv_map = out['SRVMAP'][(out['fail']==1)]
 for el in tosync_srv_map:
  	str_field = str(el).zfill(6)
  	print('mv ' + os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'c030', '*') + ' ' + os.path.join(os.environ['UCHUU'], LC_dir, str_field, 's4_c030') )
 for el in tosync_srv_map:
  	str_field = str(el).zfill(6)
  	print('mv ' + os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'eSASS', '*') + ' ' + os.path.join(os.environ['UCHUU'], LC_dir, str_field, 's4_eSASS') )

cd ~/workspace/erosim/Uchuu/LCerass
rsync -avz joco@raven.mpcdf.mpg.de:~/ptmp_joco/mpecl/comparat/data_s4_c030/00???? . # DONE
rsync -avz joco@raven.mpcdf.mpg.de:~/ptmp_joco/mpecl/comparat/data_s4_c030/0????? . # DONE

cd ~/workspace/erosim/Uchuu/LCerass
rsync -avz joco@raven.mpcdf.mpg.de:~/ptmp_joco/mpecl/comparat/data_s5_c030/00???? . # TODO
rsync -avz joco@raven.mpcdf.mpg.de:~/ptmp_joco/mpecl/comparat/data_s5_c030/0????? . # TODO

export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/esass/
sh mv_data.sh

BG model files are here :
~/workspace/erosim/simput/bkg_erosita_simput_full_sky
catalog.fits has the full map

# make a sky map with input files
python make_summarySimEvt_skymap.py # DONE, OK, all files are there.

export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/esass/
python monitor_merge.py # ONGOING

export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/esass/
python monitor_merge_noCLU.py # ONGOING

export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/esass/
python monitor_merge_noAGN.py # DONE

## in any container with python
#nohup python merge_events_onlyBG.py  > logs/merge_events_onlyBG.log  & # DONE
#
#nohup python merge_events.py 1 1     > logs/merge_events_1_1.log  & # DONE
#nohup python merge_events_noCLU.py 1 > logs/merge_events_noCLU_1.log & # DONE
#nohup python merge_events_noAGN.py 1 > logs/merge_events_noAGN_1.log & # DONE
#
#nohup python merge_events.py 2 2     > logs/merge_events_2_2.log  & # DONE
#nohup python merge_events_noCLU.py 2 > logs/merge_events_noCLU_2.log & # DONE
#nohup python merge_events_noAGN.py 2 > logs/merge_events_noAGN_2.log & # DONE
#
#nohup python merge_events.py 3 3     > logs/merge_events_3_3.log  & # ONGOING
#nohup python merge_events_noCLU.py 3 > logs/merge_events_noCLU_3.log & # DONE
#nohup python merge_events_noAGN.py 3 > logs/merge_events_noAGN_3.log & # DONE
#
#nohup python merge_events.py 4 4     > logs/merge_events_4_4.log  & # ONGOING
#nohup python merge_events_noCLU.py 4 > logs/merge_events_noCLU_4.log & # DONE
#nohup python merge_events_noAGN.py 4 > logs/merge_events_noAGN_4.log & # DONE
#
#nohup python merge_events.py 5 5     > logs/merge_events_5_5.log  & # ONGOING
#nohup python merge_events_noCLU.py 5 > logs/merge_events_noCLU_5.log & # DONE
#nohup python merge_events_noAGN.py 5 > logs/merge_events_noAGN_5.log & # DONE
#
#nohup python merge_events.py 6 6     > logs/merge_events_6_6.log  & # ONGOING
#nohup python merge_events_noCLU.py 6 > logs/merge_events_noCLU_6.log & # DONE
#nohup python merge_events_noAGN.py 6 > logs/merge_events_noAGN_6.log & # DONE
#
#nohup python merge_events.py 7 7     > logs/merge_events_7_7.log  & # ONGOING
#nohup python merge_events_noCLU.py 7 > logs/merge_events_noCLU_7.log & # DONE
#nohup python merge_events_noAGN.py 7 > logs/merge_events_noAGN_7.log & # DONE
#
#nohup python merge_events.py 8 8     > logs/merge_events_8_8.log  & # ONGOING
#nohup python merge_events_noCLU.py 8 > logs/merge_events_noCLU_8.log & # DONE
#nohup python merge_events_noAGN.py 8 > logs/merge_events_noAGN_8.log & # DONE


GE_e4_merge_AGNseed00?_SimBKG : AGN + BKG
GE_e4_merge_AGNseed00?_SimBKG_CLUseed00? : AGN + BKG + Cluster
GE_e4_merge_SimBKG : only BKG
GE_e4_merge_SimBKG_CLUseed00? : BKG + Cluster

## write eSASS commands in the relevant folders
#nohup python eRASSX_write_scripts.py GE_e4_merge_SimBKG                       > logs/eRASSX_write_scripts_GE_e4_merge_SimBKG.log                       & # DONE
#
#nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed001_SimBKG            > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed001_SimBKG.log            & # DONE
#nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed001_SimBKG_CLUseed001 > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed001_SimBKG_CLUseed001.log & # DONE
#nohup python eRASSX_write_scripts.py GE_e4_merge_SimBKG_CLUseed001            > logs/eRASSX_write_scripts_GE_e4_merge_SimBKG_CLUseed001.log            & # DONE
#
#nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed002_SimBKG            > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed002_SimBKG.log            & # DONE
#nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed002_SimBKG_CLUseed002 > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed002_SimBKG_CLUseed002.log & # DONE
#nohup python eRASSX_write_scripts.py GE_e4_merge_SimBKG_CLUseed002            > logs/eRASSX_write_scripts_GE_e4_merge_SimBKG_CLUseed002.log            & # DONE
#
#nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed003_SimBKG            > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed003_SimBKG.log            & # DONE
#nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed003_SimBKG_CLUseed003 > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed003_SimBKG_CLUseed003.log & # DONE
#nohup python eRASSX_write_scripts.py GE_e4_merge_SimBKG_CLUseed003            > logs/eRASSX_write_scripts_GE_e4_merge_SimBKG_CLUseed003.log            & # DONE
#
#nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed004_SimBKG            > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed004_SimBKG.log            & # DONE
#nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed004_SimBKG_CLUseed004 > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed004_SimBKG_CLUseed004.log & # DONE
#nohup python eRASSX_write_scripts.py GE_e4_merge_SimBKG_CLUseed004            > logs/eRASSX_write_scripts_GE_e4_merge_SimBKG_CLUseed004.log            & # DONE
#
#nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed005_SimBKG            > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed005_SimBKG.log            & # DONE
#nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed005_SimBKG_CLUseed005 > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed005_SimBKG_CLUseed005.log & # DONE
#nohup python eRASSX_write_scripts.py GE_e4_merge_SimBKG_CLUseed005            > logs/eRASSX_write_scripts_GE_e4_merge_SimBKG_CLUseed005.log            & # DONE
#
#nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed006_SimBKG            > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed006_SimBKG.log            & # DONE
#nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed006_SimBKG_CLUseed006 > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed006_SimBKG_CLUseed006.log & # DONE
#nohup python eRASSX_write_scripts.py GE_e4_merge_SimBKG_CLUseed006            > logs/eRASSX_write_scripts_GE_e4_merge_SimBKG_CLUseed006.log            & # DONE
#
#nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed007_SimBKG            > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed007_SimBKG.log            & # DONE
#nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed007_SimBKG_CLUseed007 > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed007_SimBKG_CLUseed007.log & # DONE
#nohup python eRASSX_write_scripts.py GE_e4_merge_SimBKG_CLUseed007            > logs/eRASSX_write_scripts_GE_e4_merge_SimBKG_CLUseed007.log            & # DONE
#
#nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed008_SimBKG            > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed008_SimBKG.log            & # DONE
#nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed008_SimBKG_CLUseed008 > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed008_SimBKG_CLUseed008.log & # DONE
#nohup python eRASSX_write_scripts.py GE_e4_merge_SimBKG_CLUseed008            > logs/eRASSX_write_scripts_GE_e4_merge_SimBKG_CLUseed008.log            & # DONE

# TODO

#monitor
# monitor progress
export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/esass/
python monitor.py # TODO
# runs summary_file_run.sh

#nohup python make_summary_skymap.py GE_e4_merge_SimBKG                       > logs/summary_sky_map_GE_e4_merge_SimBKG.log            & # TODO
#
#nohup python make_summary_skymap.py GE_e4_merge_AGNseed001_SimBKG            > logs/summary_sky_map_GE_e4_merge_AGNseed001_SimBKG.log            & # TODO
#nohup python make_summary_skymap.py GE_e4_merge_AGNseed001_SimBKG_CLUseed001 > logs/summary_sky_map_GE_e4_merge_AGNseed001_SimBKG_CLUseed001.log & # TODO
#nohup python make_summary_skymap.py GE_e4_merge_SimBKG_CLUseed001            > logs/summary_sky_map_GE_e4_merge_SimBKG_CLUseed001.log            & # TODO
#
#nohup python make_summary_skymap.py GE_e4_merge_AGNseed002_SimBKG            > logs/summary_sky_map_GE_e4_merge_AGNseed002_SimBKG.log            & # TODO
#nohup python make_summary_skymap.py GE_e4_merge_AGNseed002_SimBKG_CLUseed002 > logs/summary_sky_map_GE_e4_merge_AGNseed002_SimBKG_CLUseed002.log & # TODO
#nohup python make_summary_skymap.py GE_e4_merge_SimBKG_CLUseed002            > logs/summary_sky_map_GE_e4_merge_SimBKG_CLUseed002.log            & # TODO
#
#nohup python make_summary_skymap.py GE_e4_merge_AGNseed003_SimBKG            > logs/summary_sky_map_GE_e4_merge_AGNseed003_SimBKG.log            & # TODO
#nohup python make_summary_skymap.py GE_e4_merge_AGNseed003_SimBKG_CLUseed003 > logs/summary_sky_map_GE_e4_merge_AGNseed003_SimBKG_CLUseed003.log & # TODO
#nohup python make_summary_skymap.py GE_e4_merge_SimBKG_CLUseed003            > logs/summary_sky_map_GE_e4_merge_SimBKG_CLUseed003.log            & # TODO
#
#nohup python make_summary_skymap.py GE_e4_merge_AGNseed004_SimBKG            > logs/summary_skymap_GE_e4_merge_AGNseed004_SimBKG.log            & # TODO
#nohup python make_summary_skymap.py GE_e4_merge_AGNseed004_SimBKG_CLUseed004 > logs/summary_skymap_GE_e4_merge_AGNseed004_SimBKG_CLUseed004.log & # TODO
#nohup python make_summary_skymap.py GE_e4_merge_SimBKG_CLUseed004            > logs/summary_skymap_GE_e4_merge_SimBKG_CLUseed004.log            & # TODO
#
#nohup python make_summary_skymap.py GE_e4_merge_AGNseed005_SimBKG            > logs/summary_skymap_GE_e4_merge_AGNseed005_SimBKG.log            & # TODO
#nohup python make_summary_skymap.py GE_e4_merge_AGNseed005_SimBKG_CLUseed005 > logs/summary_skymap_GE_e4_merge_AGNseed005_SimBKG_CLUseed005.log & # TODO
#nohup python make_summary_skymap.py GE_e4_merge_SimBKG_CLUseed005            > logs/summary_skymap_GE_e4_merge_SimBKG_CLUseed005.log            & # TODO
#
#nohup python make_summary_skymap.py GE_e4_merge_AGNseed006_SimBKG            > logs/summary_skymap_GE_e4_merge_AGNseed006_SimBKG.log            & # TODO
#nohup python make_summary_skymap.py GE_e4_merge_AGNseed006_SimBKG_CLUseed006 > logs/summary_skymap_GE_e4_merge_AGNseed006_SimBKG_CLUseed006.log & # TODO
#nohup python make_summary_skymap.py GE_e4_merge_SimBKG_CLUseed006            > logs/summary_skymap_GE_e4_merge_SimBKG_CLUseed006.log            & # TODO
#
#nohup python make_summary_skymap.py GE_e4_merge_AGNseed007_SimBKG            > logs/summary_skymap_GE_e4_merge_AGNseed007_SimBKG.log            & # TODO
#nohup python make_summary_skymap.py GE_e4_merge_AGNseed007_SimBKG_CLUseed007 > logs/summary_skymap_GE_e4_merge_AGNseed007_SimBKG_CLUseed007.log & # TODO
#nohup python make_summary_skymap.py GE_e4_merge_SimBKG_CLUseed007            > logs/summary_skymap_GE_e4_merge_SimBKG_CLUseed007.log            & # TODO
#
#nohup python make_summary_skymap.py GE_e4_merge_AGNseed008_SimBKG            > logs/summary_skymap_GE_e4_merge_AGNseed008_SimBKG.log            & # TODO
#nohup python make_summary_skymap.py GE_e4_merge_AGNseed008_SimBKG_CLUseed008 > logs/summary_skymap_GE_e4_merge_AGNseed008_SimBKG_CLUseed008.log & # TODO
#nohup python make_summary_skymap.py GE_e4_merge_SimBKG_CLUseed008            > logs/summary_skymap_GE_e4_merge_SimBKG_CLUseed008.log            & # TODO

# TODO

# write all commands for the list of folders above.
export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/esass/
python write_exec_loop.py # > exec.sh # TODO

# TODO

source activate heasoft
[ -r /home/idies/.healpix/3_50_Linux/config ] && . /home/idies/.healpix/3_50_Linux/config
source /opt/esass/bin/esass-init.sh
export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/esass/
python monitor_esass_run_001.py #

source activate heasoft
[ -r /home/idies/.healpix/3_50_Linux/config ] && . /home/idies/.healpix/3_50_Linux/config
source /opt/esass/bin/esass-init.sh
export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/esass/
python monitor_esass_run_002.py

-rw-r--r--. 1 idies idies 9.9K Jul  7 14:21 GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_0000.sh
-rw-r--r--. 1 idies idies 9.9K Jul  7 14:21 GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_0030.sh
-rw-r--r--. 1 idies idies 2.2K Jul  7 14:21 GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_0060.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:21 GE_e4_merge_AGNseed001_SimBKG_processing_0000.sh
-rw-r--r--. 1 idies idies 8.1K Jul  7 14:21 GE_e4_merge_AGNseed001_SimBKG_processing_0030.sh
-rw-r--r--. 1 idies idies 9.9K Jul  7 14:22 GE_e4_merge_AGNseed002_SimBKG_CLUseed002_processing_0000.sh
-rw-r--r--. 1 idies idies 9.9K Jul  7 14:22 GE_e4_merge_AGNseed002_SimBKG_CLUseed002_processing_0030.sh
-rw-r--r--. 1 idies idies 1.2K Jul  7 14:22 GE_e4_merge_AGNseed002_SimBKG_CLUseed002_processing_0060.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:22 GE_e4_merge_AGNseed002_SimBKG_processing_0000.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:22 GE_e4_merge_AGNseed002_SimBKG_processing_0030.sh
-rw-r--r--. 1 idies idies 7.5K Jul  7 14:22 GE_e4_merge_AGNseed002_SimBKG_processing_0060.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:21 GE_e4_merge_SimBKG_CLUseed001_processing_0000.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:22 GE_e4_merge_SimBKG_CLUseed002_processing_0000.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:22 GE_e4_merge_SimBKG_CLUseed002_processing_0030.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:22 GE_e4_merge_SimBKG_CLUseed002_processing_0060.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:22 GE_e4_merge_SimBKG_CLUseed002_processing_0090.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:22 GE_e4_merge_SimBKG_CLUseed002_processing_0120.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:22 GE_e4_merge_SimBKG_CLUseed002_processing_0150.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:22 GE_e4_merge_SimBKG_CLUseed002_processing_0180.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:22 GE_e4_merge_SimBKG_CLUseed002_processing_0210.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:22 GE_e4_merge_SimBKG_CLUseed002_processing_0240.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:22 GE_e4_merge_SimBKG_CLUseed002_processing_0270.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:22 GE_e4_merge_SimBKG_CLUseed002_processing_0300.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:22 GE_e4_merge_SimBKG_CLUseed002_processing_0330.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:22 GE_e4_merge_SimBKG_CLUseed002_processing_0360.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:22 GE_e4_merge_SimBKG_CLUseed002_processing_0390.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:22 GE_e4_merge_SimBKG_CLUseed002_processing_0420.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:22 GE_e4_merge_SimBKG_CLUseed002_processing_0450.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:22 GE_e4_merge_SimBKG_CLUseed002_processing_0480.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:22 GE_e4_merge_SimBKG_CLUseed002_processing_0510.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:22 GE_e4_merge_SimBKG_CLUseed002_processing_0540.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:23 GE_e4_merge_SimBKG_CLUseed002_processing_0570.sh
#-rw-r--r--. 1 idies idies 9.3K Jul  7 14:23 GE_e4_merge_SimBKG_CLUseed002_processing_0600.sh
#-rw-r--r--. 1 idies idies 9.3K Jul  7 14:23 GE_e4_merge_SimBKG_CLUseed002_processing_0630.sh
#-rw-r--r--. 1 idies idies 9.3K Jul  7 14:23 GE_e4_merge_SimBKG_CLUseed002_processing_0660.sh
#-rw-r--r--. 1 idies idies 9.3K Jul  7 14:23 GE_e4_merge_SimBKG_CLUseed002_processing_0690.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:23 GE_e4_merge_SimBKG_CLUseed002_processing_0720.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:23 GE_e4_merge_SimBKG_CLUseed002_processing_0750.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:23 GE_e4_merge_SimBKG_CLUseed002_processing_0780.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:23 GE_e4_merge_SimBKG_CLUseed002_processing_0810.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:23 GE_e4_merge_SimBKG_CLUseed002_processing_0840.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:23 GE_e4_merge_SimBKG_CLUseed002_processing_0870.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:23 GE_e4_merge_SimBKG_CLUseed002_processing_0900.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:23 GE_e4_merge_SimBKG_CLUseed002_processing_0930.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:23 GE_e4_merge_SimBKG_CLUseed002_processing_0960.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:23 GE_e4_merge_SimBKG_CLUseed002_processing_0990.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:23 GE_e4_merge_SimBKG_CLUseed002_processing_1020.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:23 GE_e4_merge_SimBKG_CLUseed002_processing_1050.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:23 GE_e4_merge_SimBKG_CLUseed002_processing_1080.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:23 GE_e4_merge_SimBKG_CLUseed002_processing_1110.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:24 GE_e4_merge_SimBKG_CLUseed002_processing_1140.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:24 GE_e4_merge_SimBKG_CLUseed002_processing_1170.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:24 GE_e4_merge_SimBKG_CLUseed002_processing_1200.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:24 GE_e4_merge_SimBKG_CLUseed002_processing_1230.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:24 GE_e4_merge_SimBKG_CLUseed002_processing_1260.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:24 GE_e4_merge_SimBKG_CLUseed002_processing_1290.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:24 GE_e4_merge_SimBKG_CLUseed002_processing_1320.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:24 GE_e4_merge_SimBKG_CLUseed002_processing_1350.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:24 GE_e4_merge_SimBKG_CLUseed002_processing_1380.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:24 GE_e4_merge_SimBKG_CLUseed002_processing_1410.sh
-rw-r--r--. 1 idies idies 9.3K Jul  7 14:24 GE_e4_merge_SimBKG_CLUseed002_processing_1440.sh
-rw-r--r--. 1 idies idies 9.0K Jul  7 14:24 GE_e4_merge_SimBKG_CLUseed002_processing_1470.sh


# execute all commands of interest in an eSASS loaded container
# in a monitor script !
nohup sh GE_e4_merge_AGNseed002_SimBKG_CLUseed002_processing_0000.sh > logs/AGNseed002_SimBKG_CLUseed002_processing_0000.log & # ONGOING

# TODO
# start realizations 5,6,7,8
# other scaling relation sim


export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/esass/

python create_summary_files_RS.py GE_e4_merge_AGNseed001_SimBKG_CLUseed001
python create_summary_files_RS.py GE_e4_merge_SimBKG_CLUseed001
python create_summary_files_RS.py GE_e4_merge_AGNseed001_SimBKG
python create_summary_files_RS.py GE_e4_merge_SimBKG

python create_summary_files_RS.py GE_e4_merge_AGNseed002_SimBKG_CLUseed002
python create_summary_files_RS.py GE_e4_merge_SimBKG_CLUseed002
python create_summary_files_RS.py GE_e4_merge_AGNseed002_SimBKG

python create_summary_files_RS.py GE_e4_merge_AGNseed003_SimBKG_CLUseed003
python create_summary_files_RS.py GE_e4_merge_SimBKG_CLUseed003
python create_summary_files_RS.py GE_e4_merge_AGNseed003_SimBKG

python create_summary_files_RS.py GE_e4_merge_AGNseed004_SimBKG_CLUseed004
python create_summary_files_RS.py GE_e4_merge_SimBKG_CLUseed004
python create_summary_files_RS.py GE_e4_merge_AGNseed004_SimBKG

python create_summary_files_RS.py GE_e4_merge_AGNseed005_SimBKG_CLUseed005
python create_summary_files_RS.py GE_e4_merge_SimBKG_CLUseed005
python create_summary_files_RS.py GE_e4_merge_AGNseed005_SimBKG

python create_summary_files_RS.py GE_e4_merge_AGNseed006_SimBKG_CLUseed006
python create_summary_files_RS.py GE_e4_merge_SimBKG_CLUseed006
python create_summary_files_RS.py GE_e4_merge_AGNseed006_SimBKG

python create_summary_files_RS.py GE_e4_merge_AGNseed007_SimBKG_CLUseed007
python create_summary_files_RS.py GE_e4_merge_SimBKG_CLUseed007
python create_summary_files_RS.py GE_e4_merge_AGNseed007_SimBKG

python create_summary_files_RS.py GE_e4_merge_AGNseed008_SimBKG_CLUseed008
python create_summary_files_RS.py GE_e4_merge_SimBKG_CLUseed008
python create_summary_files_RS.py GE_e4_merge_AGNseed008_SimBKG

python summary_plotting_RS.py GE_e4_merge_AGNseed001_SimBKG_CLUseed001 /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_AGNseed002_SimBKG_CLUseed002 /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_AGNseed003_SimBKG_CLUseed003 /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_AGNseed004_SimBKG_CLUseed004 /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_AGNseed005_SimBKG_CLUseed005 /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_AGNseed006_SimBKG_CLUseed006 /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_AGNseed007_SimBKG_CLUseed007 /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_AGNseed008_SimBKG_CLUseed008 /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins

python summary_plotting_RS.py GE_e4_merge_SimBKG_CLUseed001 /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_SimBKG_CLUseed002 /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_SimBKG_CLUseed003 /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_SimBKG_CLUseed004 /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_SimBKG_CLUseed005 /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_SimBKG_CLUseed006 /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_SimBKG_CLUseed007 /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_SimBKG_CLUseed008 /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_AGNseed001_SimBKG /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_AGNseed002_SimBKG /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_AGNseed003_SimBKG /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_AGNseed004_SimBKG /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_AGNseed005_SimBKG /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_AGNseed006_SimBKG /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_AGNseed007_SimBKG /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_AGNseed008_SimBKG /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins














#
# test fields and scripts :
#
export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/esass/
cd $GIT_STMOD/src/esass/runs
sh GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_0000.sh  > logs/GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_0000.log & # running
sh GE_e4_merge_SimBKG_CLUseed001_processing_0000.sh             > logs/GE_e4_merge_SimBKG_CLUseed001_processing_0000.log            & # running
sh GE_e4_merge_SimBKG_processing_0000.sh                        > logs/GE_e4_merge_SimBKG_processing_0000.log                       & # running


cd /home/idies/workspace/erosim/Uchuu/LCerass/121048/GE_e4_merge_AGNseed001_SimBKG_CLUseed001/eSASS
sh 121048_pipeline_img1.sh
sh 121048_pipeline_det1.sh
sh 121048_pipeline_Src1.sh
cd /home/idies/workspace/erosim/software/st_mod/src/esass
python photon_matching_RS.py GE_e4_merge_AGNseed001_SimBKG_CLUseed001 121048


# test field
# execute eSASS commands on test field:
cd /home/idies/workspace/erosim/Uchuu/LCerass/121048/GE_e4_merge_AGNseed001_SimBKG/eSASS
sh 121048_pipeline_img1.sh # TODO
sh 121048_pipeline_det1.sh
sh 121048_pipeline_Src1.sh
cd $GIT_STMOD/src/esass
python photon_matching_RS.py GE_e4_merge_AGNseed001_SimBKG 121048

cd /home/idies/workspace/erosim/Uchuu/LCerass/121048/GE_e4_merge_AGNseed001_SimBKG_CLUseed001/eSASS
sh 121048_pipeline_img1.sh # TODO
sh 121048_pipeline_det1.sh
sh 121048_pipeline_Src1.sh
cd $GIT_STMOD/src/esass
python photon_matching_RS.py GE_e4_merge_AGNseed001_SimBKG_CLUseed001 121048

cd /home/idies/workspace/erosim/Uchuu/LCerass/121048/GE_e4_merge_SimBKG/eSASS
sh 121048_pipeline_img1.sh # TODO
sh 121048_pipeline_det1.sh
sh 121048_pipeline_Src1.sh
#python photon_matching_RS.py GE_e4_merge_SimBKG 121048

cd /home/idies/workspace/erosim/Uchuu/LCerass/121048/GE_e4_merge_SimBKG_CLUseed001/eSASS
sh 121048_pipeline_img1.sh # TODO
sh 121048_pipeline_det1.sh
sh 121048_pipeline_Src1.sh
cd $GIT_STMOD/src/esass
python photon_matching_RS.py GE_e4_merge_SimBKG_CLUseed001 121048

