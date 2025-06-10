#!/bin/bash/
export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/esass/

BG model files are here :
~/workspace/erosim/simput/bkg_erosita_simput_full_sky
catalog.fits has the full map

# make a sky map with input files
python make_summarySimEvt_skymap.py # DONE, OK, all files are there.

# in any container with python
nohup python merge_events_onlyBG.py  > logs/merge_events_onlyBG.log  & # ONGOING

nohup python merge_events.py 1 1     > logs/merge_events_1_1.log  & # ONGOING
nohup python merge_events_noCLU.py 1 > logs/merge_events_noCLU_1.log & # DONE
nohup python merge_events_noAGN.py 1 > logs/merge_events_noAGN_1.log & # ONGOING

nohup python merge_events.py 2 2     > logs/merge_events_2_2.log  & # ONGOING
nohup python merge_events_noCLU.py 2 > logs/merge_events_noCLU_2.log & # DONE
nohup python merge_events_noAGN.py 2 > logs/merge_events_noAGN_2.log & # ONGOING

nohup python merge_events.py 3 3     > logs/merge_events_3_3.log  & # ONGOING
nohup python merge_events_noCLU.py 3 > logs/merge_events_noCLU_3.log & # DONE
nohup python merge_events_noAGN.py 3 > logs/merge_events_noAGN_3.log & # ONGOING

nohup python merge_events.py 4 4     > logs/merge_events_4_4.log  & # ONGOING
nohup python merge_events_noCLU.py 4 > logs/merge_events_noCLU_4.log & # DONE
nohup python merge_events_noAGN.py 4 > logs/merge_events_noAGN_4.log & # ONGOING

nohup python merge_events.py 5 5     > logs/merge_events_5_5.log  & # ONGOING
nohup python merge_events_noCLU.py 5 > logs/merge_events_noCLU_5.log & # DONE
nohup python merge_events_noAGN.py 5 > logs/merge_events_noAGN_5.log & # ONGOING

nohup python merge_events.py 6 6     > logs/merge_events_6_6.log  & # ONGOING
nohup python merge_events_noCLU.py 6 > logs/merge_events_noCLU_6.log & # DONE
nohup python merge_events_noAGN.py 6 > logs/merge_events_noAGN_6.log & # ONGOING

nohup python merge_events.py 7 7     > logs/merge_events_7_7.log  & # ONGOING
nohup python merge_events_noCLU.py 7 > logs/merge_events_noCLU_7.log & # DONE
nohup python merge_events_noAGN.py 7 > logs/merge_events_noAGN_7.log & # ONGOING

nohup python merge_events.py 8 8     > logs/merge_events_8_8.log  & # ONGOING
nohup python merge_events_noCLU.py 8 > logs/merge_events_noCLU_8.log & # DONE
nohup python merge_events_noAGN.py 8 > logs/merge_events_noAGN_8.log & # ONGOING


GE_e4_merge_AGNseed00?_SimBKG : AGN + BKG
GE_e4_merge_AGNseed00?_SimBKG_CLUseed00? : AGN + BKG + Cluster
GE_e4_merge_SimBKG : only BKG
GE_e4_merge_SimBKG_CLUseed00? : BKG + Cluster

# write eSASS commands in the relevant folders
nohup python eRASSX_write_scripts.py GE_e4_merge_SimBKG                       > logs/eRASSX_write_scripts_GE_e4_merge_SimBKG.log                       & # DONE

nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed001_SimBKG            > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed001_SimBKG.log            & # DONE
nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed001_SimBKG_CLUseed001 > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed001_SimBKG_CLUseed001.log & # DONE
nohup python eRASSX_write_scripts.py GE_e4_merge_SimBKG_CLUseed001            > logs/eRASSX_write_scripts_GE_e4_merge_SimBKG_CLUseed001.log            & # DONE

nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed002_SimBKG            > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed002_SimBKG.log            & # DONE
nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed002_SimBKG_CLUseed002 > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed002_SimBKG_CLUseed002.log & # DONE
nohup python eRASSX_write_scripts.py GE_e4_merge_SimBKG_CLUseed002            > logs/eRASSX_write_scripts_GE_e4_merge_SimBKG_CLUseed002.log            & # DONE

nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed003_SimBKG            > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed003_SimBKG.log            & # DONE
nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed003_SimBKG_CLUseed003 > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed003_SimBKG_CLUseed003.log & # DONE
nohup python eRASSX_write_scripts.py GE_e4_merge_SimBKG_CLUseed003            > logs/eRASSX_write_scripts_GE_e4_merge_SimBKG_CLUseed003.log            & # DONE

nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed004_SimBKG            > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed004_SimBKG.log            & # DONE
nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed004_SimBKG_CLUseed004 > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed004_SimBKG_CLUseed004.log & # DONE
nohup python eRASSX_write_scripts.py GE_e4_merge_SimBKG_CLUseed004            > logs/eRASSX_write_scripts_GE_e4_merge_SimBKG_CLUseed004.log            & # DONE

nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed005_SimBKG            > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed005_SimBKG.log            & # DONE
nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed005_SimBKG_CLUseed005 > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed005_SimBKG_CLUseed005.log & # DONE
nohup python eRASSX_write_scripts.py GE_e4_merge_SimBKG_CLUseed005            > logs/eRASSX_write_scripts_GE_e4_merge_SimBKG_CLUseed005.log            & # DONE

nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed006_SimBKG            > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed006_SimBKG.log            & # DONE
nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed006_SimBKG_CLUseed006 > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed006_SimBKG_CLUseed006.log & # DONE
nohup python eRASSX_write_scripts.py GE_e4_merge_SimBKG_CLUseed006            > logs/eRASSX_write_scripts_GE_e4_merge_SimBKG_CLUseed006.log            & # DONE

nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed007_SimBKG            > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed007_SimBKG.log            & # DONE
nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed007_SimBKG_CLUseed007 > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed007_SimBKG_CLUseed007.log & # DONE
nohup python eRASSX_write_scripts.py GE_e4_merge_SimBKG_CLUseed007            > logs/eRASSX_write_scripts_GE_e4_merge_SimBKG_CLUseed007.log            & # DONE

nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed008_SimBKG            > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed008_SimBKG.log            & # DONE
nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed008_SimBKG_CLUseed008 > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed008_SimBKG_CLUseed008.log & # DONE
nohup python eRASSX_write_scripts.py GE_e4_merge_SimBKG_CLUseed008            > logs/eRASSX_write_scripts_GE_e4_merge_SimBKG_CLUseed008.log            & # DONE



# monitor progress
export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/esass/

nohup python make_summary_skymap.py GE_e4_merge_SimBKG                       > logs/summary_sky_map_GE_e4_merge_SimBKG.log            & # TODO

nohup python make_summary_skymap.py GE_e4_merge_AGNseed001_SimBKG            > logs/summary_sky_map_GE_e4_merge_AGNseed001_SimBKG.log            & # TODO
nohup python make_summary_skymap.py GE_e4_merge_AGNseed001_SimBKG_CLUseed001 > logs/summary_sky_map_GE_e4_merge_AGNseed001_SimBKG_CLUseed001.log & # TODO
nohup python make_summary_skymap.py GE_e4_merge_SimBKG_CLUseed001            > logs/summary_sky_map_GE_e4_merge_SimBKG_CLUseed001.log            & # TODO

nohup python make_summary_skymap.py GE_e4_merge_AGNseed002_SimBKG            > logs/summary_sky_map_GE_e4_merge_AGNseed002_SimBKG.log            & # TODO
nohup python make_summary_skymap.py GE_e4_merge_AGNseed002_SimBKG_CLUseed002 > logs/summary_sky_map_GE_e4_merge_AGNseed002_SimBKG_CLUseed002.log & # TODO
nohup python make_summary_skymap.py GE_e4_merge_SimBKG_CLUseed002            > logs/summary_sky_map_GE_e4_merge_SimBKG_CLUseed002.log            & # TODO

nohup python make_summary_skymap.py GE_e4_merge_AGNseed003_SimBKG            > logs/summary_sky_map_GE_e4_merge_AGNseed003_SimBKG.log            & # TODO
nohup python make_summary_skymap.py GE_e4_merge_AGNseed003_SimBKG_CLUseed003 > logs/summary_sky_map_GE_e4_merge_AGNseed003_SimBKG_CLUseed003.log & # TODO
nohup python make_summary_skymap.py GE_e4_merge_SimBKG_CLUseed003            > logs/summary_sky_map_GE_e4_merge_SimBKG_CLUseed003.log            & # TODO

nohup python make_summary_skymap.py GE_e4_merge_AGNseed004_SimBKG            > logs/summary_skymap_GE_e4_merge_AGNseed004_SimBKG.log            & # TODO
nohup python make_summary_skymap.py GE_e4_merge_AGNseed004_SimBKG_CLUseed004 > logs/summary_skymap_GE_e4_merge_AGNseed004_SimBKG_CLUseed004.log & # TODO
nohup python make_summary_skymap.py GE_e4_merge_SimBKG_CLUseed004            > logs/summary_skymap_GE_e4_merge_SimBKG_CLUseed004.log            & # TODO

nohup python make_summary_skymap.py GE_e4_merge_AGNseed005_SimBKG            > logs/summary_skymap_GE_e4_merge_AGNseed005_SimBKG.log            & # TODO
nohup python make_summary_skymap.py GE_e4_merge_AGNseed005_SimBKG_CLUseed005 > logs/summary_skymap_GE_e4_merge_AGNseed005_SimBKG_CLUseed005.log & # TODO
nohup python make_summary_skymap.py GE_e4_merge_SimBKG_CLUseed005            > logs/summary_skymap_GE_e4_merge_SimBKG_CLUseed005.log            & # TODO

nohup python make_summary_skymap.py GE_e4_merge_AGNseed006_SimBKG            > logs/summary_skymap_GE_e4_merge_AGNseed006_SimBKG.log            & # TODO
nohup python make_summary_skymap.py GE_e4_merge_AGNseed006_SimBKG_CLUseed006 > logs/summary_skymap_GE_e4_merge_AGNseed006_SimBKG_CLUseed006.log & # TODO
nohup python make_summary_skymap.py GE_e4_merge_SimBKG_CLUseed006            > logs/summary_skymap_GE_e4_merge_SimBKG_CLUseed006.log            & # TODO

nohup python make_summary_skymap.py GE_e4_merge_AGNseed007_SimBKG            > logs/summary_skymap_GE_e4_merge_AGNseed007_SimBKG.log            & # TODO
nohup python make_summary_skymap.py GE_e4_merge_AGNseed007_SimBKG_CLUseed007 > logs/summary_skymap_GE_e4_merge_AGNseed007_SimBKG_CLUseed007.log & # TODO
nohup python make_summary_skymap.py GE_e4_merge_SimBKG_CLUseed007            > logs/summary_skymap_GE_e4_merge_SimBKG_CLUseed007.log            & # TODO

nohup python make_summary_skymap.py GE_e4_merge_AGNseed008_SimBKG            > logs/summary_skymap_GE_e4_merge_AGNseed008_SimBKG.log            & # TODO
nohup python make_summary_skymap.py GE_e4_merge_AGNseed008_SimBKG_CLUseed008 > logs/summary_skymap_GE_e4_merge_AGNseed008_SimBKG_CLUseed008.log & # TODO
nohup python make_summary_skymap.py GE_e4_merge_SimBKG_CLUseed008            > logs/summary_skymap_GE_e4_merge_SimBKG_CLUseed008.log            & # TODO


# write all commands for the list of folders above.
python write_exec_loop.py # > exec.sh # TODO

export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/esass/runs

# execute all commands of interest in an eSASS loaded container

nohup sh GE_e4_merge_AGNseed002_SimBKG_CLUseed002_processing_0000.sh > logs/AGNseed002_SimBKG_CLUseed002_processing_0000.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_CLUseed002_processing_0050.sh > logs/AGNseed002_SimBKG_CLUseed002_processing_0050.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_CLUseed002_processing_0100.sh > logs/AGNseed002_SimBKG_CLUseed002_processing_0100.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_CLUseed002_processing_0150.sh > logs/AGNseed002_SimBKG_CLUseed002_processing_0150.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_CLUseed002_processing_0200.sh > logs/AGNseed002_SimBKG_CLUseed002_processing_0200.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_CLUseed002_processing_0250.sh > logs/AGNseed002_SimBKG_CLUseed002_processing_0250.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_CLUseed002_processing_0300.sh > logs/AGNseed002_SimBKG_CLUseed002_processing_0300.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_CLUseed002_processing_0350.sh > logs/AGNseed002_SimBKG_CLUseed002_processing_0350.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_CLUseed002_processing_0400.sh > logs/AGNseed002_SimBKG_CLUseed002_processing_0400.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_CLUseed002_processing_0450.sh > logs/AGNseed002_SimBKG_CLUseed002_processing_0450.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_CLUseed002_processing_0500.sh > logs/AGNseed002_SimBKG_CLUseed002_processing_0500.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_CLUseed002_processing_0550.sh > logs/AGNseed002_SimBKG_CLUseed002_processing_0550.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_CLUseed002_processing_0600.sh > logs/AGNseed002_SimBKG_CLUseed002_processing_0600.log & # ONGOING
nohup sh GE_e4_merge_AGNseed002_SimBKG_CLUseed002_processing_0650.sh > logs/AGNseed002_SimBKG_CLUseed002_processing_0650.log & # ONGOING

nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0000.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0000.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0050.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0050.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0100.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0100.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0150.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0150.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0200.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0200.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0250.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0250.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0300.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0300.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0350.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0350.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0400.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0400.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0450.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0450.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0500.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0500.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0550.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0550.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0600.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0600.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0650.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0650.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0700.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0700.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0750.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0750.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0800.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0800.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0850.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0850.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0900.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0900.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_0950.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_0950.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_1000.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_1000.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_1050.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_1050.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_1100.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_1100.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_1150.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_1150.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_1200.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_1200.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_1250.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_1250.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_1300.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_1300.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_1350.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_1350.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_1400.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_1400.log & # ONGOING
nohup sh GE_e4_merge_AGNseed003_SimBKG_CLUseed003_processing_1450.sh > logs/AGNseed003_SimBKG_CLUseed003_processing_1450.log & # ONGOING

export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/esass/runs

nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0000.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0000.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0050.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0050.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0100.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0100.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0150.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0150.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0200.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0200.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0250.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0250.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0300.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0300.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0350.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0350.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0400.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0400.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0450.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0450.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0500.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0500.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0550.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0550.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0600.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0600.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0650.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0650.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0700.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0700.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0750.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0750.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0800.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0800.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0850.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0850.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0900.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0900.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_0950.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_0950.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_1000.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_1000.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_1050.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_1050.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_1100.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_1100.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_1150.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_1150.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_1200.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_1200.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_1250.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_1250.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_1300.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_1300.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_1350.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_1350.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_1400.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_1400.log & # ONGOING
nohup sh GE_e4_merge_AGNseed004_SimBKG_CLUseed004_processing_1450.sh > logs/AGNseed004_SimBKG_CLUseed004_processing_1450.log & # ONGOING


export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/esass/runs

nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0000.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0000.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0050.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0050.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0100.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0100.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0150.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0150.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0200.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0200.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0250.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0250.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0300.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0300.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0350.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0350.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0400.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0400.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0450.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0450.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0500.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0500.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0550.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0550.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0600.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0600.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0650.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0650.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0700.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0700.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0750.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0750.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0800.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0800.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0850.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0850.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0900.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0900.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_0950.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_0950.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_1000.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_1000.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_1050.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_1050.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_1100.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_1100.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_1150.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_1150.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_1200.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_1200.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_1250.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_1250.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_1300.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_1300.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_1350.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_1350.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_1400.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_1400.log & # ONGOING
nohup sh GE_e4_merge_AGNseed005_SimBKG_CLUseed005_processing_1450.sh > logs/AGNseed005_SimBKG_CLUseed005_processing_1450.log & # ONGOING


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

python summary_plotting_RS.py GE_e4_merge_AGNseed001_SimBKG_CLUseed001 /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_AGNseed002_SimBKG_CLUseed002 /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_AGNseed003_SimBKG_CLUseed003 /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_AGNseed004_SimBKG_CLUseed004 /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins

python summary_plotting_RS.py GE_e4_merge_SimBKG_CLUseed001 /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins
python summary_plotting_RS.py GE_e4_merge_AGNseed001_SimBKG /home/idies/workspace/erosim/software/st_mod_data/data/validation/xray_twins














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

