#!/bin/bash/
export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/esass/

BG model files are here :
~/workspace/erosim/simput/bkg_erosita_simput_full_sky
catalog.fits has the full map

# in any container with python
nohup python merge_events_onlyBG.py  > logs/merge_events_onlyBG.log  & # RUNNING job

nohup python merge_events_e5.py 1 1     > logs/merge_events_1_1.log  & # TODO
nohup python merge_events_noCLU.py 1 > logs/merge_events_noCLU_1.log & # TODO
nohup python merge_events_noAGN.py 1 > logs/merge_events_noAGN_1.log & # TODO

nohup python merge_events_e5.py 2 2     > logs/merge_events_2_2.log  & # TODO
nohup python merge_events_noCLU.py 2 > logs/merge_events_noCLU_2.log & # TODO
nohup python merge_events_noAGN.py 2 > logs/merge_events_noAGN_2.log & # TODO

nohup python merge_events_e5.py 3 3     > logs/merge_events_3_3.log  & # TODO
nohup python merge_events_noCLU.py 3 > logs/merge_events_noCLU_3.log & # TODO
nohup python merge_events_noAGN.py 3 > logs/merge_events_noAGN_3.log & # TODO

export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/esass/
python merge_events_onlyBG.py

export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/esass/
python merge_events_e5.py 1 1

export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/esass/
python merge_events_noCLU.py 1

export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/esass/
python merge_events_noAGN.py 1

export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/esass/
python merge_events_e5.py 2 2

export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/esass/
python merge_events_noCLU.py 2

export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/esass/
python merge_events_noAGN.py 2

export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/esass/
python merge_events_e5.py 3 3

export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/esass/
python merge_events_noAGN.py 3

export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/esass/
python merge_events_noAGN.py 3



GE_e4_merge_AGNseed001_SimBKG : AGN + BKG
GE_e4_merge_AGNseed001_SimBKG_CLUseed001 : AGN + BKG + Cluster
GE_e4_merge_SimBKG : only BKG
GE_e4_merge_SimBKG_CLUseed001 : BKG + Cluster

# write eSASS commands in the relevant folders
nohup python eRASSX_write_scripts.py GE_e4_merge_SimBKG                       > logs/eRASSX_write_scripts_GE_e4_merge_SimBKG.log                       & # TODO

nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed001_SimBKG            > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed001_SimBKG.log            & # TODO
nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed001_SimBKG_CLUseed001 > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed001_SimBKG_CLUseed001.log & # TODO
nohup python eRASSX_write_scripts.py GE_e4_merge_SimBKG_CLUseed001            > logs/eRASSX_write_scripts_GE_e4_merge_SimBKG_CLUseed001.log            & # TODO

nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed002_SimBKG            > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed002_SimBKG.log            & # TODO
nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed002_SimBKG_CLUseed002 > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed002_SimBKG_CLUseed002.log & # TODO
nohup python eRASSX_write_scripts.py GE_e4_merge_SimBKG_CLUseed002            > logs/eRASSX_write_scripts_GE_e4_merge_SimBKG_CLUseed002.log            & # TODO

nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed003_SimBKG            > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed003_SimBKG.log            & # TODO
nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed003_SimBKG_CLUseed003 > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed003_SimBKG_CLUseed003.log & # TODO
nohup python eRASSX_write_scripts.py GE_e4_merge_SimBKG_CLUseed003            > logs/eRASSX_write_scripts_GE_e4_merge_SimBKG_CLUseed003.log            & # TODO



# monitor progress

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


# write all commands for the list of folders above.
python write_exec_loop.py # > exec.sh


export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/esass/
cd $GIT_STMOD/src/esass/runs
sh GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_0000.sh  > logs/GE_e4_merge_AGNseed001_SimBKG_CLUseed001_processing_0000.log & # running
sh GE_e4_merge_SimBKG_CLUseed001_processing_0000.sh             > logs/GE_e4_merge_SimBKG_CLUseed001_processing_0000.log            & # running
sh GE_e4_merge_SimBKG_processing_0000.sh                        > logs/GE_e4_merge_SimBKG_processing_0000.log                       & # running



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

