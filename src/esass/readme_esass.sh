#!/bin/bash/
export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/esass/

# in the container with sixte 2.7
nohup python merge_events.py 1 1 > logs/merge_events_1_1.log & # ONGOING
nohup python merge_events_noCLU.py 1 > logs/merge_events_noCLU_1.log & # ONGOING
nohup python merge_events_noAGN.py 1 > logs/merge_events_noAGN_1.log & # ONGOING
nohup python merge_events_onlyBG.py> logs/merge_events_onlyBG.log & # ONGOING

GE_e4_merge_AGNseed001_SimBKG
GE_e4_merge_AGNseed001_SimBKG_CLUseed001

# write eSASS commands
python eRASSX_write_scripts.py GE_e4_merge_AGNseed001_SimBKG #
python eRASSX_write_scripts.py GE_e4_merge_AGNseed001_SimBKG_CLUseed001 #

# execute eSASS commands on test field:
cd /home/idies/workspace/erosim/Uchuu/LCerass/121048/GE_e4_merge_AGNseed001_SimBKG/eSASS
sh 121048_pipeline_img1.sh # DONE
sh 121048_pipeline_det1.sh
sh 121048_pipeline_Src1.sh

cd /home/idies/workspace/erosim/Uchuu/LCerass/121048/GE_e4_merge_AGNseed001_SimBKG_CLUseed001eSASS
sh 121048_pipeline_img1.sh # DONE
sh 121048_pipeline_det1.sh
sh 121048_pipeline_Src1.sh
