#!/bin/bash/
export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/sixte_GAS/

# in the container with sixte 2.7
python merge_events.py

# in the container with sixte 2.7
# python calibrate_events_UCHUU.py # ONGOING for eRO_DE AGN + CLU
# in the container with esass
# python calibrate_events_UCHUU_part2.py

# write eSASS commands
python eRASSX_write_scripts.py

# execute eSASS commands
cd /home/idies/workspace/erosim/Uchuu/LCerass/121048/sim_evt_e4_merge/eSASS
sh 121048_pipeline_img1.sh # DONE
sh 121048_pipeline_det1.sh
sh 121048_pipeline_Src1.sh
