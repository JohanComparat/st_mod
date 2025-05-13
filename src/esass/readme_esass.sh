#!/bin/bash/
export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/esass/

# in the container with sixte 2.7
nohup python merge_events.py 1 1     > logs/merge_events_1_1.log & # ONGOING
nohup python merge_events_noCLU.py 1 > logs/merge_events_noCLU_1.log & # ONGOING
nohup python merge_events_noAGN.py 1 > logs/merge_events_noAGN_1.log & # ONGOING
nohup python merge_events_onlyBG.py  > logs/merge_events_onlyBG.log & # ONGOING

(base) idies@b631d6bc716a:~/workspace/erosim/software/st_mod/src/esass$ tail logs/merge_events_1_1.log
0.019 0.0183
SRVMAP OWNER   RA_MIN     RA_MAX     DE_MIN     DE_MAX     RA_CEN     DE_CEN    ELON_CEN   ELAT_CEN   GLON_CEN   GLAT_CEN    X_MIN      Y_MIN    N_NBRS FIELD1 FIELD2 FIELD3 FIELD4 FIELD5 FIELD6 FIELD7 FIELD8 FIELD9
                deg        deg        deg        deg        deg        deg        deg        deg        deg        deg       arcmin     arcmin
------ ----- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ------ ------ ------ ------ ------ ------ ------ ------ ------ ------
176054     0 174.545455 178.181818  34.500000  37.500000 176.363636  36.006729 160.813929  31.262582 176.504078  73.412500  18.110413  17.606671      7 172051 176051 180051 173054 180054 175057 178057      0      0
Traceback (most recent call last):
  File "/home/idies/workspace/erosim/software/st_mod/src/esass/merge_events.py", line 137, in <module>
    id_B = np.random.choice(np.arange(len(bg_tm)), size = N_ev_B, replace = False)
  File "mtrand.pyx", line 984, in numpy.random.mtrand.RandomState.choice
ValueError: Cannot take a larger sample than population when 'replace=False'


GE_e4_merge_AGNseed001_SimBKG : AGN + BKG
GE_e4_merge_AGNseed001_SimBKG_CLUseed001 : AGN + BKG + Cluster
GE_e4_merge_SimBKG : only BKG
GE_e4_merge_SimBKG_CLUseed001 : BKG + Cluster

# write eSASS commands
nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed001_SimBKG            > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed001_SimBKG.log            &
nohup python eRASSX_write_scripts.py GE_e4_merge_AGNseed001_SimBKG_CLUseed001 > logs/eRASSX_write_scripts_GE_e4_merge_AGNseed001_SimBKG_CLUseed001.log &
nohup python eRASSX_write_scripts.py GE_e4_merge_SimBKG                       > logs/eRASSX_write_scripts_GE_e4_merge_SimBKG.log                       &
nohup python eRASSX_write_scripts.py GE_e4_merge_SimBKG_CLUseed001            > logs/eRASSX_write_scripts_GE_e4_merge_SimBKG_CLUseed001.log            &

# execute eSASS commands on test field:
cd /home/idies/workspace/erosim/Uchuu/LCerass/121048/GE_e4_merge_AGNseed001_SimBKG/eSASS
sh 121048_pipeline_img1.sh # DONE
sh 121048_pipeline_det1.sh
sh 121048_pipeline_Src1.sh

cd /home/idies/workspace/erosim/Uchuu/LCerass/121048/GE_e4_merge_AGNseed001_SimBKG_CLUseed001/eSASS
sh 121048_pipeline_img1.sh # DONE
sh 121048_pipeline_det1.sh
sh 121048_pipeline_Src1.sh

cd /home/idies/workspace/erosim/Uchuu/LCerass/121048/GE_e4_merge_SimBKG/eSASS
sh 121048_pipeline_img1.sh # DONE
sh 121048_pipeline_det1.sh
sh 121048_pipeline_Src1.sh

cd /home/idies/workspace/erosim/Uchuu/LCerass/121048/GE_e4_merge_SimBKG_CLUseed001/eSASS
sh 121048_pipeline_img1.sh # DONE
sh 121048_pipeline_det1.sh
sh 121048_pipeline_Src1.sh
