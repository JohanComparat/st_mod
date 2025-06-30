#!/bin/bash

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
