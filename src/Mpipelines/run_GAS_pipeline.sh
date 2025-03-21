#!/bin/bash

export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'

# conda activate stmod
cd $GIT_STMOD/src/Mpipelines

cd ~/workspace/erosim/sixte_output_events
rsync -avz joco@raven.mpcdf.mpg.de:/ptmp/joco/erosim/sixte_output_events/* .

#
# tabulates a file with profiles inside
#
# python GAS_setup.py z0p14 FullSky # DONE only run once and do not rerun.

python GAS_pipeline.py z0p00 FullSky # DONE
python GAS_pipeline.py z0p02 FullSky # DONE
python GAS_pipeline.py z0p05 FullSky # DONE
python GAS_pipeline.py z0p09 FullSky # DONE
python GAS_pipeline.py z0p14 FullSky # DONE
python GAS_pipeline.py z0p19 FullSky # DONE
python GAS_pipeline.py z0p25 FullSky # DONE
python GAS_pipeline.py z0p30 FullSky # DONE
python GAS_pipeline.py z0p36 FullSky # DONE
python GAS_pipeline.py z0p43 FullSky # DONE
python GAS_pipeline.py z0p49 FullSky # DONE
python GAS_pipeline.py z0p56 FullSky # DONE
python GAS_pipeline.py z0p63 FullSky # DONE
python GAS_pipeline.py z0p70 FullSky # DONE
python GAS_pipeline.py z0p78 FullSky # DONE
python GAS_pipeline.py z0p86 FullSky # DONE
python GAS_pipeline.py z0p94 FullSky # DONE
python GAS_pipeline.py z1p03 FullSky # DONE
python GAS_pipeline.py z1p12 FullSky # DONE
python GAS_pipeline.py z1p22 FullSky # DONE
python GAS_pipeline.py z1p32 FullSky # DONE
python GAS_pipeline.py z1p43 FullSky # DONE
python GAS_pipeline.py z1p54 FullSky # DONE

#
# REST is NOT NEEDED
#
#
python GAS_pipeline.py z0p00 LC0002 # DONE
python GAS_pipeline.py z0p02 LC0002 # DONE
python GAS_pipeline.py z0p05 LC0002 # DONE
python GAS_pipeline.py z0p09 LC0002 # DONE
python GAS_pipeline.py z0p14 LC0002 # DONE
python GAS_pipeline.py z0p19 LC0002 # DONE
python GAS_pipeline.py z0p25 LC0002 # DONE
python GAS_pipeline.py z0p30 LC0002 # DONE
python GAS_pipeline.py z0p36 LC0002 # DONE
python GAS_pipeline.py z0p43 LC0002 # DONE
python GAS_pipeline.py z0p49 LC0002 # DONE
python GAS_pipeline.py z0p56 LC0002 # DONE
python GAS_pipeline.py z0p63 LC0002 # DONE
python GAS_pipeline.py z0p70 LC0002 # DONE
python GAS_pipeline.py z0p78 LC0002 # DONE
python GAS_pipeline.py z0p86 LC0002 # DONE
python GAS_pipeline.py z0p94 LC0002 # DONE
python GAS_pipeline.py z1p03 LC0002 # DONE
python GAS_pipeline.py z1p12 LC0002 # DONE
python GAS_pipeline.py z1p22 LC0002 # DONE
python GAS_pipeline.py z1p32 LC0002 # DONE
python GAS_pipeline.py z1p43 LC0002 # DONE
python GAS_pipeline.py z1p54 LC0002 # DONE
python GAS_pipeline.py z1p65 LC0002 # DONE
python GAS_pipeline.py z1p77 LC0002 # DONE
python GAS_pipeline.py z1p90 LC0002 # DONE
python GAS_pipeline.py z2p03 LC0002 # DONE
python GAS_pipeline.py z2p17 LC0002 # DONE
python GAS_pipeline.py z2p31 LC0002 # DONE
python GAS_pipeline.py z2p46 LC0002 # DONE
python GAS_pipeline.py z2p62 LC0002 # DONE
python GAS_pipeline.py z2p78 LC0002 # DONE
python GAS_pipeline.py z2p95 LC0002 # DONE
python GAS_pipeline.py z3p13 LC0002 # DONE
python GAS_pipeline.py z3p32 LC0002 # DONE
python GAS_pipeline.py z3p61 LC0002 # DONE
python GAS_pipeline.py z3p93 LC0002 # DONE
python GAS_pipeline.py z4p27 LC0002 # DONE
python GAS_pipeline.py z4p63 LC0002 # DONE
python GAS_pipeline.py z5p15 LC0002 # DONE
python GAS_pipeline.py z5p73 LC0002 # DONE

python GAS_pipeline.py z0p00 LC0060 # DONE
python GAS_pipeline.py z0p02 LC0060 # DONE
python GAS_pipeline.py z0p05 LC0060 # DONE
python GAS_pipeline.py z0p09 LC0060 # DONE
python GAS_pipeline.py z0p14 LC0060 # DONE
python GAS_pipeline.py z0p19 LC0060 # DONE
python GAS_pipeline.py z0p25 LC0060 # DONE
python GAS_pipeline.py z0p30 LC0060 # DONE
python GAS_pipeline.py z0p36 LC0060 # DONE
python GAS_pipeline.py z0p43 LC0060 # DONE
python GAS_pipeline.py z0p49 LC0060 # DONE
python GAS_pipeline.py z0p56 LC0060 # DONE
python GAS_pipeline.py z0p63 LC0060 # DONE
python GAS_pipeline.py z0p70 LC0060 # DONE
python GAS_pipeline.py z0p78 LC0060 # DONE
python GAS_pipeline.py z0p86 LC0060 # DONE
python GAS_pipeline.py z0p94 LC0060 # DONE
python GAS_pipeline.py z1p03 LC0060 # DONE
python GAS_pipeline.py z1p12 LC0060 # DONE
python GAS_pipeline.py z1p22 LC0060 # DONE
python GAS_pipeline.py z1p32 LC0060 # DONE
python GAS_pipeline.py z1p43 LC0060 # DONE
python GAS_pipeline.py z1p54 LC0060 # DONE
python GAS_pipeline.py z1p65 LC0060 # DONE
python GAS_pipeline.py z1p77 LC0060 # DONE
python GAS_pipeline.py z1p90 LC0060 # DONE
python GAS_pipeline.py z2p03 LC0060 # DONE
python GAS_pipeline.py z2p17 LC0060 # DONE
python GAS_pipeline.py z2p31 LC0060 # DONE
python GAS_pipeline.py z2p46 LC0060 # DONE
python GAS_pipeline.py z2p62 LC0060 # DONE
python GAS_pipeline.py z2p78 LC0060 # DONE
python GAS_pipeline.py z2p95 LC0060 # DONE
python GAS_pipeline.py z3p13 LC0060 # DONE
python GAS_pipeline.py z3p32 LC0060 # DONE
python GAS_pipeline.py z3p61 LC0060 # DONE
python GAS_pipeline.py z3p93 LC0060 # DONE

python GAS_pipeline.py z0p00 LC1800 # DONE
python GAS_pipeline.py z0p02 LC1800 # DONE
python GAS_pipeline.py z0p05 LC1800 # DONE
python GAS_pipeline.py z0p09 LC1800 # DONE
python GAS_pipeline.py z0p14 LC1800 # DONE
python GAS_pipeline.py z0p19 LC1800 # DONE
python GAS_pipeline.py z0p25 LC1800 # DONE
python GAS_pipeline.py z0p30 LC1800 # DONE
python GAS_pipeline.py z0p36 LC1800 # DONE
python GAS_pipeline.py z0p43 LC1800 # DONE
python GAS_pipeline.py z0p49 LC1800 # DONE
python GAS_pipeline.py z0p56 LC1800 # DONE
python GAS_pipeline.py z0p63 LC1800 # DONE
python GAS_pipeline.py z0p70 LC1800 # DONE
python GAS_pipeline.py z0p78 LC1800 # DONE
python GAS_pipeline.py z0p86 LC1800 # DONE
python GAS_pipeline.py z0p94 LC1800 # DONE
python GAS_pipeline.py z1p03 LC1800 # DONE
python GAS_pipeline.py z1p12 LC1800 # DONE
python GAS_pipeline.py z1p22 LC1800 # DONE
python GAS_pipeline.py z1p32 LC1800 # DONE
python GAS_pipeline.py z1p43 LC1800 # DONE
python GAS_pipeline.py z1p54 LC1800 # DONE
python GAS_pipeline.py z1p65 LC1800 # DONE
python GAS_pipeline.py z1p77 LC1800 # DONE
python GAS_pipeline.py z1p90 LC1800 # DONE
python GAS_pipeline.py z2p03 LC1800 # DONE
