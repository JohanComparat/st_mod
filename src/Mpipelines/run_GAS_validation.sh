#!/bin/bash
conda activate stmod
cd $GIT_STMOD/src/Mpipelines

# nohup python GAS_validation_logNlogS.py LC0002    > log.log &
# nohup python GAS_validation_logNlogS.py LC0060    > log.log &
# nohup python GAS_validation_logNlogS.py LC1800 8  > log.log &
nohup python GAS_validation_logNlogS.py LC1800 12 > log.log &
nohup python GAS_validation_logNlogS.py LC1800 17 > log.log &
nohup python GAS_validation_logNlogS.py LC1800 22 > log.log &
nohup python GAS_validation_logNlogS.py LC1800 27 > log.log &
# nohup python GAS_validation_logNlogS.py FullSky 5 > log.log &
nohup python GAS_validation_logNlogS.py FullSky 6 > log.log &
nohup python GAS_validation_logNlogS.py FullSky 7 > log.log &
nohup python GAS_validation_logNlogS.py FullSky 8 > log.log &

git add ../../data/validation/validation_GAS/logNlogS
git commit -m"logNlogS"
git push
[1]   Done                    nohup python GAS_validation_logNlogS.py LC0002 > log.log
[2]   Done                    nohup python GAS_validation_logNlogS.py LC0060 > log.log
[3]   Done                    nohup python GAS_validation_logNlogS.py LC1800 8 > log.log
[4]   Exit 1                  nohup python GAS_validation_logNlogS.py LC1800 12 > log.log
[5]   Exit 1                  nohup python GAS_validation_logNlogS.py LC1800 17 > log.log
[6]   Exit 1                  nohup python GAS_validation_logNlogS.py LC1800 22 > log.log
[7]   Exit 1                  nohup python GAS_validation_logNlogS.py LC1800 27 > log.log
[8]   Done                    nohup python GAS_validation_logNlogS.py FullSky 5 > log.log
[9]   Exit 1                  nohup python GAS_validation_logNlogS.py FullSky 6 > log.log
[10]-  Exit 1                  nohup python GAS_validation_logNlogS.py FullSky 7 > log.log
[11]+  Exit 1                  nohup python GAS_validation_logNlogS.py FullSky 8 > log.log

# on the full sky

# python GAS_pipeline.py z0p00 FullSky
# python GAS_pipeline.py z0p02 FullSky
# python GAS_validation_ScalingRelation.py z0p00 FullSky
# python GAS_validation_ScalingRelation.py z0p02 FullSky
# python GAS_pipeline.py z0p05 FullSky
# python GAS_validation_ScalingRelation.py z0p05 FullSky
# python GAS_pipeline.py z0p09 FullSky
# python GAS_validation_ScalingRelation.py z0p09 FullSky
#
# python GAS_pipeline.py z0p00 FullSky
# python GAS_pipeline.py z0p02 FullSky
python GAS_validation_ScalingRelation.py z0p00 FullSky
python GAS_validation_ScalingRelation.py z0p02 FullSky
python GAS_validation_ScalingRelation.py z0p05 FullSky
python GAS_validation_ScalingRelation.py z0p09 FullSky
python GAS_validation_ScalingRelation.py z0p14 FullSky
python GAS_validation_ScalingRelation.py z0p19 FullSky
python GAS_validation_ScalingRelation.py z0p25 FullSky
python GAS_validation_ScalingRelation.py z0p30 FullSky
python GAS_validation_ScalingRelation.py z0p36 FullSky
python GAS_validation_ScalingRelation.py z0p43 FullSky
python GAS_validation_ScalingRelation.py z0p49 FullSky
python GAS_validation_ScalingRelation.py z0p56 FullSky
python GAS_validation_ScalingRelation.py z0p63 FullSky
python GAS_validation_ScalingRelation.py z0p70 FullSky
python GAS_validation_ScalingRelation.py z0p78 FullSky

#



python GAS_validation_ScalingRelation.py z0p09 LC1800
python GAS_validation_ScalingRelation.py z0p14 LC1800
python GAS_validation_ScalingRelation.py z0p19 LC1800
python GAS_validation_ScalingRelation.py z0p25 LC1800
python GAS_validation_ScalingRelation.py z0p30 LC1800
python GAS_validation_ScalingRelation.py z0p36 LC1800
python GAS_validation_ScalingRelation.py z0p43 LC1800
python GAS_validation_ScalingRelation.py z0p49 LC1800
python GAS_validation_ScalingRelation.py z0p56 LC1800
python GAS_validation_ScalingRelation.py z0p63 LC1800
python GAS_validation_ScalingRelation.py z0p70 LC1800
python GAS_validation_ScalingRelation.py z0p78 LC1800
python GAS_validation_ScalingRelation.py z0p86 LC1800
python GAS_validation_ScalingRelation.py z0p94 LC1800

python GAS_validation_ScalingRelation.py z1p03 LC0060
python GAS_validation_ScalingRelation.py z1p12 LC0060
python GAS_validation_ScalingRelation.py z1p22 LC0060
python GAS_validation_ScalingRelation.py z1p32 LC0060
python GAS_validation_ScalingRelation.py z1p43 LC0060
python GAS_validation_ScalingRelation.py z1p54 LC0060
python GAS_validation_ScalingRelation.py z1p65 LC0060
python GAS_validation_ScalingRelation.py z1p77 LC0060
python GAS_validation_ScalingRelation.py z1p90 LC0060
python GAS_validation_ScalingRelation.py z2p03 LC0060
python GAS_validation_ScalingRelation.py z2p17 LC0060
python GAS_validation_ScalingRelation.py z2p31 LC0060
python GAS_validation_ScalingRelation.py z2p46 LC0060
python GAS_validation_ScalingRelation.py z2p62 LC0060
python GAS_validation_ScalingRelation.py z2p78 LC0060
python GAS_validation_ScalingRelation.py z2p95 LC0060
python GAS_validation_ScalingRelation.py z3p13 LC0060
python GAS_validation_ScalingRelation.py z3p32 LC0060
python GAS_validation_ScalingRelation.py z3p61 LC0060
python GAS_validation_ScalingRelation.py z3p93 LC0060

python GAS_validation_ScalingRelation.py z4p27 LC0002
python GAS_validation_ScalingRelation.py z4p63 LC0002
python GAS_validation_ScalingRelation.py z5p15 LC0002
python GAS_validation_ScalingRelation.py z5p73 LC0002



python GAS_validation_ScalingRelation.py z0p00 LC0060
python GAS_validation_ScalingRelation.py z0p02 LC0060
python GAS_validation_ScalingRelation.py z0p05 LC0060
python GAS_validation_ScalingRelation.py z0p09 LC0060
python GAS_validation_ScalingRelation.py z0p14 LC0060
python GAS_validation_ScalingRelation.py z0p19 LC0060
python GAS_validation_ScalingRelation.py z0p25 LC0060
python GAS_validation_ScalingRelation.py z0p30 LC0060
python GAS_validation_ScalingRelation.py z0p36 LC0060
python GAS_validation_ScalingRelation.py z0p43 LC0060
python GAS_validation_ScalingRelation.py z0p49 LC0060
python GAS_validation_ScalingRelation.py z0p56 LC0060
python GAS_validation_ScalingRelation.py z0p63 LC0060
python GAS_validation_ScalingRelation.py z0p70 LC0060
python GAS_validation_ScalingRelation.py z0p78 LC0060
python GAS_validation_ScalingRelation.py z0p86 LC0060
python GAS_validation_ScalingRelation.py z0p94 LC0060
python GAS_validation_ScalingRelation.py z1p03 LC0060
python GAS_validation_ScalingRelation.py z1p12 LC0060
python GAS_validation_ScalingRelation.py z1p22 LC0060
python GAS_validation_ScalingRelation.py z1p32 LC0060
python GAS_validation_ScalingRelation.py z1p43 LC0060
python GAS_validation_ScalingRelation.py z1p54 LC0060
python GAS_validation_ScalingRelation.py z1p65 LC0060
python GAS_validation_ScalingRelation.py z1p77 LC0060
python GAS_validation_ScalingRelation.py z1p90 LC0060
python GAS_validation_ScalingRelation.py z2p03 LC0060
python GAS_validation_ScalingRelation.py z2p17 LC0060
python GAS_validation_ScalingRelation.py z2p31 LC0060
python GAS_validation_ScalingRelation.py z2p46 LC0060
python GAS_validation_ScalingRelation.py z2p62 LC0060
python GAS_validation_ScalingRelation.py z2p78 LC0060
python GAS_validation_ScalingRelation.py z2p95 LC0060
python GAS_validation_ScalingRelation.py z3p13 LC0060
python GAS_validation_ScalingRelation.py z3p32 LC0060
python GAS_validation_ScalingRelation.py z3p61 LC0060
python GAS_validation_ScalingRelation.py z3p93 LC0060
