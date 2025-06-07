#!/bin/bash
conda activate stmod
cd $GIT_STMOD/src/Mpipelines

export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'

conda activate stmod
cd $GIT_STMOD/src/Mpipelines
conda activate clustering

nohup python GasGal_validation_XCORR.py > GasGal_validation_XCORR.log &

nohup python GAS_validation_WPRP.py z0p14 > logs/wprpGASz0p14.log & # DONE
nohup python GAS_validation_WPRP.py z0p19 > logs/wprpGASz0p19.log & # DONE
nohup python GAS_validation_WPRP.py z0p25 > logs/wprpGASz0p25.log & # DONE

python GAS_validation_ScalingRelation.py z0p00 FullSky # DONE
python GAS_validation_ScalingRelation.py z0p02 FullSky # DONE
python GAS_validation_ScalingRelation.py z0p05 FullSky # DONE
python GAS_validation_ScalingRelation.py z0p09 FullSky # DONE
python GAS_validation_ScalingRelation.py z0p14 FullSky # DONE
python GAS_validation_ScalingRelation.py z0p19 FullSky # DONE
python GAS_validation_ScalingRelation.py z0p25 FullSky # DONE
python GAS_validation_ScalingRelation.py z0p30 FullSky # DONE
python GAS_validation_ScalingRelation.py z0p36 FullSky # DONE
python GAS_validation_ScalingRelation.py z0p43 FullSky # DONE
python GAS_validation_ScalingRelation.py z0p49 FullSky # DONE
python GAS_validation_ScalingRelation.py z0p56 FullSky # DONE
python GAS_validation_ScalingRelation.py z0p63 FullSky # DONE
python GAS_validation_ScalingRelation.py z0p70 FullSky # DONE
python GAS_validation_ScalingRelation.py z0p78 FullSky # DONE
python GAS_validation_ScalingRelation.py z0p86 FullSky # DONE
python GAS_validation_ScalingRelation.py z0p94 FullSky # DONE
python GAS_validation_ScalingRelation.py z1p03 FullSky # DONE
python GAS_validation_ScalingRelation.py z1p12 FullSky # DONE
python GAS_validation_ScalingRelation.py z1p22 FullSky # DONE
python GAS_validation_ScalingRelation.py z1p32 FullSky # DONE
python GAS_validation_ScalingRelation.py z1p43 FullSky # DONE
python GAS_validation_ScalingRelation.py z1p54 FullSky # DONE

python GAS_validation_logNlogS.py LC0002 # DONE
python GAS_validation_logNlogS.py LC0060 # DONE
python GAS_validation_logNlogS.py LC1800 8  # DONE
python GAS_validation_logNlogS.py LC1800 12 # DONE
python GAS_validation_logNlogS.py LC1800 17 # DONE
python GAS_validation_logNlogS.py LC1800 22 # DONE
python GAS_validation_logNlogS.py LC1800 27 # DONE
python GAS_validation_logNlogS.py FullSky 3 # DONE
python GAS_validation_logNlogS.py FullSky 5 # DONE
python GAS_validation_logNlogS.py FullSky 8 # DONE
python GAS_validation_logNlogS.py FullSky 12 # DONE

python GAS_validation_ScalingRelation.py z0p09 LC1800 # DONE
python GAS_validation_ScalingRelation.py z0p14 LC1800 # DONE
python GAS_validation_ScalingRelation.py z0p19 LC1800 # DONE
python GAS_validation_ScalingRelation.py z0p25 LC1800 # DONE
python GAS_validation_ScalingRelation.py z0p30 LC1800 # DONE
python GAS_validation_ScalingRelation.py z0p36 LC1800 # DONE
python GAS_validation_ScalingRelation.py z0p43 LC1800 # DONE
python GAS_validation_ScalingRelation.py z0p49 LC1800 # DONE
python GAS_validation_ScalingRelation.py z0p56 LC1800 # DONE
python GAS_validation_ScalingRelation.py z0p63 LC1800 # DONE
python GAS_validation_ScalingRelation.py z0p70 LC1800 # DONE
python GAS_validation_ScalingRelation.py z0p78 LC1800 # DONE
python GAS_validation_ScalingRelation.py z0p86 LC1800 # DONE
python GAS_validation_ScalingRelation.py z0p94 LC1800 # DONE

python GAS_validation_ScalingRelation.py z1p03 LC0060 # DONE
python GAS_validation_ScalingRelation.py z1p12 LC0060 # DONE
python GAS_validation_ScalingRelation.py z1p22 LC0060 # DONE
python GAS_validation_ScalingRelation.py z1p32 LC0060 # DONE
python GAS_validation_ScalingRelation.py z1p43 LC0060 # DONE
python GAS_validation_ScalingRelation.py z1p54 LC0060 # DONE
python GAS_validation_ScalingRelation.py z1p65 LC0060 # DONE
python GAS_validation_ScalingRelation.py z1p77 LC0060 # DONE
python GAS_validation_ScalingRelation.py z1p90 LC0060 # DONE
python GAS_validation_ScalingRelation.py z2p03 LC0060 # DONE
python GAS_validation_ScalingRelation.py z2p17 LC0060 # DONE
python GAS_validation_ScalingRelation.py z2p31 LC0060 # DONE
python GAS_validation_ScalingRelation.py z2p46 LC0060 # DONE
python GAS_validation_ScalingRelation.py z2p62 LC0060 # DONE
python GAS_validation_ScalingRelation.py z2p78 LC0060 # DONE
python GAS_validation_ScalingRelation.py z2p95 LC0060 # DONE
python GAS_validation_ScalingRelation.py z3p13 LC0060 # DONE
python GAS_validation_ScalingRelation.py z3p32 LC0060 # DONE
python GAS_validation_ScalingRelation.py z3p61 LC0060 # DONE
python GAS_validation_ScalingRelation.py z3p93 LC0060 # DONE

python GAS_validation_ScalingRelation.py z4p27 LC0002 # DONE
python GAS_validation_ScalingRelation.py z4p63 LC0002 # DONE
python GAS_validation_ScalingRelation.py z5p15 LC0002 # DONE

#
# more to run, but probably not useful !

python GAS_validation_ScalingRelation.py z0p00 LC0002 # DONE
python GAS_validation_ScalingRelation.py z0p02 LC0002 # DONE
python GAS_validation_ScalingRelation.py z0p05 LC0002 # DONE
python GAS_validation_ScalingRelation.py z0p09 LC0002 # DONE
python GAS_validation_ScalingRelation.py z0p14 LC0002 # DONE
python GAS_validation_ScalingRelation.py z0p19 LC0002 # DONE
python GAS_validation_ScalingRelation.py z0p25 LC0002 # DONE
python GAS_validation_ScalingRelation.py z0p30 LC0002 # DONE
python GAS_validation_ScalingRelation.py z0p36 LC0002 # DONE
python GAS_validation_ScalingRelation.py z0p43 LC0002 # DONE
python GAS_validation_ScalingRelation.py z0p49 LC0002 # DONE
python GAS_validation_ScalingRelation.py z0p56 LC0002 # DONE
python GAS_validation_ScalingRelation.py z0p63 LC0002 # DONE
python GAS_validation_ScalingRelation.py z0p70 LC0002 # DONE
python GAS_validation_ScalingRelation.py z0p78 LC0002 # DONE
python GAS_validation_ScalingRelation.py z0p86 LC0002 # DONE
python GAS_validation_ScalingRelation.py z0p94 LC0002 # DONE
python GAS_validation_ScalingRelation.py z1p03 LC0002 # DONE
python GAS_validation_ScalingRelation.py z1p12 LC0002 # DONE
python GAS_validation_ScalingRelation.py z1p22 LC0002 # DONE
python GAS_validation_ScalingRelation.py z1p32 LC0002 # DONE
python GAS_validation_ScalingRelation.py z1p43 LC0002 # DONE
python GAS_validation_ScalingRelation.py z1p54 LC0002 # DONE
python GAS_validation_ScalingRelation.py z1p65 LC0002 # DONE
python GAS_validation_ScalingRelation.py z1p77 LC0002 # DONE
python GAS_validation_ScalingRelation.py z1p90 LC0002 # DONE
python GAS_validation_ScalingRelation.py z2p03 LC0002 # DONE
python GAS_validation_ScalingRelation.py z2p17 LC0002 # DONE
python GAS_validation_ScalingRelation.py z2p31 LC0002 # DONE
python GAS_validation_ScalingRelation.py z2p46 LC0002 # DONE
python GAS_validation_ScalingRelation.py z2p62 LC0002 # DONE
python GAS_validation_ScalingRelation.py z2p78 LC0002 # DONE
python GAS_validation_ScalingRelation.py z2p95 LC0002 # DONE
python GAS_validation_ScalingRelation.py z3p13 LC0002 # DONE
python GAS_validation_ScalingRelation.py z3p32 LC0002 # DONE
python GAS_validation_ScalingRelation.py z3p61 LC0002 # DONE
python GAS_validation_ScalingRelation.py z3p93 LC0002 # DONE
python GAS_validation_ScalingRelation.py z4p27 LC0002 # DONE
python GAS_validation_ScalingRelation.py z4p63 LC0002 # DONE
python GAS_validation_ScalingRelation.py z5p15 LC0002 # DONE
python GAS_validation_ScalingRelation.py z5p73 LC0002 # DONE

python GAS_validation_ScalingRelation.py z0p00 LC0060 # DONE
python GAS_validation_ScalingRelation.py z0p02 LC0060 # DONE
python GAS_validation_ScalingRelation.py z0p05 LC0060 # DONE
python GAS_validation_ScalingRelation.py z0p09 LC0060 # DONE
python GAS_validation_ScalingRelation.py z0p14 LC0060 # DONE
python GAS_validation_ScalingRelation.py z0p19 LC0060 # DONE
python GAS_validation_ScalingRelation.py z0p25 LC0060 # DONE
python GAS_validation_ScalingRelation.py z0p30 LC0060 # DONE
python GAS_validation_ScalingRelation.py z0p36 LC0060 # DONE
python GAS_validation_ScalingRelation.py z0p43 LC0060 # DONE
python GAS_validation_ScalingRelation.py z0p49 LC0060 # DONE
python GAS_validation_ScalingRelation.py z0p56 LC0060 # DONE
python GAS_validation_ScalingRelation.py z0p63 LC0060 # DONE
python GAS_validation_ScalingRelation.py z0p70 LC0060 # DONE
python GAS_validation_ScalingRelation.py z0p78 LC0060 # DONE
python GAS_validation_ScalingRelation.py z0p86 LC0060 # DONE
python GAS_validation_ScalingRelation.py z0p94 LC0060 # DONE
python GAS_validation_ScalingRelation.py z1p03 LC0060 # DONE
python GAS_validation_ScalingRelation.py z1p12 LC0060 # DONE
python GAS_validation_ScalingRelation.py z1p22 LC0060 # DONE
python GAS_validation_ScalingRelation.py z1p32 LC0060 # DONE
python GAS_validation_ScalingRelation.py z1p43 LC0060 # DONE
python GAS_validation_ScalingRelation.py z1p54 LC0060 # DONE
python GAS_validation_ScalingRelation.py z1p65 LC0060 # DONE
python GAS_validation_ScalingRelation.py z1p77 LC0060 # DONE
python GAS_validation_ScalingRelation.py z1p90 LC0060 # DONE
python GAS_validation_ScalingRelation.py z2p03 LC0060 # DONE
python GAS_validation_ScalingRelation.py z2p17 LC0060 # DONE
python GAS_validation_ScalingRelation.py z2p31 LC0060 # DONE
python GAS_validation_ScalingRelation.py z2p46 LC0060 # DONE
python GAS_validation_ScalingRelation.py z2p62 LC0060 # DONE
python GAS_validation_ScalingRelation.py z2p78 LC0060 # DONE
python GAS_validation_ScalingRelation.py z2p95 LC0060 # DONE
python GAS_validation_ScalingRelation.py z3p13 LC0060 # DONE
python GAS_validation_ScalingRelation.py z3p32 LC0060 # DONE
python GAS_validation_ScalingRelation.py z3p61 LC0060 # DONE
python GAS_validation_ScalingRelation.py z3p93 LC0060 # DONE


python GAS_validation_ScalingRelation.py z0p00 LC1800 # DONE
python GAS_validation_ScalingRelation.py z0p02 LC1800 # DONE
python GAS_validation_ScalingRelation.py z0p05 LC1800 # DONE
python GAS_validation_ScalingRelation.py z0p09 LC1800 # DONE
python GAS_validation_ScalingRelation.py z0p14 LC1800 # DONE
python GAS_validation_ScalingRelation.py z0p19 LC1800 # DONE
python GAS_validation_ScalingRelation.py z0p25 LC1800 # DONE
python GAS_validation_ScalingRelation.py z0p30 LC1800 # DONE
python GAS_validation_ScalingRelation.py z0p36 LC1800 # DONE
python GAS_validation_ScalingRelation.py z0p43 LC1800 # DONE
python GAS_validation_ScalingRelation.py z0p49 LC1800 # DONE
python GAS_validation_ScalingRelation.py z0p56 LC1800 # DONE
python GAS_validation_ScalingRelation.py z0p63 LC1800 # DONE
python GAS_validation_ScalingRelation.py z0p70 LC1800 # DONE
python GAS_validation_ScalingRelation.py z0p78 LC1800 # DONE
python GAS_validation_ScalingRelation.py z0p86 LC1800 # DONE
python GAS_validation_ScalingRelation.py z0p94 LC1800 # DONE
python GAS_validation_ScalingRelation.py z1p03 LC1800 # DONE
python GAS_validation_ScalingRelation.py z1p12 LC1800 # DONE
python GAS_validation_ScalingRelation.py z1p22 LC1800 # DONE
python GAS_validation_ScalingRelation.py z1p32 LC1800 # DONE
python GAS_validation_ScalingRelation.py z1p43 LC1800 # DONE
python GAS_validation_ScalingRelation.py z1p54 LC1800 # DONE
python GAS_validation_ScalingRelation.py z1p65 LC1800 # DONE
python GAS_validation_ScalingRelation.py z1p77 LC1800 # DONE
python GAS_validation_ScalingRelation.py z1p90 LC1800 # DONE
python GAS_validation_ScalingRelation.py z2p03 LC1800 # DONE
