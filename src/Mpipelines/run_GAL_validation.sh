#!/bin/bash
conda activate stmod
cd $GIT_STMOD/src/Mpipelines
python GAL_validation.py z0p00
python GAL_validation.py z0p02
python GAL_validation.py z0p05
python GAL_validation.py z0p09
python GAL_validation.py z0p14
python GAL_validation.py z0p19
python GAL_validation.py z0p25
python GAL_validation.py z0p30
python GAL_validation.py z0p36
python GAL_validation.py z0p43
python GAL_validation.py z0p49
python GAL_validation.py z0p56
python GAL_validation.py z0p63
python GAL_validation.py z0p70
python GAL_validation.py z0p78
python GAL_validation.py z0p86
python GAL_validation.py z0p94
python GAL_validation.py z1p03
python GAL_validation.py z1p12
python GAL_validation.py z1p22
python GAL_validation.py z1p32
python GAL_validation.py z1p43
python GAL_validation.py z1p54
python GAL_validation.py z1p65
python GAL_validation.py z1p77
python GAL_validation.py z1p90
python GAL_validation.py z2p03
python GAL_validation.py z2p17
python GAL_validation.py z2p31
python GAL_validation.py z2p46
python GAL_validation.py z2p62
python GAL_validation.py z2p78
python GAL_validation.py z2p95
python GAL_validation.py z3p13
python GAL_validation.py z3p32
python GAL_validation.py z3p61
python GAL_validation.py z3p93
python GAL_validation.py z4p27
python GAL_validation.py z4p63
python GAL_validation.py z5p15
python GAL_validation.py z5p73

# python GAL_validation_WPRP.py z0p02
python GAL_validation_WPRP.py z0p05
python GAL_validation_WPRP.py z0p09
python GAL_validation_WPRP.py z0p14
python GAL_validation_WPRP.py z0p19
python GAL_validation_WPRP.py z0p25
python GAL_validation_WPRP.py z0p30
python GAL_validation_WPRP.py z0p36
python GAL_validation_WPRP.py z0p43
python GAL_validation_WPRP.py z0p49
python GAL_validation_WPRP.py z0p56

cd ~/stack_models/data/validation/validation_GAL/WPRP
git add *
git commit -m"wprp"
git push
