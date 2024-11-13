#!/bin/bash
conda activate stmod
cd $GIT_STMOD/src/Mpipelines


python GAL_validation_zMsSFR.py z0p00 LC0060
python GAL_validation_zMsSFR.py z0p02 LC0060
python GAL_validation_zMsSFR.py z0p05 LC0060
python GAL_validation_zMsSFR.py z0p09 LC0060
python GAL_validation_zMsSFR.py z0p14 LC0060
python GAL_validation_zMsSFR.py z0p19 LC0060
python GAL_validation_zMsSFR.py z0p25 LC0060
python GAL_validation_zMsSFR.py z0p30 LC0060
python GAL_validation_zMsSFR.py z0p36 LC0060
python GAL_validation_zMsSFR.py z0p43 LC0060
python GAL_validation_zMsSFR.py z0p49 LC0060
python GAL_validation_zMsSFR.py z0p56 LC0060
python GAL_validation_zMsSFR.py z0p63 LC0060
python GAL_validation_zMsSFR.py z0p70 LC0060
python GAL_validation_zMsSFR.py z0p78 LC0060
python GAL_validation_zMsSFR.py z0p86 LC0060
python GAL_validation_zMsSFR.py z0p94 LC0060

python GAL_validation_zMsSFR.py z0p00 LC1800
python GAL_validation_zMsSFR.py z0p02 LC1800
python GAL_validation_zMsSFR.py z0p05 LC1800
python GAL_validation_zMsSFR.py z0p09 LC1800
python GAL_validation_zMsSFR.py z0p14 LC1800
python GAL_validation_zMsSFR.py z0p19 LC1800
python GAL_validation_zMsSFR.py z0p25 LC1800
python GAL_validation_zMsSFR.py z0p30 LC1800
python GAL_validation_zMsSFR.py z0p36 LC1800
python GAL_validation_zMsSFR.py z0p43 LC1800
python GAL_validation_zMsSFR.py z0p49 LC1800
python GAL_validation_zMsSFR.py z0p56 LC1800
python GAL_validation_zMsSFR.py z0p63 LC1800
python GAL_validation_zMsSFR.py z0p70 LC1800
python GAL_validation_zMsSFR.py z0p78 LC1800
python GAL_validation_zMsSFR.py z0p86 LC1800
python GAL_validation_zMsSFR.py z0p94 LC1800
# python GAL_validation_zMsSFR.py z1p03 LC1800
# python GAL_validation_zMsSFR.py z1p12 LC1800
# python GAL_validation_zMsSFR.py z1p22 LC1800
# python GAL_validation_zMsSFR.py z1p32 LC1800
# python GAL_validation_zMsSFR.py z1p43 LC1800
# python GAL_validation_zMsSFR.py z1p54 LC1800
# python GAL_validation_zMsSFR.py z1p65 LC1800
# python GAL_validation_zMsSFR.py z1p77 LC1800
# python GAL_validation_zMsSFR.py z1p90 LC1800
# python GAL_validation_zMsSFR.py z2p03 LC1800
# python GAL_validation_zMsSFR.py z2p17 LC1800
# python GAL_validation_zMsSFR.py z2p31 LC1800
# python GAL_validation_zMsSFR.py z2p46 LC1800
# python GAL_validation_zMsSFR.py z2p62 LC1800
# python GAL_validation_zMsSFR.py z2p78 LC1800
# python GAL_validation_zMsSFR.py z2p95 LC1800
# python GAL_validation_zMsSFR.py z3p13 LC1800
# python GAL_validation_zMsSFR.py z3p32 LC1800
# python GAL_validation_zMsSFR.py z3p61 LC1800
# python GAL_validation_zMsSFR.py z3p93 LC1800
# python GAL_validation_zMsSFR.py z4p27 LC1800
# python GAL_validation_zMsSFR.py z4p63 LC1800
# python GAL_validation_zMsSFR.py z5p15 LC1800
# python GAL_validation_zMsSFR.py z5p73 LC1800

# cd $UCHUU/FullSky/z0p05/replication_0.0_0.0_0.0/
# rsync -avz comparat@ds43:/data56s/comparat/UCHUU/FullSky/z0p05/replication_0.0_0.0_0.0/zMsSFRmatch_mags.fits .

cd $GIT_STMOD_DATA/data/validation/validation_GAL
rsync -avz comparat@ds43:/home/comparat/st_mod_data/data/validation/validation_GAL/GP_* .
