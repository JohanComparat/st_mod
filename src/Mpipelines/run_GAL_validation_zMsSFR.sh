#!/bin/bash
conda activate stmod
cd $GIT_STMOD/src/Mpipelines


python GAL_validation_zMsSFR.py z0p00 FullSky
python GAL_validation_zMsSFR.py z0p02 FullSky
python GAL_validation_zMsSFR.py z0p05 FullSky
python GAL_validation_zMsSFR.py z0p09 FullSky
python GAL_validation_zMsSFR.py z0p14 FullSky
python GAL_validation_zMsSFR.py z0p19 FullSky
python GAL_validation_zMsSFR.py z0p25 FullSky
python GAL_validation_zMsSFR.py z0p30 FullSky
python GAL_validation_zMsSFR.py z0p36 FullSky
python GAL_validation_zMsSFR.py z0p43 FullSky
python GAL_validation_zMsSFR.py z0p49 FullSky
python GAL_validation_zMsSFR.py z0p56 FullSky
python GAL_validation_zMsSFR.py z0p63 FullSky
python GAL_validation_zMsSFR.py z0p70 FullSky
python GAL_validation_zMsSFR.py z0p78 FullSky
python GAL_validation_zMsSFR.py z0p86 FullSky
python GAL_validation_zMsSFR.py z0p94 FullSky

python GAL_validation_zMsSFR.py z0p00 LC0020
python GAL_validation_zMsSFR.py z0p02 LC0020
python GAL_validation_zMsSFR.py z0p05 LC0020
python GAL_validation_zMsSFR.py z0p09 LC0020
python GAL_validation_zMsSFR.py z0p14 LC0020
python GAL_validation_zMsSFR.py z0p19 LC0020
python GAL_validation_zMsSFR.py z0p25 LC0020
python GAL_validation_zMsSFR.py z0p30 LC0020
python GAL_validation_zMsSFR.py z0p36 LC0020
python GAL_validation_zMsSFR.py z0p43 LC0020
python GAL_validation_zMsSFR.py z0p49 LC0020
python GAL_validation_zMsSFR.py z0p56 LC0020
python GAL_validation_zMsSFR.py z0p63 LC0020
python GAL_validation_zMsSFR.py z0p70 LC0020
python GAL_validation_zMsSFR.py z0p78 LC0020
python GAL_validation_zMsSFR.py z0p86 LC0020
python GAL_validation_zMsSFR.py z0p94 LC0020
python GAL_validation_zMsSFR.py z1p03 LC0020
python GAL_validation_zMsSFR.py z1p12 LC0020
python GAL_validation_zMsSFR.py z1p22 LC0020
python GAL_validation_zMsSFR.py z1p32 LC0020
python GAL_validation_zMsSFR.py z1p43 LC0020
python GAL_validation_zMsSFR.py z1p54 LC0020
python GAL_validation_zMsSFR.py z1p65 LC0020
python GAL_validation_zMsSFR.py z1p77 LC0020
python GAL_validation_zMsSFR.py z1p90 LC0020
python GAL_validation_zMsSFR.py z2p03 LC0020
python GAL_validation_zMsSFR.py z2p17 LC0020
python GAL_validation_zMsSFR.py z2p31 LC0020
python GAL_validation_zMsSFR.py z2p46 LC0020
python GAL_validation_zMsSFR.py z2p62 LC0020
python GAL_validation_zMsSFR.py z2p78 LC0020
python GAL_validation_zMsSFR.py z2p95 LC0020
python GAL_validation_zMsSFR.py z3p13 LC0020
python GAL_validation_zMsSFR.py z3p32 LC0020
python GAL_validation_zMsSFR.py z3p61 LC0020
python GAL_validation_zMsSFR.py z3p93 LC0020
python GAL_validation_zMsSFR.py z4p27 LC0020
python GAL_validation_zMsSFR.py z4p63 LC0020
python GAL_validation_zMsSFR.py z5p15 LC0020
python GAL_validation_zMsSFR.py z5p73 LC0020

# cd $UCHUU/FullSky/z0p05/replication_0.0_0.0_0.0/
# rsync -avz comparat@ds43:/data56s/comparat/UCHUU/FullSky/z0p05/replication_0.0_0.0_0.0/zMsSFRmatch_mags.fits .

cd $GIT_STMOD_DATA/data/validation/validation_GAL
rsync -avz comparat@ds43:/home/comparat/st_mod_data/data/validation/validation_GAL/GP_* .
