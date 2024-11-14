#!/bin/bash
conda activate stmod
cd $GIT_STMOD/src/Mpipelines

# attempt with GAMA
# creates gal_GP_model_GAMAtraining.pkl
python GAL_setup_zMsSFR_GAMA.py
# uses gal_GP_model_GAMAtraining.pkl on COSMOS to see the outcome
python GAL_setup_zMsSFR_GAMA_test_COSMOS.py
# limited in prediction power since the E(B-V) is not available in the data.

# attempt with COSMOS
# create gal_GP_model_COSMOStraining.pkl
python GAL_setup_zMsSFR_COSMOS.py
# good performances

# attempt with LS10, but does not have SFR !
# create gal_GP_model_COSMOStraining.pkl
# python GAL_setup_zMsSFR_LS10.py

# in ds43 screens, verify it rusn without memory issues, otherwise, loop over small parts of the files.

python GAL_pipeline_zMsSfr.py z0p00 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z0p02 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z0p05 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z0p09 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z0p14 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z0p19 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z0p25 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z0p30 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z0p36 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z0p43 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z0p49 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z0p56 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z0p63 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z0p70 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z0p78 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z0p86 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z0p94 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z1p03 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z1p12 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z1p22 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z1p32 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z1p43 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z1p54 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z1p65 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z1p77 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z1p90 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z2p03 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z2p17 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z2p31 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z2p46 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z2p62 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z2p78 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z2p95 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z3p13 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z3p32 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z3p61 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z3p93 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z4p27 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z4p63 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z5p15 LC0002 # DONE
python GAL_pipeline_zMsSfr.py z5p73 LC0002 # DONE

conda activate stmod
cd $GIT_STMOD/src/Mpipelines


python GAL_pipeline_zMsSfr.py z0p00 FullSky # DONE
python GAL_pipeline_zMsSfr.py z0p02 FullSky # DONE
python GAL_pipeline_zMsSfr.py z0p05 FullSky # ds43_3
python GAL_pipeline_zMsSfr.py z0p09 FullSky # ds43_3
python GAL_pipeline_zMsSfr.py z0p14 FullSky # ds43
python GAL_pipeline_zMsSfr.py z0p19 FullSky # ds43
python GAL_pipeline_zMsSfr.py z0p25 FullSky # ds52 #
python GAL_pipeline_zMsSfr.py z0p30 FullSky # ds52 #
python GAL_pipeline_zMsSfr.py z0p36 FullSky # ds52 #
python GAL_pipeline_zMsSfr.py z0p43 FullSky # ds52 #
python GAL_pipeline_zMsSfr.py z0p49 FullSky # ds52
python GAL_pipeline_zMsSfr.py z0p56 FullSky # ds52
python GAL_pipeline_zMsSfr.py z0p63 FullSky # ds52
python GAL_pipeline_zMsSfr.py z0p70 FullSky # ds52
python GAL_pipeline_zMsSfr.py z0p78 FullSky # ds52 2
python GAL_pipeline_zMsSfr.py z0p86 FullSky # ds52 2
python GAL_pipeline_zMsSfr.py z0p94 FullSky # ds52 2
python GAL_pipeline_zMsSfr.py z1p03 FullSky # ds52 2
python GAL_pipeline_zMsSfr.py z1p12 FullSky # ds52 2
python GAL_pipeline_zMsSfr.py z1p22 FullSky #ds43
python GAL_pipeline_zMsSfr.py z1p32 FullSky #ds43
python GAL_pipeline_zMsSfr.py z1p43 FullSky #ds43
python GAL_pipeline_zMsSfr.py z1p54 FullSky #ds43
python GAL_pipeline_zMsSfr.py z1p65 FullSky #ds43
python GAL_pipeline_zMsSfr.py z1p77 FullSky #ds43
python GAL_pipeline_zMsSfr.py z1p90 FullSky
python GAL_pipeline_zMsSfr.py z2p03 FullSky
python GAL_pipeline_zMsSfr.py z2p17 FullSky
python GAL_pipeline_zMsSfr.py z2p31 FullSky
python GAL_pipeline_zMsSfr.py z2p46 FullSky
python GAL_pipeline_zMsSfr.py z2p62 FullSky
python GAL_pipeline_zMsSfr.py z2p78 FullSky
python GAL_pipeline_zMsSfr.py z2p95 FullSky
python GAL_pipeline_zMsSfr.py z3p13 FullSky
python GAL_pipeline_zMsSfr.py z3p32 FullSky
python GAL_pipeline_zMsSfr.py z3p61 FullSky
python GAL_pipeline_zMsSfr.py z3p93 FullSky
python GAL_pipeline_zMsSfr.py z4p27 FullSky
python GAL_pipeline_zMsSfr.py z4p63 FullSky
python GAL_pipeline_zMsSfr.py z5p15 FullSky
python GAL_pipeline_zMsSfr.py z5p73 FullSky


conda activate stmod
cd $GIT_STMOD/src/Mpipelines


python GAL_pipeline_zMsSfr.py z0p00 LC0060
python GAL_pipeline_zMsSfr.py z0p02 LC0060
python GAL_pipeline_zMsSfr.py z0p05 LC0060
python GAL_pipeline_zMsSfr.py z0p09 LC0060
python GAL_pipeline_zMsSfr.py z0p14 LC0060
python GAL_pipeline_zMsSfr.py z0p19 LC0060
python GAL_pipeline_zMsSfr.py z0p25 LC0060
python GAL_pipeline_zMsSfr.py z0p30 LC0060
python GAL_pipeline_zMsSfr.py z0p36 LC0060
python GAL_pipeline_zMsSfr.py z0p43 LC0060
python GAL_pipeline_zMsSfr.py z0p49 LC0060
python GAL_pipeline_zMsSfr.py z0p56 LC0060
python GAL_pipeline_zMsSfr.py z0p63 LC0060
python GAL_pipeline_zMsSfr.py z0p70 LC0060
python GAL_pipeline_zMsSfr.py z0p78 LC0060
python GAL_pipeline_zMsSfr.py z0p86 LC0060
python GAL_pipeline_zMsSfr.py z0p94 LC0060
python GAL_pipeline_zMsSfr.py z1p03 LC0060
python GAL_pipeline_zMsSfr.py z1p12 LC0060
python GAL_pipeline_zMsSfr.py z1p22 LC0060
python GAL_pipeline_zMsSfr.py z1p32 LC0060
python GAL_pipeline_zMsSfr.py z1p43 LC0060
python GAL_pipeline_zMsSfr.py z1p54 LC0060
python GAL_pipeline_zMsSfr.py z1p65 LC0060
python GAL_pipeline_zMsSfr.py z1p77 LC0060
python GAL_pipeline_zMsSfr.py z1p90 LC0060
python GAL_pipeline_zMsSfr.py z2p03 LC0060
python GAL_pipeline_zMsSfr.py z2p17 LC0060
python GAL_pipeline_zMsSfr.py z2p31 LC0060
python GAL_pipeline_zMsSfr.py z2p46 LC0060
python GAL_pipeline_zMsSfr.py z2p62 LC0060
python GAL_pipeline_zMsSfr.py z2p78 LC0060
python GAL_pipeline_zMsSfr.py z2p95 LC0060
python GAL_pipeline_zMsSfr.py z3p13 LC0060
python GAL_pipeline_zMsSfr.py z3p32 LC0060
python GAL_pipeline_zMsSfr.py z3p61 LC0060
python GAL_pipeline_zMsSfr.py z3p93 LC0060

conda activate stmod
cd $GIT_STMOD/src/Mpipelines


python GAL_pipeline_zMsSfr.py z0p00 LC1800
python GAL_pipeline_zMsSfr.py z0p02 LC1800
python GAL_pipeline_zMsSfr.py z0p05 LC1800
python GAL_pipeline_zMsSfr.py z0p09 LC1800
python GAL_pipeline_zMsSfr.py z0p14 LC1800
python GAL_pipeline_zMsSfr.py z0p19 LC1800
python GAL_pipeline_zMsSfr.py z0p25 LC1800
python GAL_pipeline_zMsSfr.py z0p30 LC1800
python GAL_pipeline_zMsSfr.py z0p36 LC1800
python GAL_pipeline_zMsSfr.py z0p43 LC1800
python GAL_pipeline_zMsSfr.py z0p49 LC1800
python GAL_pipeline_zMsSfr.py z0p56 LC1800
python GAL_pipeline_zMsSfr.py z0p63 LC1800
python GAL_pipeline_zMsSfr.py z0p70 LC1800
python GAL_pipeline_zMsSfr.py z0p78 LC1800
python GAL_pipeline_zMsSfr.py z0p86 LC1800
python GAL_pipeline_zMsSfr.py z0p94 LC1800
python GAL_pipeline_zMsSfr.py z1p03 LC1800
python GAL_pipeline_zMsSfr.py z1p12 LC1800
python GAL_pipeline_zMsSfr.py z1p22 LC1800
python GAL_pipeline_zMsSfr.py z1p32 LC1800
python GAL_pipeline_zMsSfr.py z1p43 LC1800
python GAL_pipeline_zMsSfr.py z1p54 LC1800
python GAL_pipeline_zMsSfr.py z1p65 LC1800
python GAL_pipeline_zMsSfr.py z1p77 LC1800
python GAL_pipeline_zMsSfr.py z1p90 LC1800
python GAL_pipeline_zMsSfr.py z2p03 LC1800
