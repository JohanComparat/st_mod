#!/bin/bash
conda activate stmod
cd $GIT_STMOD/src/Mpipelines



# The GaussianProcess ML tool fails to extrapolate where the UniverseMachine model predicts the existence of galaxies that have no counterpart in current observations.

# attempt with GAMA
# creates gal_GP_model_GAMAtraining.pkl
# python GAL_setup_zMsSFR_GAMA.py
# uses gal_GP_model_GAMAtraining.pkl on COSMOS to see the outcome
# python GAL_setup_zMsSFR_GAMA_test_COSMOS.py
# limited in prediction power since the E(B-V) is not available in the data.

# attempt with COSMOS
# create gal_GP_model_COSMOStraining.pkl
# python GAL_setup_zMsSFR_COSMOS.py
# good performances

# attempt with LS10, but does not have SFR ! and limited to too low redshift !
# create gal_GP_model_COSMOStraining.pkl
# python GAL_setup_zMsSFR_LS10.py

# classifies galaxies in red sequence and blue cloud :
python GAL_setup_zMsSFR.py # DONE



# nohup python GAL_pipeline_zMsSfr.py z0p00 LC0002 > logs/GAL_pipeline_zMsSfr_z0p00_LC0002.log & # no galaxies
# nohup python GAL_pipeline_zMsSfr.py z0p02 LC0002 > logs/GAL_pipeline_zMsSfr_z0p02_LC0002.log & # DONE
# nohup python GAL_pipeline_zMsSfr.py z0p05 LC0002 > logs/GAL_pipeline_zMsSfr_z0p05_LC0002.log & # DONE
# nohup python GAL_pipeline_zMsSfr.py z0p09 LC0002 > logs/GAL_pipeline_zMsSfr_z0p09_LC0002.log & # DONE
# nohup python GAL_pipeline_zMsSfr.py z0p14 LC0002 > logs/GAL_pipeline_zMsSfr_z0p14_LC0002.log & # DONE
# nohup python GAL_pipeline_zMsSfr.py z0p19 LC0002 > logs/GAL_pipeline_zMsSfr_z0p19_LC0002.log & # DONE
# nohup python GAL_pipeline_zMsSfr.py z0p25 LC0002 > logs/GAL_pipeline_zMsSfr_z0p25_LC0002.log & # DONE
# nohup python GAL_pipeline_zMsSfr.py z0p30 LC0002 > logs/GAL_pipeline_zMsSfr_z0p30_LC0002.log & # DONE
# nohup python GAL_pipeline_zMsSfr.py z0p36 LC0002 > logs/GAL_pipeline_zMsSfr_z0p36_LC0002.log & # DONE
# nohup python GAL_pipeline_zMsSfr.py z0p43 LC0002 > logs/GAL_pipeline_zMsSfr_z0p43_LC0002.log & # DONE
# nohup python GAL_pipeline_zMsSfr.py z0p49 LC0002 > logs/GAL_pipeline_zMsSfr_z0p49_LC0002.log & # DONE
# nohup python GAL_pipeline_zMsSfr.py z0p56 LC0002 > logs/GAL_pipeline_zMsSfr_z0p56_LC0002.log & # DONE
# nohup python GAL_pipeline_zMsSfr.py z0p63 LC0002 > logs/GAL_pipeline_zMsSfr_z0p63_LC0002.log & # DONE
# nohup python GAL_pipeline_zMsSfr.py z0p70 LC0002 > logs/GAL_pipeline_zMsSfr_z0p70_LC0002.log & # DONE
# nohup python GAL_pipeline_zMsSfr.py z0p78 LC0002 > logs/GAL_pipeline_zMsSfr_z0p78_LC0002.log & # DONE
# nohup python GAL_pipeline_zMsSfr.py z0p86 LC0002 > logs/GAL_pipeline_zMsSfr_z0p86_LC0002.log & # DONE
# nohup python GAL_pipeline_zMsSfr.py z0p94 LC0002 > logs/GAL_pipeline_zMsSfr_z0p94_LC0002.log & # DONE
#
# nohup python GAL_pipeline_zMsSfr.py z0p00 LC0060 > logs/GAL_pipeline_zMsSfr_z0p00_LC0060.log & # DONE
# nohup python GAL_pipeline_zMsSfr.py z0p02 LC0060 > logs/GAL_pipeline_zMsSfr_z0p02_LC0060.log & # DONE
# nohup python GAL_pipeline_zMsSfr.py z0p05 LC0060 > logs/GAL_pipeline_zMsSfr_z0p05_LC0060.log & # DONE
# nohup python GAL_pipeline_zMsSfr.py z0p09 LC0060 > logs/GAL_pipeline_zMsSfr_z0p09_LC0060.log & # DONE
# nohup python GAL_pipeline_zMsSfr.py z0p14 LC0060 > logs/GAL_pipeline_zMsSfr_z0p14_LC0060.log & # DONE
# nohup python GAL_pipeline_zMsSfr.py z0p19 LC0060 > logs/GAL_pipeline_zMsSfr_z0p19_LC0060.log & # DONE
# nohup python GAL_pipeline_zMsSfr.py z0p25 LC0060 > logs/GAL_pipeline_zMsSfr_z0p25_LC0060.log & # DONE
# nohup python GAL_pipeline_zMsSfr.py z0p30 LC0060 > logs/GAL_pipeline_zMsSfr_z0p30_LC0060.log & # DONE
# nohup python GAL_pipeline_zMsSfr.py z0p36 LC0060 > logs/GAL_pipeline_zMsSfr_z0p36_LC0060.log & # DONE
# nohup python GAL_pipeline_zMsSfr.py z0p43 LC0060 > logs/GAL_pipeline_zMsSfr_z0p43_LC0060.log & # DONE
# nohup python GAL_pipeline_zMsSfr.py z0p49 LC0060 > logs/GAL_pipeline_zMsSfr_z0p49_LC0060.log & # DONE
# nohup python GAL_pipeline_zMsSfr.py z0p56 LC0060 > logs/GAL_pipeline_zMsSfr_z0p56_LC0060.log & # DONE
# nohup python GAL_pipeline_zMsSfr.py z0p63 LC0060 > logs/GAL_pipeline_zMsSfr_z0p63_LC0060.log & # DONE
# nohup python GAL_pipeline_zMsSfr.py z0p70 LC0060 > logs/GAL_pipeline_zMsSfr_z0p70_LC0060.log & # DONE
# nohup python GAL_pipeline_zMsSfr.py z0p78 LC0060 > logs/GAL_pipeline_zMsSfr_z0p78_LC0060.log & # DONE
# nohup python GAL_pipeline_zMsSfr.py z0p86 LC0060 > logs/GAL_pipeline_zMsSfr_z0p86_LC0060.log & # DONE
# nohup python GAL_pipeline_zMsSfr.py z0p94 LC0060 > logs/GAL_pipeline_zMsSfr_z0p94_LC0060.log & # DONE

# nohup python GAL_pipeline_zMsSfr.py z0p00 LC1800 > logs/GAL_pipeline_zMsSfr_z0p00_LC1800.log & # DONE
# nohup python GAL_pipeline_zMsSfr.py z0p02 LC1800 > logs/GAL_pipeline_zMsSfr_z0p02_LC1800.log & # DONE
# nohup python GAL_pipeline_zMsSfr.py z0p05 LC1800 > logs/GAL_pipeline_zMsSfr_z0p05_LC1800.log & # DONE
# nohup python GAL_pipeline_zMsSfr.py z0p09 LC1800 > logs/GAL_pipeline_zMsSfr_z0p09_LC1800.log & # DONE
nohup python GAL_pipeline_zMsSfr.py z0p14 LC1800 > logs/GAL_pipeline_zMsSfr_z0p14_LC1800.log & # DONE
nohup python GAL_pipeline_zMsSfr.py z0p19 LC1800 > logs/GAL_pipeline_zMsSfr_z0p19_LC1800.log & # DONE
nohup python GAL_pipeline_zMsSfr.py z0p25 LC1800 > logs/GAL_pipeline_zMsSfr_z0p25_LC1800.log & # DONE
nohup python GAL_pipeline_zMsSfr.py z0p30 LC1800 > logs/GAL_pipeline_zMsSfr_z0p30_LC1800.log & # DONE
nohup python GAL_pipeline_zMsSfr.py z0p36 LC1800 > logs/GAL_pipeline_zMsSfr_z0p36_LC1800.log & # DONE
nohup python GAL_pipeline_zMsSfr.py z0p43 LC1800 > logs/GAL_pipeline_zMsSfr_z0p43_LC1800.log & # DONE
nohup python GAL_pipeline_zMsSfr.py z0p49 LC1800 > logs/GAL_pipeline_zMsSfr_z0p49_LC1800.log & # DONE
nohup python GAL_pipeline_zMsSfr.py z0p56 LC1800 > logs/GAL_pipeline_zMsSfr_z0p56_LC1800.log & # DONE
nohup python GAL_pipeline_zMsSfr.py z0p63 LC1800 > logs/GAL_pipeline_zMsSfr_z0p63_LC1800.log & # DONE
nohup python GAL_pipeline_zMsSfr.py z0p70 LC1800 > logs/GAL_pipeline_zMsSfr_z0p70_LC1800.log & # DONE
nohup python GAL_pipeline_zMsSfr.py z0p78 LC1800 > logs/GAL_pipeline_zMsSfr_z0p78_LC1800.log & # DONE
nohup python GAL_pipeline_zMsSfr.py z0p86 LC1800 > logs/GAL_pipeline_zMsSfr_z0p86_LC1800.log & # DONE
nohup python GAL_pipeline_zMsSfr.py z0p94 LC1800 > logs/GAL_pipeline_zMsSfr_z0p94_LC1800.log & # DONE


# nohup python GAL_pipeline_zMsSfr.py z0p00 FullSky > logs/GAL_pipeline_zMsSfr_z0p00_FullSky.log & # TODO
# nohup python GAL_pipeline_zMsSfr.py z0p02 FullSky > logs/GAL_pipeline_zMsSfr_z0p02_FullSky.log & # TODO
# nohup python GAL_pipeline_zMsSfr.py z0p05 FullSky > logs/GAL_pipeline_zMsSfr_z0p05_FullSky.log & # TODO
# nohup python GAL_pipeline_zMsSfr.py z0p09 FullSky > logs/GAL_pipeline_zMsSfr_z0p09_FullSky.log & # TODO
nohup python GAL_pipeline_zMsSfr.py z0p14 FullSky > logs/GAL_pipeline_zMsSfr_z0p14_FullSky.log & # TODO
nohup python GAL_pipeline_zMsSfr.py z0p19 FullSky > logs/GAL_pipeline_zMsSfr_z0p19_FullSky.log & # TODO
nohup python GAL_pipeline_zMsSfr.py z0p25 FullSky > logs/GAL_pipeline_zMsSfr_z0p25_FullSky.log & # TODO
nohup python GAL_pipeline_zMsSfr.py z0p30 FullSky > logs/GAL_pipeline_zMsSfr_z0p30_FullSky.log & # TODO
nohup python GAL_pipeline_zMsSfr.py z0p36 FullSky > logs/GAL_pipeline_zMsSfr_z0p36_FullSky.log & # TODO
nohup python GAL_pipeline_zMsSfr.py z0p43 FullSky > logs/GAL_pipeline_zMsSfr_z0p43_FullSky.log & # TODO
nohup python GAL_pipeline_zMsSfr.py z0p49 FullSky > logs/GAL_pipeline_zMsSfr_z0p49_FullSky.log & # TODO
nohup python GAL_pipeline_zMsSfr.py z0p56 FullSky > logs/GAL_pipeline_zMsSfr_z0p56_FullSky.log & # TODO
nohup python GAL_pipeline_zMsSfr.py z0p63 FullSky > logs/GAL_pipeline_zMsSfr_z0p63_FullSky.log & # TODO
nohup python GAL_pipeline_zMsSfr.py z0p70 FullSky > logs/GAL_pipeline_zMsSfr_z0p70_FullSky.log & # TODO
nohup python GAL_pipeline_zMsSfr.py z0p78 FullSky > logs/GAL_pipeline_zMsSfr_z0p78_FullSky.log & # TODO
nohup python GAL_pipeline_zMsSfr.py z0p86 FullSky > logs/GAL_pipeline_zMsSfr_z0p86_FullSky.log & # TODO
nohup python GAL_pipeline_zMsSfr.py z0p94 FullSky > logs/GAL_pipeline_zMsSfr_z0p94_FullSky.log & # TODO

