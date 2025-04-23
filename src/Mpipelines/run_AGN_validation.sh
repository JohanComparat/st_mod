#!/bin/bash
conda activate stmod
cd $GIT_STMOD/src/Mpipelines
# A minimum of 3GB memory (so 2 CPU per task) is required to load observed catalogues
# all below worked and finished :


nohup python AGN_validation.py z0p00 FullSky > logs/valid_testz0p00.log & # DONE
nohup python AGN_validation.py z0p02 FullSky > logs/valid_testz0p02.log & # DONE
nohup python AGN_validation.py z0p05 FullSky > logs/valid_testz0p05.log & # DONE
nohup python AGN_validation.py z0p09 FullSky > logs/valid_testz0p09.log & # DONE
nohup python AGN_validation.py z0p14 FullSky > logs/valid_testz0p14.log & # DONE
nohup python AGN_validation.py z0p19 FullSky > logs/valid_testz0p19.log & # DONE
nohup python AGN_validation.py z0p25 FullSky > logs/valid_testz0p25.log & # DONE
nohup python AGN_validation.py z0p30 FullSky > logs/valid_testz0p30.log & # DONE
nohup python AGN_validation.py z0p36 FullSky > logs/valid_testz0p36.log & # DONE
nohup python AGN_validation.py z0p43 FullSky > logs/valid_testz0p43.log & # DONE
nohup python AGN_validation.py z0p49 FullSky > logs/valid_testz0p49.log & # DONE
nohup python AGN_validation.py z0p56 FullSky > logs/valid_testz0p56.log & # DONE
nohup python AGN_validation.py z0p63 FullSky > logs/valid_testz0p63.log & # DONE
nohup python AGN_validation.py z0p70 FullSky > logs/valid_testz0p70.log & # DONE
nohup python AGN_validation.py z0p78 FullSky > logs/valid_testz0p78.log & # DONE
nohup python AGN_validation.py z0p86 FullSky > logs/valid_testz0p86.log & # DONE
nohup python AGN_validation.py z0p94 FullSky > logs/valid_testz0p94.log & # DONE
nohup python AGN_validation.py z1p03 FullSky > logs/valid_testz1p03.log & # DONE
nohup python AGN_validation.py z1p12 FullSky > logs/valid_testz1p12.log & # DONE
nohup python AGN_validation.py z1p22 FullSky > logs/valid_testz1p22.log & # DONE
nohup python AGN_validation.py z1p32 FullSky > logs/valid_testz1p32.log & # DONE
nohup python AGN_validation.py z1p43 FullSky > logs/valid_testz1p43.log & # DONE
nohup python AGN_validation.py z1p54 FullSky > logs/valid_testz1p54.log & # DONE
nohup python AGN_validation.py z1p65 FullSky > logs/valid_testz1p65.log & # DONE
nohup python AGN_validation.py z1p77 FullSky > logs/valid_testz1p77.log & # DONE
nohup python AGN_validation.py z1p90 FullSky > logs/valid_testz1p90.log & # DONE
nohup python AGN_validation.py z2p03 FullSky > logs/valid_testz2p03.log & # DONE
nohup python AGN_validation.py z2p17 FullSky > logs/valid_testz2p17.log & # DONE
nohup python AGN_validation.py z2p31 FullSky > logs/valid_testz2p31.log & # DONE
nohup python AGN_validation.py z2p46 FullSky > logs/valid_testz2p46.log & # DONE
nohup python AGN_validation.py z2p62 FullSky > logs/valid_testz2p62.log & # DONE
nohup python AGN_validation.py z2p78 FullSky > logs/valid_testz2p78.log & # DONE
nohup python AGN_validation.py z2p95 FullSky > logs/valid_testz2p95.log & # DONE
nohup python AGN_validation.py z3p13 FullSky > logs/valid_testz3p13.log & # DONE
nohup python AGN_validation.py z3p32 FullSky > logs/valid_testz3p32.log & # DONE
nohup python AGN_validation.py z3p61 FullSky > logs/valid_testz3p61.log & # DONE
nohup python AGN_validation.py z3p93 FullSky > logs/valid_testz3p93.log & # DONE



# conda activate clustering
# cd $GIT_STMOD/src/Mpipelines
# python AGN_validation_WPRP.py z0p00
# python AGN_validation_WPRP.py z0p02
# python AGN_validation_WPRP.py z0p05
# python AGN_validation_WPRP.py z0p09
# python AGN_validation_WPRP.py z0p14
# python AGN_validation_WPRP.py z0p19
# python AGN_validation_WPRP.py z0p25
# python AGN_validation_WPRP.py z0p30
# python AGN_validation_WPRP.py z0p36
# python AGN_validation_WPRP.py z0p43
# python AGN_validation_WPRP.py z0p49
# python AGN_validation_WPRP.py z0p56
# python AGN_validation_WPRP.py z0p63
# python AGN_validation_WPRP.py z0p70
# python AGN_validation_WPRP.py z0p78
# python AGN_validation_WPRP.py z0p86
# python AGN_validation_WPRP.py z0p94
# python AGN_validation_WPRP.py z1p03
# python AGN_validation_WPRP.py z1p12
# python AGN_validation_WPRP.py z1p22
# python AGN_validation_WPRP.py z1p32
# python AGN_validation_WPRP.py z1p43
# python AGN_validation_WPRP.py z1p54
#
#
# conda activate clustering
# cd $GIT_STMOD/src/Mpipelines
# python AGN_validation_WTH.py z0p05
# python AGN_validation_WTH.py z0p09
# python AGN_validation_WTH.py z0p14
# python AGN_validation_WTH.py z0p19
# python AGN_validation_WTH.py z0p25
# python AGN_validation_WTH.py z0p30
# python AGN_validation_WTH.py z0p36
# python AGN_validation_WTH.py z0p43
# python AGN_validation_WTH.py z0p49
# python AGN_validation_WTH.py z0p56
# python AGN_validation_WTH.py z0p63
# python AGN_validation_WTH.py z0p70
# python AGN_validation_WTH.py z0p78
# python AGN_validation_WTH.py z0p86
# python AGN_validation_WTH.py z0p94
# python AGN_validation_WTH.py z1p03
# python AGN_validation_WTH.py z1p12
# python AGN_validation_WTH.py z1p22
# python AGN_validation_WTH.py z1p32
# python AGN_validation_WTH.py z1p43
# python AGN_validation_WTH.py z1p54


python AGN_validation_logNlogS.py z0p00 FullSky
python AGN_validation_logNlogS.py z0p02 FullSky
python AGN_validation_logNlogS.py z0p05 FullSky
python AGN_validation_logNlogS.py z0p09 FullSky
python AGN_validation_logNlogS.py z0p14 FullSky
python AGN_validation_logNlogS.py z0p19 FullSky
python AGN_validation_logNlogS.py z0p25 FullSky
python AGN_validation_logNlogS.py z0p30 FullSky
python AGN_validation_logNlogS.py z0p36 FullSky
python AGN_validation_logNlogS.py z0p43 FullSky
python AGN_validation_logNlogS.py z0p49 FullSky
python AGN_validation_logNlogS.py z0p56 FullSky
python AGN_validation_logNlogS.py z0p63 FullSky
python AGN_validation_logNlogS.py z0p70 FullSky
python AGN_validation_logNlogS.py z0p78 FullSky
python AGN_validation_logNlogS.py z0p86 FullSky
python AGN_validation_logNlogS.py z0p94 FullSky
python AGN_validation_logNlogS.py z1p03 FullSky
python AGN_validation_logNlogS.py z1p12 FullSky
python AGN_validation_logNlogS.py z1p22 FullSky
python AGN_validation_logNlogS.py z1p32 FullSky
python AGN_validation_logNlogS.py z1p43 FullSky
python AGN_validation_logNlogS.py z1p54 FullSky
python AGN_validation_logNlogS.py z1p65 FullSky
python AGN_validation_logNlogS.py z1p77 FullSky
python AGN_validation_logNlogS.py z1p90 FullSky
python AGN_validation_logNlogS.py z2p03 FullSky
python AGN_validation_logNlogS.py z2p17 FullSky
python AGN_validation_logNlogS.py z2p31 FullSky
python AGN_validation_logNlogS.py z2p46 FullSky
python AGN_validation_logNlogS.py z2p62 FullSky
python AGN_validation_logNlogS.py z2p78 FullSky
python AGN_validation_logNlogS.py z2p95 FullSky
python AGN_validation_logNlogS.py z3p13 FullSky
python AGN_validation_logNlogS.py z3p32 FullSky
python AGN_validation_logNlogS.py z3p61 FullSky
python AGN_validation_logNlogS.py z3p93 FullSky



python AGN_validation_logNlogS.py z0p00 LC0002
python AGN_validation_logNlogS.py z0p02 LC0002
python AGN_validation_logNlogS.py z0p05 LC0002
python AGN_validation_logNlogS.py z0p09 LC0002
python AGN_validation_logNlogS.py z0p14 LC0002
python AGN_validation_logNlogS.py z0p19 LC0002
python AGN_validation_logNlogS.py z0p25 LC0002
python AGN_validation_logNlogS.py z0p30 LC0002
python AGN_validation_logNlogS.py z0p36 LC0002
python AGN_validation_logNlogS.py z0p43 LC0002
python AGN_validation_logNlogS.py z0p49 LC0002
python AGN_validation_logNlogS.py z0p56 LC0002
python AGN_validation_logNlogS.py z0p63 LC0002
python AGN_validation_logNlogS.py z0p70 LC0002
python AGN_validation_logNlogS.py z0p78 LC0002
python AGN_validation_logNlogS.py z0p86 LC0002
python AGN_validation_logNlogS.py z0p94 LC0002
python AGN_validation_logNlogS.py z1p03 LC0002
python AGN_validation_logNlogS.py z1p12 LC0002
python AGN_validation_logNlogS.py z1p22 LC0002
python AGN_validation_logNlogS.py z1p32 LC0002
python AGN_validation_logNlogS.py z1p43 LC0002
python AGN_validation_logNlogS.py z1p54 LC0002
python AGN_validation_logNlogS.py z1p65 LC0002
python AGN_validation_logNlogS.py z1p77 LC0002
python AGN_validation_logNlogS.py z1p90 LC0002
python AGN_validation_logNlogS.py z2p03 LC0002
python AGN_validation_logNlogS.py z2p17 LC0002
python AGN_validation_logNlogS.py z2p31 LC0002
python AGN_validation_logNlogS.py z2p46 LC0002
python AGN_validation_logNlogS.py z2p62 LC0002
python AGN_validation_logNlogS.py z2p78 LC0002
python AGN_validation_logNlogS.py z2p95 LC0002
python AGN_validation_logNlogS.py z3p13 LC0002
python AGN_validation_logNlogS.py z3p32 LC0002
python AGN_validation_logNlogS.py z3p61 LC0002
python AGN_validation_logNlogS.py z3p93 LC0002
python AGN_validation_logNlogS.py z4p27 LC0002
python AGN_validation_logNlogS.py z4p63 LC0002
python AGN_validation_logNlogS.py z5p15 LC0002
python AGN_validation_logNlogS.py z5p73 LC0002


python AGN_validation_logNlogS.py z0p00 LC0060
python AGN_validation_logNlogS.py z0p02 LC0060
python AGN_validation_logNlogS.py z0p05 LC0060
python AGN_validation_logNlogS.py z0p09 LC0060
python AGN_validation_logNlogS.py z0p14 LC0060
python AGN_validation_logNlogS.py z0p19 LC0060
python AGN_validation_logNlogS.py z0p25 LC0060
python AGN_validation_logNlogS.py z0p30 LC0060
python AGN_validation_logNlogS.py z0p36 LC0060
python AGN_validation_logNlogS.py z0p43 LC0060
python AGN_validation_logNlogS.py z0p49 LC0060
python AGN_validation_logNlogS.py z0p56 LC0060
python AGN_validation_logNlogS.py z0p63 LC0060
python AGN_validation_logNlogS.py z0p70 LC0060
python AGN_validation_logNlogS.py z0p78 LC0060
python AGN_validation_logNlogS.py z0p86 LC0060
python AGN_validation_logNlogS.py z0p94 LC0060
python AGN_validation_logNlogS.py z1p03 LC0060
python AGN_validation_logNlogS.py z1p12 LC0060
python AGN_validation_logNlogS.py z1p22 LC0060
python AGN_validation_logNlogS.py z1p32 LC0060
python AGN_validation_logNlogS.py z1p43 LC0060
python AGN_validation_logNlogS.py z1p54 LC0060
python AGN_validation_logNlogS.py z1p65 LC0060
python AGN_validation_logNlogS.py z1p77 LC0060
python AGN_validation_logNlogS.py z1p90 LC0060
python AGN_validation_logNlogS.py z2p03 LC0060
python AGN_validation_logNlogS.py z2p17 LC0060
python AGN_validation_logNlogS.py z2p31 LC0060
python AGN_validation_logNlogS.py z2p46 LC0060
python AGN_validation_logNlogS.py z2p62 LC0060
python AGN_validation_logNlogS.py z2p78 LC0060
python AGN_validation_logNlogS.py z2p95 LC0060
python AGN_validation_logNlogS.py z3p13 LC0060
python AGN_validation_logNlogS.py z3p32 LC0060
python AGN_validation_logNlogS.py z3p61 LC0060
python AGN_validation_logNlogS.py z3p93 LC0060



python AGN_validation_logNlogS.py z0p00 LC1800
python AGN_validation_logNlogS.py z0p02 LC1800
python AGN_validation_logNlogS.py z0p05 LC1800
python AGN_validation_logNlogS.py z0p09 LC1800
python AGN_validation_logNlogS.py z0p14 LC1800
python AGN_validation_logNlogS.py z0p19 LC1800
python AGN_validation_logNlogS.py z0p25 LC1800
python AGN_validation_logNlogS.py z0p30 LC1800
python AGN_validation_logNlogS.py z0p36 LC1800
python AGN_validation_logNlogS.py z0p43 LC1800
python AGN_validation_logNlogS.py z0p49 LC1800
python AGN_validation_logNlogS.py z0p56 LC1800
python AGN_validation_logNlogS.py z0p63 LC1800
python AGN_validation_logNlogS.py z0p70 LC1800
python AGN_validation_logNlogS.py z0p78 LC1800
python AGN_validation_logNlogS.py z0p86 LC1800
python AGN_validation_logNlogS.py z0p94 LC1800
python AGN_validation_logNlogS.py z1p03 LC1800
python AGN_validation_logNlogS.py z1p12 LC1800
python AGN_validation_logNlogS.py z1p22 LC1800
python AGN_validation_logNlogS.py z1p32 LC1800
python AGN_validation_logNlogS.py z1p43 LC1800
python AGN_validation_logNlogS.py z1p54 LC1800
python AGN_validation_logNlogS.py z1p65 LC1800
python AGN_validation_logNlogS.py z1p77 LC1800
python AGN_validation_logNlogS.py z1p90 LC1800
python AGN_validation_logNlogS.py z2p03 LC1800



# python AGN_pipeline.py z0p19 LC0002  0.8 0
# python AGN_pipeline.py z0p19 LC0060  0.8 0
# python AGN_pipeline.py z0p19 LC1800  0.8 0
# python AGN_pipeline.py z0p19 FullSky 0.8 0

# python AGN_validation.py z0p19 LC0002
# python AGN_validation.py z0p19 LC0060
# python AGN_validation.py z0p19 LC1800
# python AGN_validation.py z0p19 FullSky
#
# python AGN_validation_logNlogS.py z0p19 LC0002
# python AGN_validation_logNlogS.py z0p19 LC0060
# python AGN_validation_logNlogS.py z0p19 LC1800
# python AGN_validation_logNlogS.py z0p19 FullSky

python AGN_validation_MERGE_logNlogS.py
nohup python AGN_validation_PLOT_logNlogS.py FullSky > plot1.log &
nohup python AGN_validation_PLOT_logNlogS.py LC0002  > plot2.log &
nohup python AGN_validation_PLOT_logNlogS.py LC0060  > plot3.log &
nohup python AGN_validation_PLOT_logNlogS.py LC1800  > plot4.log &











python AGN_validation.py z0p00 LC0002
python AGN_validation.py z0p02 LC0002
python AGN_validation.py z0p05 LC0002
python AGN_validation.py z0p09 LC0002
python AGN_validation.py z0p14 LC0002
python AGN_validation.py z0p19 LC0002
python AGN_validation.py z0p25 LC0002
python AGN_validation.py z0p30 LC0002
python AGN_validation.py z0p36 LC0002
python AGN_validation.py z0p43 LC0002
python AGN_validation.py z0p49 LC0002
python AGN_validation.py z0p56 LC0002
python AGN_validation.py z0p63 LC0002
python AGN_validation.py z0p70 LC0002
python AGN_validation.py z0p78 LC0002
python AGN_validation.py z0p86 LC0002
python AGN_validation.py z0p94 LC0002
python AGN_validation.py z1p03 LC0002
python AGN_validation.py z1p12 LC0002
python AGN_validation.py z1p22 LC0002
python AGN_validation.py z1p32 LC0002
python AGN_validation.py z1p43 LC0002
python AGN_validation.py z1p54 LC0002
python AGN_validation.py z1p65 LC0002
python AGN_validation.py z1p77 LC0002
python AGN_validation.py z1p90 LC0002
python AGN_validation.py z2p03 LC0002
python AGN_validation.py z2p17 LC0002
python AGN_validation.py z2p31 LC0002
python AGN_validation.py z2p46 LC0002
python AGN_validation.py z2p62 LC0002
python AGN_validation.py z2p78 LC0002
python AGN_validation.py z2p95 LC0002
python AGN_validation.py z3p13 LC0002
python AGN_validation.py z3p32 LC0002
python AGN_validation.py z3p61 LC0002
python AGN_validation.py z3p93 LC0002
python AGN_validation.py z4p27 LC0002
python AGN_validation.py z4p63 LC0002
python AGN_validation.py z5p15 LC0002
python AGN_validation.py z5p73 LC0002

python AGN_validation.py z0p00 LC0060
python AGN_validation.py z0p02 LC0060
python AGN_validation.py z0p05 LC0060
python AGN_validation.py z0p09 LC0060
python AGN_validation.py z0p14 LC0060
python AGN_validation.py z0p19 LC0060
python AGN_validation.py z0p25 LC0060
python AGN_validation.py z0p30 LC0060
python AGN_validation.py z0p36 LC0060
python AGN_validation.py z0p43 LC0060
python AGN_validation.py z0p49 LC0060
python AGN_validation.py z0p56 LC0060
python AGN_validation.py z0p63 LC0060
python AGN_validation.py z0p70 LC0060
python AGN_validation.py z0p78 LC0060
python AGN_validation.py z0p86 LC0060
python AGN_validation.py z0p94 LC0060
python AGN_validation.py z1p03 LC0060
python AGN_validation.py z1p12 LC0060
python AGN_validation.py z1p22 LC0060
python AGN_validation.py z1p32 LC0060
python AGN_validation.py z1p43 LC0060
python AGN_validation.py z1p54 LC0060
python AGN_validation.py z1p65 LC0060
python AGN_validation.py z1p77 LC0060
python AGN_validation.py z1p90 LC0060
python AGN_validation.py z2p03 LC0060
python AGN_validation.py z2p17 LC0060
python AGN_validation.py z2p31 LC0060
python AGN_validation.py z2p46 LC0060
python AGN_validation.py z2p62 LC0060
python AGN_validation.py z2p78 LC0060
python AGN_validation.py z2p95 LC0060
python AGN_validation.py z3p13 LC0060
python AGN_validation.py z3p32 LC0060
python AGN_validation.py z3p61 LC0060
python AGN_validation.py z3p93 LC0060

python AGN_validation.py z0p00 LC1800
python AGN_validation.py z0p02 LC1800
python AGN_validation.py z0p05 LC1800
python AGN_validation.py z0p09 LC1800
python AGN_validation.py z0p14 LC1800
python AGN_validation.py z0p19 LC1800
python AGN_validation.py z0p25 LC1800
python AGN_validation.py z0p30 LC1800
python AGN_validation.py z0p36 LC1800
python AGN_validation.py z0p43 LC1800
python AGN_validation.py z0p49 LC1800
python AGN_validation.py z0p56 LC1800
python AGN_validation.py z0p63 LC1800
python AGN_validation.py z0p70 LC1800
python AGN_validation.py z0p78 LC1800
python AGN_validation.py z0p86 LC1800
python AGN_validation.py z0p94 LC1800
python AGN_validation.py z1p03 LC1800
python AGN_validation.py z1p12 LC1800
python AGN_validation.py z1p22 LC1800
python AGN_validation.py z1p32 LC1800
python AGN_validation.py z1p43 LC1800
python AGN_validation.py z1p54 LC1800
python AGN_validation.py z1p65 LC1800
python AGN_validation.py z1p77 LC1800
python AGN_validation.py z1p90 LC1800
python AGN_validation.py z2p03 LC1800
