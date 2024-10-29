#!/bin/bash
conda activate stmod
cd $GIT_STMOD/src/Mpipelines
# A minimum of 3GB memory (so 2 CPU per task) is required to load observed catalogues
# all below worked and finished :
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 2 --mem 6gb  -t  10 sbatch.cmd python GAL_pipeline.py z0p00 FullSky # <1 min, 2.9GB, CPU 3%
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 2 --mem 6gb  -t  10 sbatch.cmd python GAL_pipeline.py z0p02 FullSky # <1 min, 2.6GB, CPU 20%
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 2 --mem 6gb  -t  10 sbatch.cmd python GAL_pipeline.py z0p05 FullSky # 2 min, 2.7GB, CPU 21%
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 2 --mem 6gb  -t  30 sbatch.cmd python GAL_pipeline.py z0p09 FullSky # 12 min, 3.4GB, CPU 24%
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 2 --mem 6gb  -t  60 sbatch.cmd python GAL_pipeline.py z0p14 FullSky # 33 min, 4.7GB, CPU 24%
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 4 --mem 12gb -t  80 sbatch.cmd python GAL_pipeline.py z0p19 FullSky # 60 min, 7.1GB, CPU 24%
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 4 --mem 12gb -t 100 sbatch.cmd python GAL_pipeline.py z0p25 FullSky # 76 min, 9.5GB, CPU 24%
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 100 sbatch.cmd python GAL_pipeline.py z0p30 FullSky #   min,  GB, CPU
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z0p36 FullSky #   min,  GB, CPU
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z0p43 FullSky #   min,  GB, CPU
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z0p49 FullSky #   min,  GB, CPU
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z0p56 FullSky #   min,  GB, CPU
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z0p63 FullSky #   min,  GB, CPU
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z0p70 FullSky #   min,  GB, CPU
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z0p78 FullSky #   min,  GB, CPU
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z0p86 FullSky #   min,  GB, CPU
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z0p94 FullSky #   min,  GB, CPU
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z1p03 FullSky #   min,  GB, CPU
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z1p12 FullSky #   min,  GB, CPU
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z1p22 FullSky #   min,  GB, CPU
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z1p32 FullSky #   min,  GB, CPU
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z1p43 FullSky #   min,  GB, CPU
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z1p54 FullSky #   min,  GB, CPU
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z1p65 FullSky #   min,  GB, CPU
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z1p77 FullSky #   min,  GB, CPU
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z1p90 FullSky #   min,  GB, CPU
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z2p03 FullSky #   min,  GB, CPU
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z2p17 FullSky #   min,  GB, CPU
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z2p31 FullSky #   min,  GB, CPU
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z2p46 FullSky #   min,  GB, CPU
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z2p62 FullSky #   min,  GB, CPU
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z2p78 FullSky #   min,  GB, CPU
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z2p95 FullSky #   min,  GB, CPU
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z3p13 FullSky #   min,  GB, CPU
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z3p32 FullSky #   min,  GB, CPU
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z3p61 FullSky #   min,  GB, CPU
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z3p93 FullSky #   min,  GB, CPU
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z4p27 FullSky #   min,  GB, CPU
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z4p63 FullSky #   min,  GB, CPU
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z5p15 FullSky #   min,  GB, CPU
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z5p73 FullSky #   min,  GB, CPU


sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 2 --mem 6gb  -t  10 sbatch.cmd python GAL_pipeline.py z0p00 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 2 --mem 6gb  -t  10 sbatch.cmd python GAL_pipeline.py z0p02 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 2 --mem 6gb  -t  10 sbatch.cmd python GAL_pipeline.py z0p05 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 2 --mem 6gb  -t  30 sbatch.cmd python GAL_pipeline.py z0p09 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 2 --mem 6gb  -t  60 sbatch.cmd python GAL_pipeline.py z0p14 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 4 --mem 12gb -t  80 sbatch.cmd python GAL_pipeline.py z0p19 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 4 --mem 12gb -t 100 sbatch.cmd python GAL_pipeline.py z0p25 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 100 sbatch.cmd python GAL_pipeline.py z0p30 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z0p36 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z0p43 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z0p49 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z0p56 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z0p63 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z0p70 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z0p78 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z0p86 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z0p94 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z1p03 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z1p12 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z1p22 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z1p32 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z1p43 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z1p54 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z1p65 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z1p77 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z1p90 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z2p03 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z2p17 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z2p31 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z2p46 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z2p62 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z2p78 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z2p95 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z3p13 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z3p32 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z3p61 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z3p93 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z4p27 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z4p63 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z5p15 LC0002
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z5p73 LC0002

sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 2 --mem 6gb  -t  10 sbatch.cmd python GAL_pipeline.py z0p00 LC0060
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 2 --mem 6gb  -t  10 sbatch.cmd python GAL_pipeline.py z0p02 LC0060
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 2 --mem 6gb  -t  10 sbatch.cmd python GAL_pipeline.py z0p05 LC0060
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 2 --mem 6gb  -t  30 sbatch.cmd python GAL_pipeline.py z0p09 LC0060
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 2 --mem 6gb  -t  60 sbatch.cmd python GAL_pipeline.py z0p14 LC0060
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 4 --mem 12gb -t  80 sbatch.cmd python GAL_pipeline.py z0p19 LC0060
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 4 --mem 12gb -t 100 sbatch.cmd python GAL_pipeline.py z0p25 LC0060
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 100 sbatch.cmd python GAL_pipeline.py z0p30 LC0060
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z0p36 LC0060
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z0p43 LC0060
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z0p49 LC0060
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z0p56 LC0060
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z0p63 LC0060
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z0p70 LC0060
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z0p78 LC0060
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z0p86 LC0060
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z0p94 LC0060
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z1p03 LC0060
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z1p12 LC0060
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z1p22 LC0060
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z1p32 LC0060
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z1p43 LC0060
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z1p54 LC0060
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z1p65 LC0060
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z1p77 LC0060
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z1p90 LC0060
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z2p03 LC0060
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z2p17 LC0060
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z2p31 LC0060
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z2p46 LC0060
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z2p62 LC0060
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z2p78 LC0060
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z2p95 LC0060
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z3p13 LC0060
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z3p32 LC0060
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z3p61 LC0060
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z3p93 LC0060

sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 2 --mem 6gb  -t  10 sbatch.cmd python GAL_pipeline.py z0p00 LC1800
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 2 --mem 6gb  -t  10 sbatch.cmd python GAL_pipeline.py z0p02 LC1800
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 2 --mem 6gb  -t  10 sbatch.cmd python GAL_pipeline.py z0p05 LC1800
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 2 --mem 6gb  -t  30 sbatch.cmd python GAL_pipeline.py z0p09 LC1800
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 2 --mem 6gb  -t  60 sbatch.cmd python GAL_pipeline.py z0p14 LC1800
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 4 --mem 12gb -t  80 sbatch.cmd python GAL_pipeline.py z0p19 LC1800
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 4 --mem 12gb -t 100 sbatch.cmd python GAL_pipeline.py z0p25 LC1800
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 100 sbatch.cmd python GAL_pipeline.py z0p30 LC1800
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z0p36 LC1800
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z0p43 LC1800
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z0p49 LC1800
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z0p56 LC1800
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z0p63 LC1800
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z0p70 LC1800
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z0p78 LC1800
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z0p86 LC1800
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z0p94 LC1800
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z1p03 LC1800
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z1p12 LC1800
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z1p22 LC1800
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z1p32 LC1800
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z1p43 LC1800
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z1p54 LC1800
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z1p65 LC1800
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z1p77 LC1800
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z1p90 LC1800
sbatch -N 1 --ntasks-per-node 1 --cpus-per-task 8 --mem 24gb -t 120 sbatch.cmd python GAL_pipeline.py z2p03 LC1800
