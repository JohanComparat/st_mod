#!/bin/bash -l
#
# Single-core example job script for MPCDF Raven.
# In addition to the Python example show here, the script
# is valid for any single-threaded program, including
# sequential Matlab, Mathematica, Julia, and similar cases.
#
#SBATCH -J stmGAL
#SBATCH -o ./GAL_logs/out.%j
#SBATCH -e ./GAL_logs/err.%j
#SBATCH -D ./

# Set number of OMP threads to fit the number of available cpus, if applicable.
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

ARGS="$@" && echo "ARGS $ARGS"

source activate stmod
cd $GIT_STMOD/src/Mpipelines
# Run single-core program
srun $ARGS
