# stack_models

## General statements

This repository contains the models used to fit X-ray stacks around galaxies. 

Seven distinct model modules :
 1. [PSF] PSF profile 
 2. [BG] Background surface brightness level
 3. [SE] Selection effects via the stellar to halo mass relation
 4. [SATC] contamination by satellite galaxies in haloes
 5. [PS] emission for point sources ( PS_AGN, PS_XRB )
 6. [GAS] emission from hot gas
 7. [GAL] galaxies' SED for mock catalogues

# Setup and Data

In your bashrc, add two environment variables linking to where the code is stored (GIT_STMOD) and to where the data to interpret is stored (DATA_S4). In addition link to the Uchuu light cone repository (UCHUU)
```
export GIT_STMOD='/home/${USERNAME}/software/st_mod'
export GIT_STMOD_DATA='/home/${USERNAME}/software/st_mod_data'
export DATA_S4='/home/${USERNAME}/data/data_s4'
export UCHUU='/home/${USERNAME}/data/Uchuu'
```

$UCHUU should contain the 'GPX8' directory (also linked 'FullSky') with light cone shells z?p?? and their content (replication_*). The complete data set for the full sky light cone until redshift 6 is 18 TB.

In DATA_S4 stack results from Zhang et al. 2024a,b are available
It contains three directories :
 * mergedCubes : with the tabulated profiles and spectra
 * galaxy_catalogues : with the galaxy catalogues
 * simulated_catalogues : with the simulated mock catalogues

In each directory, different stacking experiments correspond to different folders :
 * Ti20_SDSS_kdgroups : spec-Z stack using the Tinker 2022 SDSS catalogue
 * DR7Q_SDSS_quasar : SDSS DR7 QSO sample (Richards et al. 2002)
 * DR9_Ebins_s4 : legacy survey DR9 sample (full and isolated)
 * SGA_Ebins : edge-on and face-on extracted legacy survey DR9 sample and SGA
 * GAMA : spec-Z stacks on eFEDS and eRASS
 * SDSS : spec-Z stacks with SDSS on eRASS
 * compilation : on the spec-Z compilation

# Documentation

To generate an automated sphinx-style documentation (web pages), execute :
```
cd $GIT_STMOD
sphinx-build -M html docs/source/ docs/build/
```

# Light cones

Light-weight light cones with full redshift reach to iterate on model fitting
 * $UCHUU/LC0002 : 2 square degrees e.g. COSMOS selected by 149.2<ra<151.1 && 1.3<Dec<3.1 up to redshift 6. Size (glist.fits files): 1.2GB
 * $UCHUU/LC0060 : 60 square degrees e.g. GAMA 9h selected by 129<ra<141 && -2<Dec<3      up to redshift 6. Size (glist.fits files) : 16G
 * $UCHUU/LC1800 : 1800 square degrees e.g. fat equatorial stripe (HSC Wide) selected by abs(dec)<2.5 & All R.A. Size (glist.fits files): 243GB
 * $UCHUU/FullSky : full sky. Size (all files) 18TB.

Each cone covers approx 30 times more area than the previous

Metadata about the area covered by each lightcone piece are in these files : $UCHUU/area_per_replica_z?p??.fits

## extraction of small light cones

Using: extract_small_LC.py, recreate LC0002, LC0060, LC1800 by running :
```
cd $GIT_STMOD/src/lightcones
sh run_small_LC.sh
```

LC0002, LC0060, LC1800, FullSky, are generic arguments (LC_dir) of the model classes to know on which lightcone they should be applied.

# conda setup and packages

Create the anaconda environment 'stmod' for stack models

First, update your conda :
```conda update -n base -c defaults conda ```

Then  ```cd $GIT_STMOD && conda env create -f environment.yml ```

To use the environment do:  ```conda activate stmod ``` (or  ```source activate stmod ``` in a sbatch command)

# Models

The models are available in the directory :
 * $GIT_STMOD/src/models/
  * AGN.py : class to handle the AGN model
  * XRB.py : class to handle the XRB model
  * GAS.py    : class to handle the GAS model
  * GAL.py    : class to handle the GAL model
  * SATC.py   : class to handle the SAT model
  * PSF.py    : class to handle the PSF model
  * SE.py     : class to handle the SE model
  * BG.py     : class to handle the BG model

Input files for the models are stored here :
 * $GIT_STMOD_DATA/data/models/model_AGN
 * $GIT_STMOD_DATA/data/models/model_GAS
 * $GIT_STMOD_DATA/data/models/model_GAL

Model Pipelines. The set of scripts to run the models on simulated data are here : $GIT_STMOD/src/Mpipelines

### Models for GAL, GAS and AGN based on Uchuu/Rockstar/UniverseMachine

On the Uchuu simulated light cone, we paint X-ray emission (Comparat et al. 2019, 2020) on the galaxy population modelled with UniverseMachine (Behroozi et al. 2019).
A calculation on the full sky up to redshift 6 takes few hours, so an MCMC on the parameters going in is not feasible.

On the smaller lightcones and for a limited redshift reach the light cone can be loaded in memory and model computed on the fly enabling parameter optimization.

### Validation data sets

The models' predictions need to agree with current measurements of the luminosity functions, logNlogS, etc
The validation data (for comparison) are stored here :
 * $GIT_STMOD_DATA/data/validation/validation_AGN
 * $GIT_STMOD_DATA/data/validation/validation_GAS
 * $GIT_STMOD_DATA/data/validation/validation_GAL


# GAL model

### setup the GAL model
Creates the look up table for the K mag absolute magnitude -- stellar mass -- redshift relation
The setup creates input data for the model recorded here :
 * $GIT_STMOD_DATA/data/models/model_GAL
```
cd $GIT_STMOD/src/Mpipelines/
python GAL_setup.py
```
### Run the GAL model
Submit all computations to the slurm sbatch system
```
cd $GIT_STMOD/src/Mpipelines/
sh run_GAL_pipelines.sh
```

### Validate the GAL model

Creates a set of fiducial comparisons (figures) to validate the model computed. It currently outputs the stellar mass function and the K-band luminosity function. Output are here : $GIT_STMOD_DATA/data/validation/validation_GAL
```
cd $GIT_STMOD/src/Mpipelines/
nohup sh run_GAL_validation.sh > GAL_logs/run_GAL_validation.log &
```
### Next developments

Next dev :
 * interface with sbatch to ease submission of jobs (match memory requirement accurately) (Mpipelines/run_GAL_pipelines.sh and sbatch.cmd)
 * split the populations into star foming and quiescent galaxies to assign SED (GAL.py).
 * more validation (other luminosity functions: r band, i band, ...) and even more when galaxy populations will be split into red sequence and blue cloud (GAL_validation.py).


# GAS model

The GAS model is taken from Comparat, Eckert et al. 2019, then extended by Seppi, Comparat et al. 2022.
External dependences are included in this package here :
 * $GIT_STMOD_DATA/data/models/model_GAS

New version of scaling relation sampling implemented to make catalogues.
Todo : tabulate profiles and sets of images


## Run the GAS model
Submit all computations to the slurm sbatch system
```
cd $GIT_STMOD/src/Mpipelines/
nohup sh run_GAS_pipeline.sh > logs/run_GAS_pipelines.log &
```
The run is finished for b=0.8 

ds43, two screens
xpsec K-correction ongoing
cd $GIT_STMOD/src/Mpipelines/
python cluster_tabulate_spectra.py
Then push tabulated results on the git
cd $GIT_STMOD
git add data/models/model_GAS/xray_k_correction
git commit -m"Xspec K-correction for APEC model"
git push


## Validate the GAS model

Creates a set of fiducial comparisons (figures) to validate the model computed.
It outputs per redshift slice :
 * scaling relation M500c - LX 500c
 * scaling relation stellar mass - LX 500c
It should also output for a full light cone
 * logN logS

Outputs are here : Output are here : $GIT_STMOD_DATA/data/validation/validation_GAS

```
cd $GIT_STMOD/src/Mpipelines/
nohup sh run_GAS_validation.sh > GAL_logs/run_GAS_validation.log &
```

### Next developments


Next dev to improve on the model
 * split red sequence - blue cloud predictions ?
 * for low masses based on the stacking results


# Active galactic nuclei model (PS_AGN)

## setup the AGN model

Tabulates the AGN model from Comparat et al. 2019 in redshift bins of dz=0.01 until redshift 6 in that folder : $UCHUU/AGN_LX_tables (600 files, 262 GB)

```
cd $GIT_STMOD/src/Mpipelines/
sh run_AGN_setup.sh
```

## Run the AGN model

Abundance matching procedure with the simulated files to create AGN catalogues.

```
cd $GIT_STMOD/src/Mpipelines/
sh run_AGN_pipeline.sh
```
```
cd $GIT_STMOD/src/Mpipelines/

nohup sh run_AGN_pipeline.sh                 > AGN_logs/run_AGN_pipeline.log &            # DONE
nohup sh run_AGN_pipeline_sigma0p4_fsat0.sh  > AGN_logs/run_AGN_pipeline_sigma0p4_fsat0 & # DONE
nohup sh run_AGN_pipeline_sigma0p4_fsat8.sh  > AGN_logs/run_AGN_pipeline_sigma0p4_fsat8 & # DONE
nohup sh run_AGN_pipeline_sigma0p6_fsat0.sh  > AGN_logs/run_AGN_pipeline_sigma0p6_fsat0 & # DONE
# 06 :
nohup sh run_AGN_pipeline_sigma0p6_fsat8.sh  > AGN_logs/run_AGN_pipeline_sigma0p6_fsat8 &
nohup sh run_AGN_pipeline_sigma0p8_fsat0.sh  > AGN_logs/run_AGN_pipeline_sigma0p8_fsat0 & # DONE
nohup sh run_AGN_pipeline_sigma0p8_fsat8.sh  > AGN_logs/run_AGN_pipeline_sigma0p8_fsat8 & # DONE
nohup sh run_AGN_pipeline_sigma0p8_fsat20.sh  > AGN_logs/run_AGN_pipeline_sigma0p8_fsat20 & # DONE
```

## Validate the AGN model


```
cd $GIT_STMOD/src/Mpipelines/
nohup sh run_AGN_validation.sh > AGN_logs/run_AGN_validation.log &
```

### Next developments

Need to create benchmark with the stacking results :
 * with and wthout masking
 * SF / QU

Stack of all CEN galaxies with M* bin split in SF/QU
 * remove XRBs
 * remove CGM (beta profile)
 * infer AGN LX


# X-ray Binaries (PS_XRB)

They are purely analytical and based on Aird et al. (2017). They rely on our knowledge of stellar mass and star formation rate of the samples.
They can be applied directly to the data samples or to the simulated equivalent.

Tables with parameters of the models are stored here :
 * $GIT_STMOD_DATA/data/models/model_XRB

# SE models

The SE model leverage the mock catalogues to infer how the sample studied samples the stellar mass function, the halo mass function and its relation : stellar to halo mass relation SMHMR.

# BG models

Inference of the background level

# SATC models

Satellite contamination of the stacked profiles




# Reading the data (X-ray stacked profiles and spectra)

All the functions to read outputs from the stacking analysis are available in the library :
 * $GIT_STMOD/src/io/lib_read_stacks.py
  * Profile.py : class to handle X-ray surface brightness profiles
  * Spectrum.py : class to handle X-ray spectra

## Visualization of the data

the readme.sh file explain what each plotting script outputs.

# Comparing Models and Data

Set of scripts to compare observations and theory : $GIT_STMOD/src/fit

##
