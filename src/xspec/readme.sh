#!/bin/bash

# tabulates LX - flux conversion as a function redshift, temperature and nH
# outputs in dir_2_result = os.path.join(os.environ['GIT_STMOD'], 'data/models/model_GAS', 'xray_k_correction')
cd $GIT_STMOD/src/xspec
python cluster_tabulate_spectra.py


