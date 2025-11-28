
import time
import sys, os, glob
from astropy.table import Table
import numpy as np

#Starting time
t0 = time.time()

#Directory for code
sys.path.append(os.path.join(os.environ['GIT_STMOD'], 'src'))
from models import GAS as GG

z_dir = sys.argv[1]
LC_dir = sys.argv[2] # 'FullSky'
C_GAS = GG.GAS(z_dir, b_HS=0.8, logM500c_min=11., logFX_min=-18, LC_dir=LC_dir)

#Iterate over replications
print('There are {0} replications to loop over.\nInitializing GAS class took {1} seconds'.format(len(C_GAS.p_2_catalogues), time.time()-t0))
is_cat = np.array([os.path.isfile(el) for el in C_GAS.p_2_catalogues])
for p2ci, p_2_catalogue in enumerate(C_GAS.p_2_catalogues[is_cat]):
	t0 = time.time()
	print('='*100)
	print('Replication {0} of {1}'.format(p2ci+1, len(C_GAS.p_2_catalogues)))
	print('Reading file {0}'.format(p_2_catalogue))
	t_sim_gal = Table.read(p_2_catalogue)
	t_sim_gal['ID_glist'] = np.arange(len(t_sim_gal))
	# select distinct haloes
	t1 = t_sim_gal[(t_sim_gal['upid']==-1)&(t_sim_gal['Mvir']>=5e11)]
	C_GAS.prepare_halo_cat(t1)
	# assigns
	# C_GAS.CAT (catalog) and N_cat, zmin, zmax
	# assigns
	# simulated objects as attribute to the class
	p_2_catalogue_out = os.path.join( os.path.dirname(p_2_catalogue), 'Xgas_bHS0.8.fits')
	C_GAS.populate_cat( p_2_profiles = os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/model_GAS', 'profiles_010z015_1e14M2e14.fits'), kt_m_slope = 0.6, kt_m_intercept = -8. )
	C_GAS.CAT.write(p_2_catalogue_out, overwrite = True)
	print('Output written to {0}\nProcessing the replication took {1} seconds'.format(p_2_catalogue_out, time.time()-t0))
	C_GAS.make_simput( p_2_catalogue_out,
				   p_2_profiles = os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/model_GAS', 'profiles_010z015_1e14M2e14.fits'),
				   dir_2_simput = os.path.join( os.path.dirname(p_2_catalogue) ),
				   simput_file_name =  'Xgas_bHS0.8_simput.fits'
				   )