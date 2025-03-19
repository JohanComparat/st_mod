
import time
t0 = time.time()
import sys, os, glob
from astropy.table import Table
import numpy as np
sys.path.append( os.path.join(os.environ['GIT_STMOD'], 'src') )
from models import GAS as GG

z_dir = sys.argv[1]
LC_dir = sys.argv[2] # 'FullSky'
C_GAS = GG.GAS(z_dir, b_HS=0.8, logM500c_min=11., logFX_min=-18, LC_dir=LC_dir)

print(len(C_GAS.p_2_catalogues), 'catalogues to loop over', ', Dt=', time.time()-t0, 'seconds')
is_cat = np.array([os.path.isfile(el) for el in C_GAS.p_2_catalogues])
for p_2_catalogue in C_GAS.p_2_catalogues[is_cat]:
	print('='*100)
	#if os.path.isfile(p_2_catalogue_out)==False:
	print('reading', p_2_catalogue)
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
	C_GAS.populate_cat( p_2_profiles = os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/model_GAS', 'profiles_010z015_1e14M2e14.fits') )
	C_GAS.CAT.write(p_2_catalogue_out, overwrite = True)
	print(p_2_catalogue_out, 'written, t=', time.time()-t0)
	C_GAS.make_simput( p_2_catalogue_out,
				   p_2_profiles = os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/model_GAS', 'profiles_010z015_1e14M2e14.fits'),
				   dir_2_simput = os.path.join( os.path.dirname(p_2_catalogue), 'simput')
				   simput_file_name =  'Xgas_bHS0.8_simput.fits'
				   )

	#p_2_catalogue_out = os.path.join( os.path.dirname(p_2_catalogue), 'Xgas_bHS0.6.fits')
	#C_GAS.b_HS = 0.6
	#C_GAS.populate_cat()
	#C_GAS.CAT.write(p_2_catalogue_out, overwrite = True)
	#print(p_2_catalogue_out, 'written, t=', time.time()-t0)

	#p_2_catalogue_out = os.path.join( os.path.dirname(p_2_catalogue), 'Xgas_bHS1.0.fits')
	#C_GAS.b_HS = 1.0
	#C_GAS.populate_cat()
	#C_GAS.CAT.write(p_2_catalogue_out, overwrite = True)
	#print(p_2_catalogue_out, 'written, t=', time.time()-t0)

