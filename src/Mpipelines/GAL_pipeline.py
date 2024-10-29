
import time
t0 = time.time()
import sys, os, glob
from astropy.table import Table
import numpy as np
sys.path.append( os.path.join(os.environ['GIT_STMOD'], 'src') )
from models import GAL as GG

z_dir = sys.argv[1]
LC_dir = sys.argv[2] # 'FullSky'
C_GAL = GG.GAL(z_dir, LC_dir=LC_dir)
C_GAL.construct_OBS_kdTree()

print(len(C_GAL.p_2_catalogues), 'catalogues to loop over', ', Dt=', time.time()-t0, 'seconds')
is_cat = np.array([os.path.isfile(el) for el in C_GAL.p_2_catalogues])
for p_2_catalogue in C_GAL.p_2_catalogues[is_cat]:
	print('='*100)

	p_2_catalogue_out = os.path.join( os.path.dirname(p_2_catalogue), 'Kmatch_mags.fits')

	#if os.path.isfile(p_2_catalogue_out)==False:

	print('reading', p_2_catalogue)
	t_sim_gal = Table.read(p_2_catalogue)
	print('N lines in catalog = ',len(t_sim_gal), ', Dt=', time.time()-t0, 'seconds')

	t_sim_gal = C_GAL.add_kmag(t_sim_gal)
	print(len(t_sim_gal),' K mag abs added, Dt=', time.time()-t0, 'seconds')

	t_out = C_GAL.add_magnitudes_direct_match_K(t_sim_gal)
	print( len(t_out),' matched magnitudes from ref cat added, Dt=', time.time()-t0, 'seconds')

	t_out.write(p_2_catalogue_out, overwrite = True)
	print(p_2_catalogue_out, 'written, t=', time.time()-t0)

