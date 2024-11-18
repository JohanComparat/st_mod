
import time
t0 = time.time()
import sys, os, glob
from astropy.table import Table, vstack
import numpy as np
sys.path.append( os.path.join(os.environ['GIT_STMOD'], 'src') )
from models import GAL as GG

z_dir = sys.argv[1]
LC_dir = sys.argv[2] # 'FullSky'
C_GAL = GG.GAL(z_dir, LC_dir=LC_dir)
C_GAL.construct_OBS_kdTree_RSBC()

print(len(C_GAL.p_2_catalogues), 'catalogues to loop over', ', Dt=', time.time()-t0, 'seconds')
is_cat = np.array([os.path.isfile(el) for el in C_GAL.p_2_catalogues])
for p_2_catalogue in C_GAL.p_2_catalogues[is_cat]:
	print('='*100)

	p_2_catalogue_out = os.path.join( os.path.dirname(p_2_catalogue), 'zMsSFRmatch_mags.fits')

	#if os.path.isfile(p_2_catalogue_out)==False:

	print('reading', p_2_catalogue)
	t_sim_gal = Table.read(p_2_catalogue)
	print('N lines in catalog = ',len(t_sim_gal), ', Dt=', time.time()-t0, 'seconds')
	# split RS, BC
	ssfr = np.log10(t_sim_gal['obs_sfr']/t_sim_gal['obs_sm'])
	sim_qu = ssfr<=-11
	sim_sf = ssfr>-11
	print('N QU = ',len(t_sim_gal[sim_qu]))
	print('N SF = ',len(t_sim_gal[sim_sf]))

	t_sim_gal = C_GAL.add_kmag(t_sim_gal)
	print(len(t_sim_gal),' K mag abs added, Dt=', time.time()-t0, 'seconds')
	t_sim_gal['ID_glist']   = np.arange(len(t_sim_gal))

	t_out_sf = C_GAL.add_magnitudes_direct_match_K(C_GAL.Tree_Obs_SF, C_GAL.t_ref_BC, t_sim_gal[sim_sf])
	t_out_qu = C_GAL.add_magnitudes_direct_match_K(C_GAL.Tree_Obs_QU, C_GAL.t_ref_RS, t_sim_gal[sim_qu])
	t_out_sf['ID_glist'] = t_sim_gal['ID_glist'][sim_sf]
	t_out_qu['ID_glist'] = t_sim_gal['ID_glist'][sim_qu]
	print( len(t_out_qu),' matched magnitudes from QU ref cat added, Dt=', time.time()-t0, 'seconds')
	print( len(t_out_sf),' matched magnitudes from SF ref cat added, Dt=', time.time()-t0, 'seconds')
	t_out = vstack(( t_out_sf, t_out_qu ))
	t_out = t_out[np.argsort(t_out['ID_glist'])]
	t_out.write(p_2_catalogue_out, overwrite = True)
	print(p_2_catalogue_out, 'written, t=', time.time()-t0)

