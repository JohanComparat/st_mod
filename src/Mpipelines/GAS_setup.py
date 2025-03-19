
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

C_GAS.draw_and_save_simulated_profiles_m14(p_2_out=os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/model_GAS', 'profiles_010z015_1e14M2e14.fits'), nsim = 100000, zmin = 0.1, zmax = 0.15 )
