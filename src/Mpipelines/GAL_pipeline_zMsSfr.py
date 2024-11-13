"""
What it does
------------

Setup of the GAL pipeline

Trains the Gaussian Proces to predict magnitudes

"""
import time
t0 = time.time()
import sys, os, glob
from astropy.table import Table
import numpy as np
sys.path.append( os.path.join(os.environ['GIT_STMOD'], 'src') )
from models import GAL as GG

z_dir = sys.argv[1]
LC_dir = sys.argv[2] # 'FullSky'
#z_dir = 'z0p09'
#LC_dir = 'FullSky'

C_GAL = GG.GAL(z_dir, LC_dir=LC_dir)

from scipy.interpolate import interp1d
from scipy.stats import norm
from astropy.table import Table
from scipy.optimize import curve_fit
import numpy as np

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import extinction
cosmoUCHUU = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
h = 0.6774
L_box = 1000.0 / h
cosmo = cosmoUCHUU

zs = np.arange(0.0000001, 7.1, 0.001)
dm_itp = interp1d(zs, cosmo.distmod(zs).value)

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import DotProduct, WhiteKernel

fig_dir = os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/model_GAL/GP_cal')
os.system('mkdir -p '+fig_dir)
model_dir = os.path.join(os.environ['GIT_STMOD_DATA'], 'data', 'models','model_GAL')
os.system('mkdir -p '+model_dir)
print('figures are written here : ', fig_dir)
print('model files are written here : ', model_dir)

path_2_GAMA = os.path.join(os.environ['DATA'], 'GAMA', 'forJohan.fits')
path_2_COSMOS = os.path.join(os.environ['DATA'], 'COSMOS', 'photoz_vers2.0_010312.fits')


def norm_var(var):
    return (var-var.min())/(var.max()-var.min())

def norm_var_v(var, v_min, v_max):
    return (var-v_min)/(v_max-v_min)

def de_norm_var(var, min_val, max_val):
    return var*(max_val-min_val)+min_val

#
# GAMA
#
t_COS = Table.read(path_2_COSMOS)

#t_COS['U']
#t_COS['G']
#t_COS['R']
#t_COS['I']
#t_COS['Z']
#t_COS['NUV']
#t_COS['J']
#t_COS['I1']
#t_COS['K']

t_COS['log10sSFR'] = t_COS['ssfr_med']
t_COS['log10sSFR'][t_COS['log10sSFR']<-15] = -15 # clips the lowest XX%

mag_min = (t_COS['U']>0) &(t_COS['G']>0) &(t_COS['R']>0) & (t_COS['I']>0) & (t_COS['Z']>0) &(t_COS['NUV']>0) &(t_COS['J']>0) &(t_COS['I1']>0) &(t_COS['K']>0)
mag_max = (t_COS['U']<26) &(t_COS['G']<26) &(t_COS['R']<26) & (t_COS['I']<26) & (t_COS['Z']<26) &(t_COS['NUV']<26) &(t_COS['J']<26) &(t_COS['I1']<26) &(t_COS['K']<26)

keep_COS = (mag_min) & (mag_max) & (t_COS['photoz']>0.01) & (t_COS['photoz']<1.) &(t_COS['mass_med']<13)&(t_COS['mass_med']>7.1)
t_COS = Table(t_COS[keep_COS])

dm_values = dm_itp(t_COS['photoz'].data.data)
t_COS['dist_mod'] = dm_values

p2_model = os.path.join(model_dir, "gal_GP_model_COSMOStraining.pkl")
print(p2_model, 'opening')
from pickle import load
with open(p2_model, "rb") as f:
    gpr = load(f)

print(len(C_GAL.p_2_catalogues), 'catalogues to loop over', ', Dt=', time.time()-t0, 'seconds')
is_cat = np.array([os.path.isfile(el) for el in C_GAL.p_2_catalogues])
for p_2_catalogue in C_GAL.p_2_catalogues[is_cat]:
	print('='*100)
	p_2_catalogue_out = os.path.join( os.path.dirname(p_2_catalogue), 'zMsSFRmatch_mags.fits')
	print('reading', p_2_catalogue)
	t_sim_gal = Table.read(p_2_catalogue)
	print('N lines in catalog = ',len(t_sim_gal), ', Dt=', time.time()-t0, 'seconds')
	print(t_sim_gal.info())
	'obs_sm'
	'obs_sfr'
	'redshift_S'
	dm_values = dm_itp(t_sim_gal['redshift_S'].data.data)
	keep_cols = [
				#'x',
				#'y',
				#'z',
				#'Mvir',
				#'icl',
				#'id',
				'obs_sfr',
				'obs_sm',
				#'obs_uv',
				#'sfr',
				#'sm',
				#'upid',
				#'A_UV',
				#'Mpeak',
				#'Vmax_Mpeak',
				#'desc_id',
				#'lvmp',
				#'vmax',
				#'vx',
				#'vy',
				#'vz',
				#'RA',
				#'DEC',
				#'g_lat',
				#'g_lon',
				#'ecl_lat',
				#'ecl_lon',
				#'redshift_R',
				'redshift_S',
				#'dL',
				#'nH',
				#'ebv',
				]
	t_sim_gal.keep_columns(keep_cols)
	t_sim_gal['dist_mod'] = dm_values

	# input features
	t_sim_gal['obs_sfr'][t_sim_gal['obs_sfr']<=0]=1e-6
	x_0 = norm_var_v(np.log10(t_sim_gal['obs_sm']), np.min(t_COS['mass_med']), np.max(t_COS['mass_med']))
	x_1 = norm_var_v(np.log10(t_sim_gal['obs_sfr']/t_sim_gal['obs_sm']), np.min(t_COS['log10sSFR'] ), np.max(t_COS['log10sSFR'] ))
	X_all = np.transpose([x_0, x_1])

	t0 = time.time()
	y_all = gpr.predict(X_all)#, return_std=True)
	print((time.time()-t0)/len(X_all))


	t_sim_gal['nuv_GP'] = de_norm_var( y_all.T[0] , np.min(t_COS['NUV'] -t_COS['dist_mod']), np.max(t_COS['NUV']-t_COS['dist_mod']) ) + t_sim_gal['dist_mod']
	t_sim_gal['u_GP'  ] = de_norm_var( y_all.T[1] , np.min(t_COS['U']   -t_COS['dist_mod']), np.max(t_COS['U']  -t_COS['dist_mod']) ) + t_sim_gal['dist_mod']
	t_sim_gal['g_GP'  ] = de_norm_var( y_all.T[2] , np.min(t_COS['G']   -t_COS['dist_mod']), np.max(t_COS['G']  -t_COS['dist_mod']) ) + t_sim_gal['dist_mod']
	t_sim_gal['r_GP'  ] = de_norm_var( y_all.T[3] , np.min(t_COS['R']   -t_COS['dist_mod']), np.max(t_COS['R']  -t_COS['dist_mod']) ) + t_sim_gal['dist_mod']
	t_sim_gal['i_GP'  ] = de_norm_var( y_all.T[4] , np.min(t_COS['I']   -t_COS['dist_mod']), np.max(t_COS['I']  -t_COS['dist_mod']) ) + t_sim_gal['dist_mod']
	t_sim_gal['z_GP'  ] = de_norm_var( y_all.T[5] , np.min(t_COS['Z']   -t_COS['dist_mod']), np.max(t_COS['Z']  -t_COS['dist_mod']) ) + t_sim_gal['dist_mod']
	t_sim_gal['j_GP'  ] = de_norm_var( y_all.T[6] , np.min(t_COS['J']   -t_COS['dist_mod']), np.max(t_COS['J']  -t_COS['dist_mod']) ) + t_sim_gal['dist_mod']
	t_sim_gal['i1_GP' ] = de_norm_var( y_all.T[7] , np.min(t_COS['I1']  -t_COS['dist_mod']), np.max(t_COS['I1'] -t_COS['dist_mod']) ) + t_sim_gal['dist_mod']
	t_sim_gal['k_GP'  ] = de_norm_var( y_all.T[8] , np.min(t_COS['K']   -t_COS['dist_mod']), np.max(t_COS['K']  -t_COS['dist_mod']) ) + t_sim_gal['dist_mod']
	t_sim_gal['ext_GP'] = de_norm_var( y_all.T[9] , np.min(t_COS['ebv'])                   , np.max(t_COS['ebv'])                   )

	t_sim_gal.write(p_2_catalogue_out, overwrite = True)
	print(p_2_catalogue_out, 'written')


