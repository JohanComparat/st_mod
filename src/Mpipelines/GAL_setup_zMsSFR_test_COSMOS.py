"""
What it does
------------

Setup of the GAL pipeline

Trains the Gaussian Proces to predict magnitudes

"""
import time
t0 = time.time()

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as plt

import sys, os
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
t_gama = Table.read(path_2_GAMA)

t_gama['u_mag'] = 8.9 - 2.5*np.log10(t_gama['flux_ut'])
t_gama['g_mag'] = 8.9 - 2.5*np.log10(t_gama['flux_gt'])
t_gama['r_mag'] = 8.9 - 2.5*np.log10(t_gama['flux_rt'])
t_gama['i_mag'] = 8.9 - 2.5*np.log10(t_gama['flux_it'])
t_gama['z_mag'] = 8.9 - 2.5*np.log10(t_gama['flux_Zt'])
t_gama['y_mag'] = 8.9 - 2.5*np.log10(t_gama['flux_Yt'])
t_gama['j_mag'] = 8.9 - 2.5*np.log10(t_gama['flux_Jt'])
t_gama['h_mag'] = 8.9 - 2.5*np.log10(t_gama['flux_Ht'])
t_gama['k_mag'] = 8.9 - 2.5*np.log10(t_gama['flux_Kt'])

t_gama['log10sSFR'] = np.log10(t_gama['SFR_50']/t_gama['StellarMass_50'])
t_gama['log10sSFR'][t_gama['log10sSFR']<-15] = -15 # clips the lowest 2.5%
t_gama['log10Mstar'] = np.log10(t_gama['StellarMass_50'])
#t_gama['Z']

mag_min = (t_gama['u_mag']>0) &(t_gama['g_mag']>0) &(t_gama['r_mag']>0) & (t_gama['i_mag']>0) & (t_gama['z_mag']>0) &(t_gama['y_mag']>0) &(t_gama['j_mag']>0) &(t_gama['h_mag']>0) &(t_gama['k_mag']>0)
mag_max = (t_gama['u_mag']<25) &(t_gama['g_mag']<25) &(t_gama['r_mag']<25) & (t_gama['i_mag']<25) & (t_gama['z_mag']<25) &(t_gama['y_mag']<25) &(t_gama['j_mag']<25) &(t_gama['h_mag']<25) &(t_gama['k_mag']<25)

keep_gama = (mag_min) & (mag_max) & (t_gama['Z']>0.01) & (t_gama['Z']<0.4) &(t_gama['log10Mstar']<13)&(t_gama['log10Mstar']>6)
t_gama = Table(t_gama[keep_gama])

dm_values = dm_itp(t_gama['Z'].data.data)
t_gama['dist_mod'] = dm_values

p2_model = os.path.join(model_dir, "gal_GP_model_GAMAtraining.pkl")
print(p2_model, 'opening')
from pickle import load
with open(p2_model, "rb") as f:
    gpr = load(f)

#
# COSMOS
#
t_cosmos = Table.read(path_2_COSMOS)
good = (t_cosmos['photoz']>0.05 )&( t_cosmos['photoz']< 6. ) & ( t_cosmos['K']>0 )& ( t_cosmos['MK']<0 )&( t_cosmos['MK']>-40 )&( t_cosmos['mass_med']<13 )&( t_cosmos['mass_med']>6 )
t_cosmos = t_cosmos[good]
dm_values = dm_itp(t_cosmos['photoz'].data.data)
t_cosmos['dist_mod'] = dm_values

# input features
x_0 = norm_var_v(t_cosmos['mass_med'], np.min(t_gama['log10Mstar']), np.max(t_gama['log10Mstar']))
x_1 = norm_var_v(t_cosmos['ssfr_med'], np.min(t_gama['log10sSFR'] ), np.max(t_gama['log10sSFR'] ))
X_all_cosmos = np.transpose([x_0, x_1])

sys.exit()
# slice X_all_cosmos in many small instances

t0 = time.time()
y_all_cosmos = gpr.predict(X_all_cosmos)#, return_std=True)
print((time.time()-t0)/len(X_all_cosmos))

t_cosmos['u_GP'] = de_norm_var( y_all_cosmos.T[0] , np.min(t_gama['u_mag']-t_gama['dist_mod']), np.max(t_gama['u_mag']-t_gama['dist_mod']) ) + t_cosmos['dist_mod']
t_cosmos['g_GP'] = de_norm_var( y_all_cosmos.T[1] , np.min(t_gama['g_mag']-t_gama['dist_mod']), np.max(t_gama['g_mag']-t_gama['dist_mod']) ) + t_cosmos['dist_mod']
t_cosmos['r_GP'] = de_norm_var( y_all_cosmos.T[2] , np.min(t_gama['r_mag']-t_gama['dist_mod']), np.max(t_gama['r_mag']-t_gama['dist_mod']) ) + t_cosmos['dist_mod']
t_cosmos['i_GP'] = de_norm_var( y_all_cosmos.T[3] , np.min(t_gama['i_mag']-t_gama['dist_mod']), np.max(t_gama['i_mag']-t_gama['dist_mod']) ) + t_cosmos['dist_mod']
t_cosmos['z_GP'] = de_norm_var( y_all_cosmos.T[4] , np.min(t_gama['z_mag']-t_gama['dist_mod']), np.max(t_gama['z_mag']-t_gama['dist_mod']) ) + t_cosmos['dist_mod']
t_cosmos['y_GP'] = de_norm_var( y_all_cosmos.T[5] , np.min(t_gama['y_mag']-t_gama['dist_mod']), np.max(t_gama['y_mag']-t_gama['dist_mod']) ) + t_cosmos['dist_mod']
t_cosmos['j_GP'] = de_norm_var( y_all_cosmos.T[6] , np.min(t_gama['j_mag']-t_gama['dist_mod']), np.max(t_gama['j_mag']-t_gama['dist_mod']) ) + t_cosmos['dist_mod']
t_cosmos['h_GP'] = de_norm_var( y_all_cosmos.T[7] , np.min(t_gama['h_mag']-t_gama['dist_mod']), np.max(t_gama['h_mag']-t_gama['dist_mod']) ) + t_cosmos['dist_mod']
t_cosmos['k_GP'] = de_norm_var( y_all_cosmos.T[8] , np.min(t_gama['k_mag']-t_gama['dist_mod']), np.max(t_gama['k_mag']-t_gama['dist_mod']) ) + t_cosmos['dist_mod']

p2_out_COSMOS = path_2_COSMOS[:-5]+'_GPmags.fits'
t_cosmos.write(p2_out_COSMOS, overwrite = True)

