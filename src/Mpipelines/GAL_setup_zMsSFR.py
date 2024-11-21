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
import calc_kcor as kk

zs = np.arange(0.0000001, 7.5, 0.001)
dm_itp = interp1d(zs, cosmo.distmod(zs).value)
# cosmology setup
z_array = np.arange(0.0001, 7.5, 0.001)
d_C = cosmo.comoving_distance(z_array)
dc_mpc = (d_C).value
dc_interpolation = interp1d(z_array, dc_mpc)
z_interpolation = interp1d(dc_mpc, z_array)
DM2_itp = interp1d(z_array, cosmo.distmod(z_array).value)

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import DotProduct, WhiteKernel

fig_dir = os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/model_GAL/GP_cal')
os.system('mkdir -p '+fig_dir)
model_dir = os.path.join(os.environ['GIT_STMOD_DATA'], 'data', 'models','model_GAL')
os.system('mkdir -p '+model_dir)
print('figures are written here : ', fig_dir)
print('model files are written here : ', model_dir)

#path_2_SDSS = os.path.join(os.environ['DATA'], 'mpecl/comparat/data_s4/galaxy_catalogues/Ti20_SDSS_kdgroups','Ti20_full_MERGED.fits')
#path_2_KIDS = os.path.join(os.environ['DATA'], 'GAMA', 'G09.GAMADR4+LegacyDR9.galreference+RM.fits')
path_2_GAMA = os.path.join(os.environ['DATA'], 'GAMA', 'forJohan.fits')
path_2_COSMOS = os.path.join(os.environ['DATA'], 'COSMOS', 'photoz_vers2.0_010312.fits')
path_2_LS10 = os.path.join(os.environ['LSDR10'], 'sweep', 'MergeALL_BGSlike_LPH.fits' )
path_2_KIDS = os.path.join(os.environ['DATA'], 'GAMA', 'G09.GAMADR4+LegacyDR9.galreference+RM.fits')

path_2_LS10_RSBC = path_2_LS10[:-5]+'_zlt0p6_RSBC.fits'
path_2_GAMA_RSBC = path_2_GAMA[:-5]+'_zlt0p6_RSBC.fits'
path_2_KIDS_RSBC = path_2_KIDS[:-5]+'_zlt0p94_RSBC.fits'
path_2_COSMOS_RSBC = path_2_COSMOS[:-5]+'_zlt0p94_RSBC.fits'

def norm_var(var):
    return (var-var.min())/(var.max()-var.min())

def norm_var_v(var, v_min, v_max):
    return (var-v_min)/(v_max-v_min)

def de_norm_var(var, min_val, max_val):
    return var*(max_val-min_val)+min_val


def get_Mabs_RS(data, str_redshift='z'):#):
    #zeropoint = 22.5
    #data['g_mag']  = zeropoint - 2.5 * np.log10(data['FLUX_G'] /data['MW_TRANSMISSION_G'])
    #data['r_mag']  = zeropoint - 2.5 * np.log10(data['FLUX_R'] /data['MW_TRANSMISSION_R'])
    #data['i_mag']  = zeropoint - 2.5 * np.log10(data['FLUX_I'] /data['MW_TRANSMISSION_I'])
    #data['z_mag']  = zeropoint - 2.5 * np.log10(data['FLUX_Z'] /data['MW_TRANSMISSION_Z'])
    DM2 = DM2_itp(data[str_redshift])
    #data['DM'] = DM2
    #kcorr_value_g = kk.calc_kcor('g', data[str_redshift], 'g - r', data['g_mag']-data['r_mag'])
    #data['mag_abs_g'] = data['g_mag'] - DM2 - kcorr_value_g
    kcorr_value_r = kk.calc_kcor('r', data[str_redshift], 'g - r', data['g_mag']-data['r_mag'])
    data['mag_abs_r'] = data['r_mag'] - DM2 - kcorr_value_r
    #kcorr_value_z = kk.calc_kcor('z', data[str_redshift], 'g - z', data['g_mag']-data['z_mag'])
    #data['mag_abs_z'] = data['z_mag'] - DM2 - kcorr_value_z

    # selection from GAMA
    #data['is_BC'] = ( data['g_mag']-data['r_mag']< - 0.055 * data['mag_abs_r']  - 0.4 )
    #data['is_GV'] = ( data['g_mag']-data['r_mag']>= - 0.055 * data['mag_abs_r']  - 0.4 ) & ( data['g_mag']-data['r_mag']<= - 0.055 * data['mag_abs_r']  - 0.35 )
    #data['is_RS'] = ( data['g_mag']-data['r_mag']> - 0.055 * data['mag_abs_r']  - 0.35 )

    # selection from redmapper
    RS_model = Table.read( os.path.join( os.environ['GIT_STMOD_DATA'], 'data', 'models', 'model_GAL', 'legacy_dr10_south_v0.3_grz_z_cal_zspec_redgals_model.fit') )
    z_RS = np.hstack(( 0., RS_model['nodes'][0] ))
    gz_RS = np.hstack(( 1.3, RS_model['meancol'][0].sum(axis=1) ))
    RS_color_gz = interp1d(z_RS, gz_RS)
    #RS_color_gz = lambda redshift : redshift * 3 + 1.3
    gz_med_RS = RS_color_gz(data[str_redshift])
    #gz_min = gz_med_RS - 2 * scat
    data['gz_med_RS'] = gz_med_RS
    data['is_RS'] = ( data['g_mag']-data['z_mag']> data['gz_med_RS'] - 0.15 )
    data['is_BC'] = ( data['g_mag']-data['z_mag']< data['gz_med_RS'] - 0.23 )
    data['is_GV'] = (~data['is_RS'])&(~data['is_BC'])
    return data

z_max_RS = 0.9499

#
# COSMOS
#
if os.path.isfile(path_2_COSMOS_RSBC):
    cosmos_RSBC = Table.read(path_2_COSMOS_RSBC)
else:
    t_cosmos = Table.read(path_2_COSMOS)
    t_cosmos['g_mag'] = t_cosmos['G']
    t_cosmos['r_mag'] = t_cosmos['R']
    t_cosmos['z_mag'] = t_cosmos['Z']
    mag_min = (t_cosmos['g_mag']>0) &(t_cosmos['r_mag']>0) &(t_cosmos['z_mag']>0)
    mag_max = (t_cosmos['g_mag']<25) &(t_cosmos['r_mag']<25) &  (t_cosmos['z_mag']<25)
    keep_cosmos = (mag_min) & (mag_max) & (t_cosmos['photoz']>0.01) & (t_cosmos['photoz']<z_max_RS)
    t_cosmos = Table(t_cosmos[keep_cosmos])
    cosmos_RSBC = get_Mabs_RS(t_cosmos, str_redshift='photoz')
    cosmos_RSBC.write(path_2_COSMOS_RSBC, overwrite = True)
    print(path_2_COSMOS_RSBC, 'written')


#
# KIDS
#
if os.path.isfile(path_2_KIDS_RSBC):
    kids_RSBC = Table.read(path_2_KIDS_RSBC)
else:
    t_kids = Table.read(path_2_KIDS)
    t_kids['g_mag'] = 8.9 - 2.5*np.log10(t_kids['flux_gt'])
    t_kids['r_mag'] = 8.9 - 2.5*np.log10(t_kids['flux_rt'])
    t_kids['z_mag'] = 8.9 - 2.5*np.log10(t_kids['flux_Zt'])
    mag_min = (t_kids['g_mag']>0) &(t_kids['r_mag']>0) &(t_kids['z_mag']>0)
    mag_max = (t_kids['g_mag']<25) &(t_kids['r_mag']<25) &  (t_kids['z_mag']<25)
    keep_kids = (mag_min) & (mag_max) & (t_kids['z_peak']>0.01) & (t_kids['z_peak']<z_max_RS)
    t_kids = Table(t_kids[keep_kids])
    kids_RSBC = get_Mabs_RS(t_kids, str_redshift='z_peak')
    kids_RSBC.write(path_2_KIDS_RSBC, overwrite = True)
    print(path_2_KIDS_RSBC, 'written')



#
# LS10
#
if os.path.isfile(path_2_LS10_RSBC):
    ls10_RSBC = Table.read(path_2_LS10_RSBC)
else:
    t_LS10 = Table.read(path_2_LS10)
    t_LS10['BEST_Z'] = t_LS10['Z_PHOT_MEAN']
    with_i = (t_LS10['Z_PHOT_MEAN_I']>-9.0)
    t_LS10['BEST_Z'][with_i] = t_LS10['Z_PHOT_MEAN_I'][with_i]
    with_S = (t_LS10['Z_SPEC']>-9.0)
    t_LS10['BEST_Z'][with_S] = t_LS10['Z_SPEC'][with_S]
    selection = (t_LS10['Z_PHOT_MEAN']<0.6)&(t_LS10['Z_PHOT_MEAN']>0.01)&(t_LS10['BEST_Z']<0.6)&(t_LS10['BEST_Z']>0.01)
    ls10 = t_LS10[selection]
    ls10_RSBC = get_Mabs_RS(ls10, str_redshift='BEST_Z')
    ls10_RSBC.write(path_2_LS10_RSBC, overwrite = True)
    print(path_2_LS10_RSBC, 'written')
#
# GAMA
#
if os.path.isfile(path_2_GAMA_RSBC):
    gama_RSBC = Table.read(path_2_GAMA_RSBC)
else:
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
    mag_min = (t_gama['u_mag']>0) &(t_gama['g_mag']>0) &(t_gama['r_mag']>0) & (t_gama['i_mag']>0) & (t_gama['z_mag']>0) &(t_gama['y_mag']>0) &(t_gama['j_mag']>0) &(t_gama['h_mag']>0) &(t_gama['k_mag']>0)
    mag_max = (t_gama['u_mag']<25) &(t_gama['g_mag']<25) &(t_gama['r_mag']<25) & (t_gama['i_mag']<25) & (t_gama['z_mag']<25) &(t_gama['y_mag']<25) &(t_gama['j_mag']<25) &(t_gama['h_mag']<25) &(t_gama['k_mag']<25)
    keep_gama = (mag_min) & (mag_max) & (t_gama['Z']>0.01) & (t_gama['Z']<0.4) &(t_gama['log10Mstar']<13)&(t_gama['log10Mstar']>6)
    t_gama = Table(t_gama[keep_gama])
    #dm_values = dm_itp(t_gama['Z'].data.data)
    #t_gama['dist_mod'] = dm_values
    gama_RSBC = get_Mabs_RS(t_gama, str_redshift='Z')
    gama_RSBC.write(path_2_GAMA_RSBC, overwrite = True)
    print(path_2_GAMA_RSBC, 'written')

s_train = np.random.random(len(t_gama)) < 0.2

# input features
x_0 = norm_var(t_gama['log10Mstar'])
x_1 = norm_var(t_gama['log10sSFR'])
X = np.transpose([x_0, x_1])[s_train]
X_verif = np.transpose([ x_0, x_1 ])[~s_train]
X_all = np.transpose([x_0, x_1])

# output features
y_u = norm_var(t_gama['u_mag']-t_gama['dist_mod'])
y_g = norm_var(t_gama['g_mag']-t_gama['dist_mod'])
y_r = norm_var(t_gama['r_mag']-t_gama['dist_mod'])
y_i = norm_var(t_gama['i_mag']-t_gama['dist_mod'])
y_z = norm_var(t_gama['z_mag']-t_gama['dist_mod'])
y_y = norm_var(t_gama['y_mag']-t_gama['dist_mod'])
y_j = norm_var(t_gama['j_mag']-t_gama['dist_mod'])
y_h = norm_var(t_gama['h_mag']-t_gama['dist_mod'])
y_k = norm_var(t_gama['k_mag']-t_gama['dist_mod'])

y = np.transpose([
                    y_u,
                    y_g,
                    y_r,
                    y_i,
                    y_z,
                    y_y,
                    y_j,
                    y_h,
                    y_k
                        ]) [s_train]

y_verif = np.transpose([
                        y_u,
                        y_g,
                        y_r,
                        y_i,
                        y_z,
                        y_y,
                        y_j,
                        y_h,
                        y_k
                            ]) [~s_train]

gpr = GaussianProcessRegressor().fit(X, y)
score = gpr.score(X, y) # 0.84
t0 = time.time()
y_pred = gpr.predict(X_verif)#, return_std=True)
print(np.min(y_verif-y_pred), np.mean(y_verif-y_pred), np.median(y_verif-y_pred), np.max(y_verif-y_pred), np.std(y_verif-y_pred) )
print((time.time()-t0)/len(y_verif))

t0 = time.time()
y_all = gpr.predict(X_all)#, return_std=True)
print((time.time()-t0)/len(X_all))

t_gama['u_GP'] = de_norm_var( y_all.T[0] , np.min(t_gama['u_mag']-t_gama['dist_mod']), np.max(t_gama['u_mag']-t_gama['dist_mod']) ) + t_gama['dist_mod']
t_gama['g_GP'] = de_norm_var( y_all.T[1] , np.min(t_gama['g_mag']-t_gama['dist_mod']), np.max(t_gama['g_mag']-t_gama['dist_mod']) ) + t_gama['dist_mod']
t_gama['r_GP'] = de_norm_var( y_all.T[2] , np.min(t_gama['r_mag']-t_gama['dist_mod']), np.max(t_gama['r_mag']-t_gama['dist_mod']) ) + t_gama['dist_mod']
t_gama['i_GP'] = de_norm_var( y_all.T[3] , np.min(t_gama['i_mag']-t_gama['dist_mod']), np.max(t_gama['i_mag']-t_gama['dist_mod']) ) + t_gama['dist_mod']
t_gama['z_GP'] = de_norm_var( y_all.T[4] , np.min(t_gama['z_mag']-t_gama['dist_mod']), np.max(t_gama['z_mag']-t_gama['dist_mod']) ) + t_gama['dist_mod']
t_gama['y_GP'] = de_norm_var( y_all.T[5] , np.min(t_gama['y_mag']-t_gama['dist_mod']), np.max(t_gama['y_mag']-t_gama['dist_mod']) ) + t_gama['dist_mod']
t_gama['j_GP'] = de_norm_var( y_all.T[6] , np.min(t_gama['j_mag']-t_gama['dist_mod']), np.max(t_gama['j_mag']-t_gama['dist_mod']) ) + t_gama['dist_mod']
t_gama['h_GP'] = de_norm_var( y_all.T[7] , np.min(t_gama['h_mag']-t_gama['dist_mod']), np.max(t_gama['h_mag']-t_gama['dist_mod']) ) + t_gama['dist_mod']
t_gama['k_GP'] = de_norm_var( y_all.T[8] , np.min(t_gama['k_mag']-t_gama['dist_mod']), np.max(t_gama['k_mag']-t_gama['dist_mod']) ) + t_gama['dist_mod']

p2_out_GAMA = path_2_GAMA[:-5]+'_GPmags.fits'
t_gama.write(p2_out_GAMA, overwrite = True)
print(p2_out_GAMA)

p2_model = os.path.join(model_dir, "gal_GP_model_GAMAtraining.pkl")
from pickle import dump
with open(p2_model, "wb") as f:
    dump(gpr, f, protocol=5)

print(p2_model, 'written')


# Here you can replace pickle with joblib or cloudpickle
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

sys.exit()

gmag_verif = de_norm_var(y_verif.T[0],np.min(t_gama['g_mag']-t_gama['dist_mod']), np.max(t_gama['g_mag']-t_gama['dist_mod']))
gmag_pred  = de_norm_var(y_pred .T[0],np.min(t_gama['g_mag']-t_gama['dist_mod']), np.max(t_gama['g_mag']-t_gama['dist_mod']))

rmag_verif = de_norm_var(y_verif.T[1],np.min(t_gama['r_mag']-t_gama['dist_mod']), np.max(t_gama['r_mag']-t_gama['dist_mod']))
rmag_pred  = de_norm_var(y_pred .T[1],np.min(t_gama['r_mag']-t_gama['dist_mod']), np.max(t_gama['r_mag']-t_gama['dist_mod']))

print(np.min(gmag_verif-gmag_pred), np.mean(gmag_verif-gmag_pred), np.median(gmag_verif-gmag_pred), np.max(gmag_verif-gmag_pred), np.std(gmag_verif-gmag_pred) )
print((time.time()-t0)/len(y_verif))

plt.figure(figsize=(6,6))
# data
plt.hist( np.log10(abs(gmag_verif-gmag_pred)), bins=50)
# models
#plt.xlabel(r'$\log_{10}(M_s/M_\odot)$')
plt.xlabel(r'$\log_{10}(|g-g_{prediction}|)$')
#plt.yscale('log')
plt.legend(frameon=False, loc=3)
#plt.ylim((-26,-15))
#plt.xlim((7.5, 12.5))
plt.tight_layout()
plt.grid()
plt.savefig(os.path.join( fig_dir, "GP_gmag_hist.png") )
plt.clf()

plt.figure(figsize=(6,6))
# data
plt.hexbin(t_gama['log10Mstar'][~s_train], np.log10(abs(gmag_verif-gmag_pred)), gridsize=20)
# models
plt.xlabel(r'$\log_{10}(M_s/M_\odot)$')
plt.ylabel(r'$\log_{10}(|g-g_{prediction}|)$')
#plt.yscale('log')
plt.legend(frameon=False, loc=3)
#plt.ylim((-26,-15))
#plt.xlim((7.5, 12.5))
plt.tight_layout()
plt.grid()
plt.savefig(os.path.join( fig_dir, "GP_Mstar_gmag.png") )
plt.clf()


plt.figure(figsize=(6,6))
# data
plt.hexbin(t_gama['log10sSFR'][~s_train], np.log10(abs(gmag_verif-gmag_pred)), gridsize=20)
# models
plt.xlabel(r'$\log_{10}(SFR/[M_\odot/yr])$')
plt.ylabel(r'$\log_{10}(|g-g_{prediction}|)$')
#plt.yscale('log')
plt.legend(frameon=False, loc=3)
#plt.ylim((-26,-15))
#plt.xlim((7.5, 12.5))
plt.tight_layout()
plt.grid()
plt.savefig(os.path.join( fig_dir, "GP_SFR_gmag.png") )
plt.clf()

plt.figure(figsize=(6,6))
# data
plt.hexbin(t_gama['dist_mod'][~s_train], np.log10(abs(gmag_verif-gmag_pred)), gridsize=20)
# models
plt.xlabel(r'Distance Modulus')
plt.ylabel(r'$\log_{10}(|g-g_{prediction}|)$')
#plt.yscale('log')
plt.legend(frameon=False, loc=3)
#plt.ylim((-26,-15))
#plt.xlim((7.5, 12.5))
plt.tight_layout()
plt.grid()
plt.savefig(os.path.join( fig_dir, "GP_distmod_gmag.png") )
plt.clf()


plt.figure(figsize=(6,6))
# data
plt.hist( np.log10(abs(rmag_verif-rmag_pred)), bins=50)
# models
#plt.xlabel(r'$\log_{10}(M_s/M_\odot)$')
plt.xlabel(r'$\log_{10}(|r-r_{prediction}|)$')
#plt.yscale('log')
plt.legend(frameon=False, loc=3)
#plt.ylim((-26,-15))
#plt.xlim((7.5, 12.5))
plt.tight_layout()
plt.grid()
plt.savefig(os.path.join( fig_dir, "GP_rmag_hist.png") )
plt.clf()

plt.figure(figsize=(6,6))
# data
plt.hexbin(t_gama['log10Mstar'][~s_train], np.log10(abs(rmag_verif-rmag_pred)), gridsize=20)
# models
plt.xlabel(r'$\log_{10}(M_s/M_\odot)$')
plt.ylabel(r'$\log_{10}(|r-r_{prediction}|)$')
#plt.yscale('log')
plt.legend(frameon=False, loc=3)
#plt.ylim((-26,-15))
#plt.xlim((7.5, 12.5))
plt.tight_layout()
plt.grid()
plt.savefig(os.path.join( fig_dir, "GP_Mstar_rmag.png") )
plt.clf()


plt.figure(figsize=(6,6))
# data
plt.hexbin(t_gama['log10sSFR'][~s_train], np.log10(abs(rmag_verif-rmag_pred)), gridsize=20)
# models
plt.xlabel(r'$\log_{10}(SFR/[M_\odot/yr])$')
plt.ylabel(r'$\log_{10}(|r-r_{prediction}|)$')
plt.legend(frameon=False, loc=3)
#plt.ylim((-26,-15))
#plt.xlim((7.5, 12.5))
plt.tight_layout()
plt.grid()
plt.savefig(os.path.join( fig_dir, "GP_SFR_rmag.png") )
plt.clf()

plt.figure(figsize=(6,6))
# data
plt.hexbin(t_gama['dist_mod'][~s_train], np.log10(abs(rmag_verif-rmag_pred)), gridsize=20)
# models
plt.xlabel(r'Distance Modulus')
plt.ylabel(r'$\log_{10}(|r-r_{prediction}|)$')
#plt.yscale('log')
plt.legend(frameon=False, loc=3)
#plt.ylim((-26,-15))
#plt.xlim((7.5, 12.5))
plt.tight_layout()
plt.grid()
plt.savefig(os.path.join( fig_dir, "GP_distmod_rmag.png") )
plt.clf()


#p_a = np.transpose(params)
#p_b = np.transpose(resid)
#OUT = np.transpose([ zall[:-1], zall[1:],  p_a[0], p_b[0], p_b[1] ])
#np.savetxt(os.path.join( model_dir, 'logMs-MK-look-up-table_2023Aug22.txt'), OUT)

#p_a = np.transpose(params2)
#p_b = np.transpose(resid2)
#OUT = np.transpose([zall[:-1], zall[1:],  p_a[0], p_a[1], p_b[0], p_b[1] ])
#np.savetxt(os.path.join( model_dir, 'logMs-MK-look-up-table-2parameters_2023Aug22.txt'), OUT)
