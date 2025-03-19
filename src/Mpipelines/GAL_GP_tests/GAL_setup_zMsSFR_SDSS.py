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

fig_dir = os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/model_GAL/GP_cal_COSMOS')
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
path_2_SDSS = os.path.join(os.environ['DATA'], 'SDSS', 'MPA_JHU', 'full_MERGED.fits' )

def norm_var(var):
    return (var-var.min())/(var.max()-var.min())

def norm_var_v(var, v_min, v_max):
    return (var-v_min)/(v_max-v_min)

def de_norm_var(var, min_val, max_val):
    return var*(max_val-min_val)+min_val

#
# SDSS
#
t_SDSS = Table.read(path_2_SDSS)

t_SDSS['U_mag'] = t_SDSS['PLUG_MAG'].T[0]
t_SDSS['G_mag'] = t_SDSS['PLUG_MAG'].T[1]
t_SDSS['R_mag'] = t_SDSS['PLUG_MAG'].T[2]
t_SDSS['I_mag'] = t_SDSS['PLUG_MAG'].T[3]
t_SDSS['Z_mag'] = t_SDSS['PLUG_MAG'].T[4]


mag_min = (t_SDSS['U_mag']>0)  &(t_SDSS['G_mag']>0)  &(t_SDSS['R_mag']>0)  & (t_SDSS['I_mag']>0)  & (t_SDSS['Z_mag']>0)
mag_max = (t_SDSS['U_mag']<26) &(t_SDSS['G_mag']<26) &(t_SDSS['R_mag']<26) & (t_SDSS['I_mag']<26) & (t_SDSS['Z_mag']<26)

keep_SDSS = (mag_min) & (mag_max) & (t_SDSS['Z']>0.01) & (t_SDSS['Z']<1.) &(t_SDSS['LGM_TOT_P50']<13)&(t_SDSS['LGM_TOT_P50']>7.1)&(t_SDSS['SPECSFR_TOT_P50']>-20)
t_SDSS = Table(t_SDSS[keep_SDSS])

dm_values = dm_itp(t_SDSS['Z'].data.data)
t_SDSS['dist_mod'] = dm_values

s_train = np.random.random(len(t_SDSS)) < 0.01

# input features
x_0 = norm_var(t_SDSS['LGM_TOT_P50'])
x_1 = norm_var(t_SDSS['SPECSFR_TOT_P50'])
X = np.transpose([x_0, x_1])[s_train]
X_verif = np.transpose([ x_0, x_1 ])[~s_train]
X_all = np.transpose([x_0, x_1])

# output features
y_u   = norm_var(t_SDSS['U_mag'  ]-t_SDSS['dist_mod'])
y_g   = norm_var(t_SDSS['G_mag'  ]-t_SDSS['dist_mod'])
y_r   = norm_var(t_SDSS['R_mag'  ]-t_SDSS['dist_mod'])
y_i   = norm_var(t_SDSS['I_mag'  ]-t_SDSS['dist_mod'])
y_z   = norm_var(t_SDSS['Z_mag'  ]-t_SDSS['dist_mod'])

y = np.transpose([
                    y_u  ,
                    y_g  ,
                    y_r  ,
                    y_i  ,
                    y_z
                        ]) [s_train]

y_verif = np.transpose([
                    y_u  ,
                    y_g  ,
                    y_r  ,
                    y_i  ,
                    y_z
                            ]) [~s_train]

gpr = GaussianProcessRegressor().fit(X, y)
score = gpr.score(X, y) # 0.84
t0 = time.time()
y_pred = gpr.predict(X_verif)#, return_std=True)
print(np.min(y_verif-y_pred), np.mean(y_verif-y_pred), np.median(y_verif-y_pred), np.max(y_verif-y_pred), np.std(y_verif-y_pred) )
print((time.time()-t0)/len(y_verif))
sys.exit()


# slice X_all in 100k
N_step = 10_000
Y_result = []
for j_slice in np.arange(0, len(x_0), N_step):
    t0 = time.time()
    Y_result.append( gpr.predict(X_all[j_slice:j_slice+N_step]) )
    print(N_step, (time.time()-t0)/N_step)


y_all = np.vstack(( Y_result ))

t_SDSS['u_GP'  ] = de_norm_var( y_all.T[0] , np.min(t_SDSS['U_mag']-t_SDSS['dist_mod']), np.max(t_SDSS['U_mag']-t_SDSS['dist_mod']) ) + t_SDSS['dist_mod']
t_SDSS['g_GP'  ] = de_norm_var( y_all.T[1] , np.min(t_SDSS['G_mag']-t_SDSS['dist_mod']), np.max(t_SDSS['G_mag']-t_SDSS['dist_mod']) ) + t_SDSS['dist_mod']
t_SDSS['r_GP'  ] = de_norm_var( y_all.T[2] , np.min(t_SDSS['R_mag']-t_SDSS['dist_mod']), np.max(t_SDSS['R_mag']-t_SDSS['dist_mod']) ) + t_SDSS['dist_mod']
t_SDSS['i_GP'  ] = de_norm_var( y_all.T[3] , np.min(t_SDSS['I_mag']-t_SDSS['dist_mod']), np.max(t_SDSS['I_mag']-t_SDSS['dist_mod']) ) + t_SDSS['dist_mod']
t_SDSS['z_GP'  ] = de_norm_var( y_all.T[4] , np.min(t_SDSS['Z_mag']-t_SDSS['dist_mod']), np.max(t_SDSS['Z_mag']-t_SDSS['dist_mod']) ) + t_SDSS['dist_mod']

p2_out_SDSS = path_2_SDSS[:-5]+'_SDSStrainedGPmags_0.fits'
t_SDSS.write(p2_out_SDSS, overwrite = True)
print(p2_out_SDSS)

p2_model = os.path.join(model_dir, "gal_GP_model_SDSStraining.pkl")
from pickle import dump
with open(p2_model, "wb") as f:
    dump(gpr, f, protocol=5)

print(p2_model, 'written')

def get_outlier_N(q1, q2, threshold):
    return len( (abs(t_SDSS[q1]-t_SDSS[q2])>threshold).nonzero()[0] )

def get_mean_diff(q1, q2):
    return np.mean(t_SDSS[q1]-t_SDSS[q2])
def get_std_diff(q1, q2):
    return np.std(t_SDSS[q1]-t_SDSS[q2])

print('band, mean(mag diff), std(mag diff), mean uncertainty (on input mag)')
print( 'U'  ,np.round(get_mean_diff('u_GP'  , 'U'  ), 4),np.round(get_std_diff('u_GP'  , 'U'  ), 4), np.round(t_SDSS['eU'  ].mean(), 4) )
print( 'G'  ,np.round(get_mean_diff('g_GP'  , 'G'  ), 4),np.round(get_std_diff('g_GP'  , 'G'  ), 4), np.round(t_SDSS['eG'  ].mean(), 4) )
print( 'R'  ,np.round(get_mean_diff('r_GP'  , 'R'  ), 4),np.round(get_std_diff('r_GP'  , 'R'  ), 4), np.round(t_SDSS['eR'  ].mean(), 4) )
print( 'I'  ,np.round(get_mean_diff('i_GP'  , 'I'  ), 4),np.round(get_std_diff('i_GP'  , 'I'  ), 4), np.round(t_SDSS['eI'  ].mean(), 4) )
print( 'Z'  ,np.round(get_mean_diff('z_GP'  , 'Z'  ), 4),np.round(get_std_diff('z_GP'  , 'Z'  ), 4), np.round(t_SDSS['eZ'  ].mean(), 4) )


N_tot = len(t_SDSS)
print('outlier fraction with Delta MAG > 1')
print( 'U'  , np.round( get_outlier_N('u_GP'  , 'U'  , 1)/N_tot * 100, 3) , '%')
print( 'G'  , np.round( get_outlier_N('g_GP'  , 'G'  , 1)/N_tot * 100, 3) , '%')
print( 'R'  , np.round( get_outlier_N('r_GP'  , 'R'  , 1)/N_tot * 100, 3) , '%')
print( 'I'  , np.round( get_outlier_N('i_GP'  , 'I'  , 1)/N_tot * 100, 3) , '%')
print( 'Z'  , np.round( get_outlier_N('z_GP'  , 'Z'  , 1)/N_tot * 100, 3) , '%')

print('outlier fraction with Delta MAG > 0.5')
print( 'U'  , np.round( get_outlier_N('u_GP'  , 'U'  , 0.5)/N_tot * 100, 3) , '%')
print( 'G'  , np.round( get_outlier_N('g_GP'  , 'G'  , 0.5)/N_tot * 100, 3) , '%')
print( 'R'  , np.round( get_outlier_N('r_GP'  , 'R'  , 0.5)/N_tot * 100, 3) , '%')
print( 'I'  , np.round( get_outlier_N('i_GP'  , 'I'  , 0.5)/N_tot * 100, 3) , '%')
print( 'Z'  , np.round( get_outlier_N('z_GP'  , 'Z'  , 0.5)/N_tot * 100, 3) , '%')

#sys.exit()

# Here you can replace pickle with joblib or cloudpickle
from pickle import load
with open(p2_model, "rb") as f:
    gpr = load(f)

#
# COSMOS
#
t_cosmos = Table.read(path_2_COSMOS)
good = (t_cosmos['photoz']>0.05 )&( t_cosmos['photoz']< 6. ) & ( t_cosmos['K']>0 )& ( t_cosmos['MK']<0 )&( t_cosmos['MK']>-40 )&( t_cosmos['LGM_TOT_P50']<13 )&( t_cosmos['LGM_TOT_P50']>6 )
t_cosmos = t_cosmos[good]
dm_values = dm_itp(t_cosmos['photoz'].data.data)
t_cosmos['dist_mod'] = dm_values

# input features
x_0 = norm_var_v(t_cosmos['mass_med'], np.min(t_SDSS['mass_med']), np.max(t_SDSS['mass_med']))
x_1 = norm_var_v(t_cosmos['ssfr_med'], np.min(t_SDSS['SPECSFR_TOT_P50'] ), np.max(t_SDSS['SPECSFR_TOT_P50'] ))
X_all_cosmos = np.transpose([x_0, x_1])

sys.exit()
# slice X_all_cosmos in many small instances

t0 = time.time()
y_all_cosmos = gpr.predict(X_all_cosmos)#, return_std=True)
print((time.time()-t0)/len(X_all_cosmos))

t_cosmos['u_GP'] = de_norm_var( y_all_cosmos.T[0] , np.min(t_SDSS['U']-t_SDSS['dist_mod']), np.max(t_SDSS['U']-t_SDSS['dist_mod']) ) + t_cosmos['dist_mod']
t_cosmos['g_GP'] = de_norm_var( y_all_cosmos.T[1] , np.min(t_SDSS['G']-t_SDSS['dist_mod']), np.max(t_SDSS['G']-t_SDSS['dist_mod']) ) + t_cosmos['dist_mod']
t_cosmos['r_GP'] = de_norm_var( y_all_cosmos.T[2] , np.min(t_SDSS['R']-t_SDSS['dist_mod']), np.max(t_SDSS['R']-t_SDSS['dist_mod']) ) + t_cosmos['dist_mod']
t_cosmos['i_GP'] = de_norm_var( y_all_cosmos.T[3] , np.min(t_SDSS['I']-t_SDSS['dist_mod']), np.max(t_SDSS['I']-t_SDSS['dist_mod']) ) + t_cosmos['dist_mod']
t_cosmos['z_GP'] = de_norm_var( y_all_cosmos.T[4] , np.min(t_SDSS['Z']-t_SDSS['dist_mod']), np.max(t_SDSS['Z']-t_SDSS['dist_mod']) ) + t_cosmos['dist_mod']
t_cosmos['y_GP'] = de_norm_var( y_all_cosmos.T[5] , np.min(t_SDSS['NUV']-t_SDSS['dist_mod']), np.max(t_SDSS['NUV']-t_SDSS['dist_mod']) ) + t_cosmos['dist_mod']
t_cosmos['j_GP'] = de_norm_var( y_all_cosmos.T[6] , np.min(t_SDSS['J']-t_SDSS['dist_mod']), np.max(t_SDSS['J']-t_SDSS['dist_mod']) ) + t_cosmos['dist_mod']
t_cosmos['h_GP'] = de_norm_var( y_all_cosmos.T[7] , np.min(t_SDSS['I1']-t_SDSS['dist_mod']), np.max(t_SDSS['I1']-t_SDSS['dist_mod']) ) + t_cosmos['dist_mod']
t_cosmos['k_GP'] = de_norm_var( y_all_cosmos.T[8] , np.min(t_SDSS['K']-t_SDSS['dist_mod']), np.max(t_SDSS['K']-t_SDSS['dist_mod']) ) + t_cosmos['dist_mod']

p2_out_SDSSMOS = path_2_COSMOS[:-5]+'_GPmags.fits'
t_cosmos.write(p2_out_SDSSMOS, overwrite = True)

#sys.exit()

#gmag_verif = de_norm_var(y_verif.T[0],np.min(t_SDSS['G']-t_SDSS['dist_mod']), np.max(t_SDSS['G']-t_SDSS['dist_mod']))
#gmag_pred  = de_norm_var(y_pred .T[0],np.min(t_SDSS['G']-t_SDSS['dist_mod']), np.max(t_SDSS['G']-t_SDSS['dist_mod']))

#rmag_verif = de_norm_var(y_verif.T[1],np.min(t_SDSS['R']-t_SDSS['dist_mod']), np.max(t_SDSS['R']-t_SDSS['dist_mod']))
#rmag_pred  = de_norm_var(y_pred .T[1],np.min(t_SDSS['R']-t_SDSS['dist_mod']), np.max(t_SDSS['R']-t_SDSS['dist_mod']))

#print(np.min(gmag_verif-gmag_pred), np.mean(gmag_verif-gmag_pred), np.median(gmag_verif-gmag_pred), np.max(gmag_verif-gmag_pred), np.std(gmag_verif-gmag_pred) )
#print((time.time()-t0)/len(y_verif))

#plt.figure(figsize=(6,6))
## data
#plt.hist( np.log10(abs(gmag_verif-gmag_pred)), bins=50)
## models
##plt.xlabel(r'$\log_{10}(M_s/M_\odot)$')
#plt.xlabel(r'$\log_{10}(|g-g_{prediction}|)$')
##plt.yscale('log')
#plt.legend(frameon=False, loc=3)
##plt.ylim((-26,-15))
##plt.xlim((7.5, 12.5))
#plt.tight_layout()
#plt.grid()
#plt.savefig(os.path.join( fig_dir, "GP_gmag_hist.png") )
#plt.clf()

#plt.figure(figsize=(6,6))
## data
#plt.hexbin(t_SDSS['mass_med'][~s_train], np.log10(abs(gmag_verif-gmag_pred)), gridsize=20)
## models
#plt.xlabel(r'$\log_{10}(M_s/M_\odot)$')
#plt.ylabel(r'$\log_{10}(|g-g_{prediction}|)$')
##plt.yscale('log')
#plt.legend(frameon=False, loc=3)
##plt.ylim((-26,-15))
##plt.xlim((7.5, 12.5))
#plt.tight_layout()
#plt.grid()
#plt.savefig(os.path.join( fig_dir, "GP_Mstar_gmag.png") )
#plt.clf()


#plt.figure(figsize=(6,6))
## data
#plt.hexbin(t_SDSS['SPECSFR_TOT_P50'][~s_train], np.log10(abs(gmag_verif-gmag_pred)), gridsize=20)
## models
#plt.xlabel(r'$\log_{10}(SFR/[M_\odot/yr])$')
#plt.ylabel(r'$\log_{10}(|g-g_{prediction}|)$')
##plt.yscale('log')
#plt.legend(frameon=False, loc=3)
##plt.ylim((-26,-15))
##plt.xlim((7.5, 12.5))
#plt.tight_layout()
#plt.grid()
#plt.savefig(os.path.join( fig_dir, "GP_SFR_gmag.png") )
#plt.clf()

#plt.figure(figsize=(6,6))
## data
#plt.hexbin(t_SDSS['dist_mod'][~s_train], np.log10(abs(gmag_verif-gmag_pred)), gridsize=20)
## models
#plt.xlabel(r'Distance Modulus')
#plt.ylabel(r'$\log_{10}(|g-g_{prediction}|)$')
##plt.yscale('log')
#plt.legend(frameon=False, loc=3)
##plt.ylim((-26,-15))
##plt.xlim((7.5, 12.5))
#plt.tight_layout()
#plt.grid()
#plt.savefig(os.path.join( fig_dir, "GP_distmod_gmag.png") )
#plt.clf()


#plt.figure(figsize=(6,6))
## data
#plt.hist( np.log10(abs(rmag_verif-rmag_pred)), bins=50)
## models
##plt.xlabel(r'$\log_{10}(M_s/M_\odot)$')
#plt.xlabel(r'$\log_{10}(|r-r_{prediction}|)$')
##plt.yscale('log')
#plt.legend(frameon=False, loc=3)
##plt.ylim((-26,-15))
##plt.xlim((7.5, 12.5))
#plt.tight_layout()
#plt.grid()
#plt.savefig(os.path.join( fig_dir, "GP_rmag_hist.png") )
#plt.clf()

#plt.figure(figsize=(6,6))
## data
#plt.hexbin(t_SDSS['mass_med'][~s_train], np.log10(abs(rmag_verif-rmag_pred)), gridsize=20)
## models
#plt.xlabel(r'$\log_{10}(M_s/M_\odot)$')
#plt.ylabel(r'$\log_{10}(|r-r_{prediction}|)$')
##plt.yscale('log')
#plt.legend(frameon=False, loc=3)
##plt.ylim((-26,-15))
##plt.xlim((7.5, 12.5))
#plt.tight_layout()
#plt.grid()
#plt.savefig(os.path.join( fig_dir, "GP_Mstar_rmag.png") )
#plt.clf()


#plt.figure(figsize=(6,6))
## data
#plt.hexbin(t_SDSS['SPECSFR_TOT_P50'][~s_train], np.log10(abs(rmag_verif-rmag_pred)), gridsize=20)
## models
#plt.xlabel(r'$\log_{10}(SFR/[M_\odot/yr])$')
#plt.ylabel(r'$\log_{10}(|r-r_{prediction}|)$')
#plt.legend(frameon=False, loc=3)
##plt.ylim((-26,-15))
##plt.xlim((7.5, 12.5))
#plt.tight_layout()
#plt.grid()
#plt.savefig(os.path.join( fig_dir, "GP_SFR_rmag.png") )
#plt.clf()

#plt.figure(figsize=(6,6))
## data
#plt.hexbin(t_SDSS['dist_mod'][~s_train], np.log10(abs(rmag_verif-rmag_pred)), gridsize=20)
## models
#plt.xlabel(r'Distance Modulus')
#plt.ylabel(r'$\log_{10}(|r-r_{prediction}|)$')
##plt.yscale('log')
#plt.legend(frameon=False, loc=3)
##plt.ylim((-26,-15))
##plt.xlim((7.5, 12.5))
#plt.tight_layout()
#plt.grid()
#plt.savefig(os.path.join( fig_dir, "GP_distmod_rmag.png") )
#plt.clf()

