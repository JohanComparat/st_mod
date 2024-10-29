"""
What it does
------------

Setup of the GAL pipeline

Measurement of the mean relation between stellar mass, K band absolute magnitude and redshift

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

fig_dir = os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/model_GAL/figures')
os.system('mkdir -p '+fig_dir)
model_dir = os.path.join(os.environ['GIT_STMOD_DATA'], 'data', 'models','model_GAL')
os.system('mkdir -p '+model_dir)
print('figures are written here : ', fig_dir)
print('model files are written here : ', model_dir)
path_2_kids = os.path.join(os.environ['DATA'], 'GAMA', 'forJohan.fits')
path_2_COSMOS = os.path.join(os.environ['DATA'], 'COSMOS', 'photoz_vers2.0_010312.fits')
#
# COSMOS
#
t_cosmos = Table.read(path_2_COSMOS)
good = (t_cosmos['photoz']>0.05 )&( t_cosmos['photoz']< 6. ) & ( t_cosmos['K']>0 )& ( t_cosmos['MK']<0 )&( t_cosmos['MK']>-40 )&( t_cosmos['mass_med']<14 )&( t_cosmos['mass_med']>6 )
t_cosmos = t_cosmos[good]
#
# GAMA
#
t_kids = Table.read(path_2_kids)
keep_kids = (8.9 - 2.5*np.log10(t_kids['flux_Kt'])>0) & (t_kids['Z']>0.01) & (t_kids['Z']<1.9)
t_kids = Table(t_kids[keep_kids])
kmag = 8.9 - 2.5*np.log10(t_kids['flux_Kt'])
zmag = 8.9 - 2.5*np.log10(t_kids['flux_Zt'])
rmag = 8.9 - 2.5*np.log10(t_kids['flux_rt'])
imag = 8.9 - 2.5*np.log10(t_kids['flux_it'])
# K band absolute magnitude :
KCORR_DATA = np.loadtxt( os.path.join( os.environ['GIT_STMOD_DATA'], 'data', 'models','model_GAL', 'VISTA_Ks_kcorrections.txt'), unpack = True)
kcorr_itp = interp1d(KCORR_DATA[0], KCORR_DATA[3])
dm_values = dm_itp(t_kids['Z'].data.data)
kmag_abs = kmag - ( dm_values + kcorr_itp(t_kids['Z'].data.data) )
t_kids['Kmag_abs'] = np.array(kmag_abs)
t_kids['dist_mod'] = dm_values
t_kids['ms'] = np.array(np.log10( t_kids['StellarMass_50'] ))

dz = 0.1
zall = np.arange(0., 2.+dz, dz)
params, covar = [], []
params2, covar2 = [], []
params_lm = []
resid, resid2 = [], []
for zmin, zmax in zip(zall[:-1], zall[1:]):
    # full data set
    s00el =  (t_cosmos['photoz']>zmin) & (t_cosmos['photoz']<=zmax)
    s00el2 = (t_kids['Z']>zmin) & (t_kids['Z']<=zmax)
    # data for fitting the relation
    sel =  (t_cosmos['photoz']>zmin) & (t_cosmos['photoz']<=zmax) &( t_cosmos['mass_med']>8+zmin )
    sel2 = (t_kids['Z']>zmin) & (t_kids['Z']<=zmax) &( t_kids['ms']>8+zmin )
    # function to fit, both slope and amplitude are free
    fun2 = lambda x, a, b : a*x+b
    # fixing the slope to avoid jumps with redshift correlated to the second parameters
    fun = lambda x, b : -2.15 * x + b
    # data to fit : by joining the 2 data sets
    x_data = np.hstack(( t_cosmos['mass_med'][sel], t_kids['ms'][sel2]))
    y_data = np.hstack(( t_cosmos['MK'][sel]      , t_kids['Kmag_abs'][sel2]))
    # fit
    popt, pcov= curve_fit(fun, x_data, y_data, p0=(0.))#, sigma=y_err)
    popt2, pcov2= curve_fit(fun2, x_data, y_data, p0=(-2.15, 0.))#, sigma=y_err)
    params.append(popt)
    covar.append(pcov)
    params2.append(popt2)
    covar2.append(pcov2)
    # alternative with linmix, but slow
    #np.random.seed(2)
    #lm = linmix.LinMix(x_data, y_data, K=2)
    #lm.run_mcmc(miniter=200, maxiter=600)
    #params_lm.append([np.median(lm.chain['beta']), np.median(lm.chain['alpha'])])
    #print(np.round(zmin,1),np.round(zmax,1),'curve_fit',np.round(popt,2), np.round(pcov[0][0],4), np.round(pcov[1][1],4))#,'linmix',params_lm[-1])# = 0.05
    print(np.round(zmin,1),np.round(zmax,1),'curve_fit',np.round(popt,2), np.round(pcov[0][0]**0.5,4))
    print(np.round(zmin,1),np.round(zmax,1),'curve_fit',np.round(popt2[1],2),np.round(popt2[0],2), np.round(pcov2[0][0]**0.5,4), np.round(pcov2[1][1]**0.5,4))
    #
    # figure
    #
    plt.figure(figsize=(6,6))
    # data
    plt.plot(t_cosmos['mass_med'][s00el], t_cosmos['MK'][s00el], 'k,')
    plt.plot(t_kids['ms'][s00el2], t_kids['Kmag_abs'][s00el2], 'g,')
    plt.plot(t_cosmos['mass_med'][sel], t_cosmos['MK'][sel], 'k+', alpha=0.5, label='COSMOS, N='+str(len(t_cosmos['mass_med'][sel])))
    plt.plot(t_kids['ms'][sel2], t_kids['Kmag_abs'][sel2], 'gx', alpha=0.5, label='GAMA, N='+str(len(t_kids['ms'][sel2])))
    # models
    xs2 = np.arange(7.5, 12.5, 0.01)
    plt.plot(xs2, fun(xs2, popt[0]), 'b--', lw=3, label=r'$MK=-2.15\log_{10}(M_s)$'+'+'+str(np.round(popt[0],2)))
    plt.plot(xs2, fun2(xs2, popt2[0], popt2[1]), 'r--', lw=2, label=r'$MK=$'+str(np.round(popt2[0],2))+r'$\log_{10}(M_s)$'+'+'+str(np.round(popt2[1],2)))
    plt.xlabel(r'$\log_{10}(M_s/M_\odot)$')
    plt.ylabel(r'MK')
    plt.legend(frameon=False, loc=3)
    plt.ylim((-26,-15))
    plt.xlim((7.5, 12.5))
    plt.grid()
    plt.title(str(np.round(zmin,2))+r'$<z<$'+str(np.round(zmax,2)) )#+ r", $F_X>1\times10^{-17}$ ")
    plt.savefig(os.path.join( fig_dir, "fit_mass_Kmag_"+str(np.round(zmin,2))+".png") )
    plt.clf()
    #
    # RESIDUALS
    #
    # data
    dy_B = y_data - fun(x_data, popt[0])
    dy_C = y_data - fun2(x_data, popt2[0], popt2[1])
    #fitting function
    funN = lambda x, loc, scale : norm.cdf(x, loc=loc, scale=scale)
    # figure
    bins=np.arange(-1.5,1.5,0.05)
    xbins=bins[1:]-0.025
    plt.figure(2, (6,6))
    # model & data 1 params
    n_t2 = plt.hist(dy_B , bins=bins, density=True, cumulative=True, histtype='step', label='residuals 1p' )[0]
    pt, ct = curve_fit(funN, xbins, n_t2, p0=(0,0.4))
    plt.plot(xbins, norm.cdf(xbins, loc=pt[0], scale=pt[1]), label='N1p('+str(np.round(pt[0],2))+', '+str(np.round(pt[1],2))+')', ls='dashed')
    # model & data 2 params
    n_t2 = plt.hist(dy_C , bins=bins, density=True, cumulative=True, histtype='step', label='residuals 2p' )[0]
    pt2, ct2 = curve_fit(funN, xbins, n_t2, p0=(0,0.4))
    plt.plot(xbins, norm.cdf(xbins, loc=pt2[0], scale=pt2[1]), label='N2p('+str(np.round(pt2[0],2))+', '+str(np.round(pt2[1],2))+')', ls='dashed')
    plt.ylabel('counts')
    plt.xlabel(r'delta mag')
    plt.grid()
    plt.legend(frameon=False, loc=0)
    plt.title(str(np.round(zmin,2))+r'$<z<$'+str(np.round(zmax,2)))
    plt.savefig(os.path.join( fig_dir, "residual_fit_mass_Kmag_"+str(np.round(zmin,2))+".png"))
    plt.clf()
    print('residuals', pt, pt2)
    resid.append( pt)# ,ct] )
    resid2.append( pt2)# ,ct2] )

p_a = np.transpose(params)
p_b = np.transpose(resid)
OUT = np.transpose([ zall[:-1], zall[1:],  p_a[0], p_b[0], p_b[1] ])
np.savetxt(os.path.join( model_dir, 'logMs-MK-look-up-table_2023Aug22.txt'), OUT)

p_a = np.transpose(params2)
p_b = np.transpose(resid2)
OUT = np.transpose([zall[:-1], zall[1:],  p_a[0], p_a[1], p_b[0], p_b[1] ])
np.savetxt(os.path.join( model_dir, 'logMs-MK-look-up-table-2parameters_2023Aug22.txt'), OUT)
