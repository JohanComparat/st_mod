
from astropy.table import Table, Column, hstack, vstack
import time
t0 = time.time()
from xsingle_cubes_lib import *
import matplotlib.ticker as ticker
import time
t0 = time.time()

import numpy as n
import numpy as np
import os, sys, glob
import astropy.io.fits as fits
from astropy.table import Table, Column
import csv
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as p
from matplotlib.colors import LogNorm
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gs
import astropy.units as u
import astropy.constants as cc
from scipy.optimize import curve_fit
def loadprofile(path):
    txt_pro=n.loadtxt(path)
    return txt_pro
#from dark_emulator import model_hod
from astropy.cosmology import FlatLambdaCDM
cosmoUNIT = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
h = 0.6774
L_box = 1000.0 / h
cosmo = cosmoUNIT
from scipy.interpolate import interp1d



# PATH TO THE CUBES
root_directory = os.path.join(os.environ['HOME'], 'coor')
output_directory = os.path.join(os.environ['YIERASS'])
bins = n.hstack(( 0., 10., 30., 50., 75., 100.,150, 200.,250,300.,350,400.,500.,600.,700.,800.,900., 1000. ,1100,1200,1300,1400,1600,1800,2000,2200,2400,2600,2800,3000))#l_bins # n.hstack(( 0., 10**n.arange(1.2, n.log10(3000.), 0.2), 3000. ))
area_shell= ( bins[1:]**2 - bins[:-1]**2 ) * n.pi
MList=n.log10(n.array([5.5e9,1.8e10,5.5e10,1.3e11,2.2e11,4.1e11]))
zlist=n.array([0.07,0.1,0.14,0.28,0.35,0.35])

Rvirlist=n.array([115.3,150.6, 215.3, 306.3, 423.8,672.7])
Rvirlist_blue=n.array([115.4,147.2, 197.9, 260.1, 337.6,563.3])
Rvirlist_red=n.array([115.0, 157.5, 231.1, 323.1, 437.7,678.5])

colors_2D_gal = ['pink','purple',p.cm.tab20b(0)]
colors_2D_sf = [p.cm.tab20b(2),p.cm.tab20c(0),p.cm.tab20b(0)]
colors_2D_qu = ['pink','orange','darkred']
def plotlx_m(title,r,mask,sel,GAL):
    data_file = os.path.join(root_directory, 'DR9_Ebins_s4', 'data', 'lx_m','LX_M_'+r+'_' + mask + '_' + sel + '_' + GAL + '.txt')
    lx_m=n.transpose(n.loadtxt(data_file))
    #([mmm, low,high,lllx, lllx_err, xrb[:, 0], xrb[:, 1], xrb[:, 2], agn[:, 0], agn[:, 1], agn_0p8[:, 0], agn_0p8[:, 1]])
    figure_name = os.path.join( root_directory, 'DR9_Ebins_s4','fig_'+mask+'_new','LX_M_'+r+'_' + mask + '_' + sel + '_' + GAL + '.pdf' )
    title_list=os.path.basename(data_file).split('_')
    fig, ax = p.subplots(figsize=[9,9])
    p.title(title, fontsize=22)
    p.tick_params(which='both',length=6, width=1.2)
    ax.errorbar(lx_m[0], lx_m[3],xerr=[lx_m[1],lx_m[2]],yerr=lx_m[4], lw=3, ls='',color='purple',zorder=20)
    ax.fill_between(lx_m[0],lx_m[3]-lx_m[4],lx_m[3]+lx_m[4],color='purple',alpha=0.5,label='Isolated galaxy')
    ax.axhline(lx_m[-1][0]*0.01,lw=3, color='purple', ls='dashed', label='1% background luminosity')
    if (r=='0_1') | (r=='0_0p5'):
        ax.errorbar(lx_m[0],lx_m[5],xerr=[lx_m[1],lx_m[2]],yerr=[lx_m[6],lx_m[7]],color='pink',ls='',lw=3,label='XRB')
        ax.errorbar(lx_m[0],lx_m[8],xerr=[lx_m[1],lx_m[2]],yerr=lx_m[9],color='gray', ls='', lw=3,label='AGN (50%)')
        ax.errorbar(lx_m[0],lx_m[10],xerr=[lx_m[1],lx_m[2]],yerr=lx_m[11],color='brown', ls='', lw=3,label='AGN (80%)')
    ax.set_xlabel(r'$M_*$ [M$_\odot$]', fontsize=22)
    ax.set_ylabel(r'$L_{X,0.5-2keV}$ [erg/s]', fontsize=22)  # ; keV^{-1}]$')# / (\Delta keV=40/1000)$ ')
    ax.legend(frameon=True, markerscale=1, fontsize=20, loc='upper left', ncol=1)  # ,title=legend_title)
    # p.xlim((E_min-0.01, E_max+0.01))
    ax.set_yscale('log')
    ax.set_xscale('log')
    p.xticks(fontsize=22)
    p.yticks(fontsize=22)
    ax.set_xlim((3e9, 1e12))
    ax.set_ylim((5e38, 5e43))
    fig.tight_layout()  # (rect=[0, 0.03, 1, 0.95])
    fig.savefig(figure_name)
    fig.clf()
plotlx_m(r"$0<R_p<R_{vir}$",'0_1','all','2D','GAL')
plotlx_m(r"$1/2R_{vir}<R_p<R_{vir}$",'0p5_1','all','2D','GAL')
plotlx_m(r"$0<R_p<1/2R_{vir}$",'0_0p5','all','2D','GAL')
plotlx_m(r"$0<R_p<R_{vir}$",'0_1','all','2D','SF')
plotlx_m(r"$1/2R_{vir}<R_p<R_{vir}$",'0p5_1','all','2D','SF')
plotlx_m(r"$0<R_p<1/2R_{vir}$",'0_0p5','all','2D','SF')
plotlx_m(r"$0<R_p<R_{vir}$",'0_1','all','2D','QU')
plotlx_m(r"$1/2R_{vir}<R_p<R_{vir}$",'0p5_1','all','2D','QU')
plotlx_m(r"$0<R_p<1/2R_{vir}$",'0_0p5','all','2D','QU')

#######GAL profile######

def plot_profile(mask,sel,GAL,colors_2D_gal):
    profile_path = glob.glob(os.path.join(root_directory, 'DR9_Ebins_s4','data','profile','profile_*_'+mask+'_'+sel+'_'+GAL+'.txt'))
    sort=n.argsort(n.array([os.path.basename(i).split('_')[1] for i in profile_path]).astype(float))
    profile_path=n.array(profile_path)[sort]
    figure_name = os.path.join( root_directory, 'DR9_Ebins_s4','fig_'+mask+'_new','profile_'+mask+'_'+sel+'_'+GAL+'_withBG.pdf' )
    fig, ax = p.subplots(figsize=[9,9])
    p.title("Isolated galaxy", fontsize=22)
    p.tick_params(which='both',length=6, width=1.2)
    i=0
    for p_2_pro in profile_path[1:4]:
        Profile=n.transpose(loadprofile(p_2_pro))
        xb = Profile[0]
        totalerr=Profile[2]
        BG_new_mean=Profile[3][0]
        BG_new_std=Profile[4][0]
        title_list=os.path.basename(p_2_pro).split('_')
        #[xb,Ix_mean,totalerr, n.ones_like(xb)*BG_new_mean, n.ones_like(xb)*BG_new_std,ps_profile_normed])
        ax.errorbar(xb, Profile[1], yerr=totalerr, xerr=[xb - bins[:-1], bins[1:] - xb], lw=3, ls='', color=colors_2D_gal[i],zorder=20, label=title_list[1]+r'<log($M_*$)<'+title_list[3])
        ax.fill_between(xb, Profile[1] -totalerr, Profile[1] + totalerr, color=colors_2D_gal[i], alpha=0.3, zorder=5)
        ax.axhline(BG_new_mean, lw=3, color=colors_2D_gal[i], ls='dashed')
        ax.fill_between(xb, BG_new_mean -BG_new_std, BG_new_mean + BG_new_std, color=colors_2D_gal[i], alpha=0.3, zorder=3)
        ax.step(bins[1:], Profile[5], lw=2, zorder=3, ls=':', color=colors_2D_gal[i])
        ax.axvline(Profile[-1][0], color=colors_2D_gal[i], ls='--')
        i+=1
    ax.axhline(1, lw=3, color=colors_2D_gal[1], ls='dashed', label='background')
    ax.axhline(1, lw=2, zorder=3, ls=':', color=colors_2D_gal[1], label='empirical PSF')

    ax.set_xlabel(r'$R_p$ [kpc]',fontsize=22)
    ax.set_ylabel(r'$S_X$ $[erg\; s^{-1}\; kpc^{-2}]$',fontsize=22) #; keV^{-1}]$')# / (\Delta keV=40/1000)$ ')
    ax.legend(frameon=True, markerscale=25 ,fontsize=21, loc=1, ncol=1)#,title=legend_title)
    p.xticks(fontsize=22)
    p.yticks(fontsize=22)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim((4, 3000))
    ax.set_ylim((7e36, 5e37))
    fig.tight_layout()#(rect=[0, 0.03, 1, 0.99])
    fig.savefig( figure_name )
    fig.clf()
    figure_name = os.path.join(root_directory, 'DR9_Ebins_s4', 'fig_' + mask + '_new',
                               'profile_' + mask + '_' + sel + '_' + GAL + '.pdf')
    fig, ax = p.subplots(figsize=[9, 9])
    p.title("Isolated galaxy, subtract BG", fontsize=22)
    p.tick_params(which='both', length=6, width=1.2)
    i = 0
    for p_2_pro in profile_path[1:4]:
        Profile = n.transpose(loadprofile(p_2_pro))
        xb = Profile[0]
        totalerr = Profile[2]
        BG_new_mean = Profile[3][0]
        BG_new_std = Profile[4][0]
        title_list = os.path.basename(p_2_pro).split('_')
        # [xb,Ix_mean,totalerr, n.ones_like(xb)*BG_new_mean, n.ones_like(xb)*BG_new_std,ps_profile_normed])
        sel_lim = (Profile[1] - BG_new_mean - totalerr > 0)
        ax.errorbar(xb[sel_lim], Profile[1][sel_lim] - BG_new_mean, yerr=totalerr[sel_lim],
                    xerr=[xb[sel_lim] - bins[:-1][sel_lim], bins[1:][sel_lim] - xb[sel_lim]], lw=2, ls='',
                    color=colors_2D_gal[i], label=title_list[1] + r'<log($M_*$)<' + title_list[3],
                    zorder=5)
        for x, y, dy, wid, xup, xlow, y_ln, x_low in zip(xb, Profile[ 1] - BG_new_mean + totalerr, totalerr / 4,
                                                         bins[1:] - bins[:-1],
                                                         xb - bins[:-1], bins[1:] - xb, Profile[ 1] - BG_new_mean,
                                                         Profile[1] - BG_new_mean - totalerr):
            if x_low <= 0:
                ax.errorbar(x, y, xerr=[[xup], [xlow]], yerr=dy, uplims=True, color=colors_2D_gal[i],
                            )
        ax.step(bins[1:],  Profile[5] - BG_new_mean, lw=2, zorder=3, ls=':', color=colors_2D_gal[i])
        ax.axhline(BG_new_mean / 100, lw=3, color=colors_2D_gal[i], alpha=0.7, ls='dashed')
        ax.axvline( Profile[-1][0], color=colors_2D_gal[i], ls='--')
        i += 1
    ax.axhline(1, lw=3, color=colors_2D_gal[1], ls='dashed', alpha=0.7, label='1% background')
    ax.axhline(1, lw=2, zorder=3, ls=':', color=colors_2D_gal[1], label='empirical PSF - background')
    ax.set_xlabel(r'$R_p$ [kpc]', fontsize=22)
    ax.set_ylabel(r'$S_X$ $[erg\; s^{-1}\; kpc^{-2}]$', fontsize=22)  # ; keV^{-1}]$')# / (\Delta keV=40/1000)$ ')
    ax.legend(frameon=True, markerscale=25, fontsize=20, loc=1, ncol=1)  # ,title=legend_title)
    # p.xlim((E_min-0.01, E_max+0.01))
    ax.set_yscale('log')
    ax.set_xscale('log')
    p.xticks(fontsize=22)
    p.yticks(fontsize=22)
    ax.set_xlim((4, 3000))
    ax.set_ylim((2e34, 5e37))
    fig.tight_layout()  # (rect=[0, 0.03, 1, 0.95])
    fig.savefig(figure_name)
    fig.clf()

plot_profile('all','2D','GAL',colors_2D_gal)
plot_profile('all','2D','SF',colors_2D_sf)
plot_profile('all','2D','QU',colors_2D_qu)


profile_path = glob.glob(os.path.join(root_directory, 'DR9_Ebins_s4','data','profile','profile_*_all_2D_GAL.txt'))
sort=n.argsort(n.array([os.path.basename(i).split('_')[1] for i in profile_path]).astype(float))
profile_path_GAL=n.array(profile_path)[sort]
profile_path = glob.glob(os.path.join(root_directory, 'DR9_Ebins_s4','data','profile','profile_*_all_2D_SF.txt'))
sort=n.argsort(n.array([os.path.basename(i).split('_')[1] for i in profile_path]).astype(float))
profile_path_SF=n.array(profile_path)[sort]
profile_path = glob.glob(os.path.join(root_directory, 'DR9_Ebins_s4','data','profile','profile_*_all_2D_QU.txt'))
sort=n.argsort(n.array([os.path.basename(i).split('_')[1] for i in profile_path]).astype(float))
profile_path_QU=n.array(profile_path)[sort]
profile_path = glob.glob(os.path.join(root_directory, 'DR9_Ebins_s4','data','profile','profile_*_all_full_GAL.txt'))
sort=n.argsort(n.array([os.path.basename(i).split('_')[1] for i in profile_path]).astype(float))
profile_path_GAL_full=n.array(profile_path)[sort]
profile_path = glob.glob(os.path.join(root_directory, 'DR9_Ebins_s4','data','profile','profile_*_all_full_SF.txt'))
sort=n.argsort(n.array([os.path.basename(i).split('_')[1] for i in profile_path]).astype(float))
profile_path_SF_full=n.array(profile_path)[sort]
profile_path = glob.glob(os.path.join(root_directory, 'DR9_Ebins_s4','data','profile','profile_*_all_full_QU.txt'))
sort=n.argsort(n.array([os.path.basename(i).split('_')[1] for i in profile_path]).astype(float))
profile_path_QU_full=n.array(profile_path)[sort]

def plot_profile_line(p_2_pro,colors_2D_gal,i):
    Profile = n.transpose(loadprofile(p_2_pro))
    xb = Profile[0]
    totalerr = Profile[2]
    BG_new_mean = Profile[3][0]
    BG_new_std = Profile[4][0]
    title_list = os.path.basename(p_2_pro).split('_')
    # [xb,Ix_mean,totalerr, n.ones_like(xb)*BG_new_mean, n.ones_like(xb)*BG_new_std,ps_profile_normed])
    ax.errorbar(xb, Profile[1], yerr=totalerr, xerr=[xb - bins[:-1], bins[1:] - xb], lw=3, ls='',
                color=colors_2D_gal[i], zorder=20, label=title_list[1] + r'<log($M_*$)<' + title_list[3]+', '+title_list[-2])
    ax.fill_between(xb, Profile[1] - totalerr, Profile[1] + totalerr, color=colors_2D_gal[i], alpha=0.3, zorder=5)
    ax.axhline(BG_new_mean, lw=3, color=colors_2D_gal[i], ls='dashed')
    ax.fill_between(xb, BG_new_mean - BG_new_std, BG_new_mean + BG_new_std, color=colors_2D_gal[i], alpha=0.3, zorder=3)
    ax.step(bins[1:], Profile[5], lw=2, zorder=3, ls=':', color=colors_2D_gal[i])
    ax.axvline(Profile[-1][0], color=colors_2D_gal[i], ls='--')

def plot_profile_line_sub(p_2_pro,colors_2D_gal,i):
    Profile = n.transpose(loadprofile(p_2_pro))
    xb = Profile[0]
    totalerr = Profile[2]
    BG_new_mean = Profile[3][0]
    BG_new_std = Profile[4][0]
    title_list = os.path.basename(p_2_pro).split('_')
    # [xb,Ix_mean,totalerr, n.ones_like(xb)*BG_new_mean, n.ones_like(xb)*BG_new_std,ps_profile_normed])
    sel_lim = (Profile[1] - BG_new_mean - totalerr > 0)
    ax.errorbar(xb[sel_lim], Profile[1][sel_lim] - BG_new_mean, yerr=totalerr[sel_lim],
                xerr=[xb[sel_lim] - bins[:-1][sel_lim], bins[1:][sel_lim] - xb[sel_lim]], lw=2, ls='',
                color=colors_2D_gal[i], label=title_list[1] + r'<log($M_*$)<' + title_list[3]+', '+title_list[-2],
                zorder=5)
    for x, y, dy, wid, xup, xlow, y_ln, x_low in zip(xb, Profile[ 1] - BG_new_mean + totalerr, totalerr / 4,
                                                     bins[1:] - bins[:-1],
                                                     xb - bins[:-1], bins[1:] - xb, Profile[ 1] - BG_new_mean,
                                                     Profile[1] - BG_new_mean - totalerr):
        if x_low <= 0:
            ax.errorbar(x, y, xerr=[[xup], [xlow]], yerr=dy, uplims=True, color=colors_2D_gal[i],
                        )
    ax.step(bins[1:],  Profile[5] - BG_new_mean, lw=2, zorder=3, ls=':', color=colors_2D_gal[i])
    ax.axhline(BG_new_mean / 100, lw=3, color=colors_2D_gal[i], alpha=0.7, ls='dashed')
    ax.axvline( Profile[-1][0], color=colors_2D_gal[i], ls='--')


figure_name = os.path.join( root_directory, 'DR9_Ebins_s4','fig_all_new','profile_all_M31_withBG.pdf' )
fig, ax = p.subplots(figsize=[9,9])
p.title("M31", fontsize=22)
p.tick_params(which='both',length=6, width=1.2)

plot_profile_line(profile_path_GAL[3],colors_2D_gal,0)
plot_profile_line(profile_path_GAL_full[3],colors_2D_gal,1)

ax.axhline(1, lw=3, color=colors_2D_gal[1], ls='dashed', label='background')
ax.axhline(1, lw=2, zorder=3, ls=':', color=colors_2D_gal[1], label='empirical PSF')

ax.set_xlabel(r'$R_p$ [kpc]',fontsize=22)
ax.set_ylabel(r'$S_X$ $[erg\; s^{-1}\; kpc^{-2}]$',fontsize=22) #; keV^{-1}]$')# / (\Delta keV=40/1000)$ ')
ax.legend(frameon=True, markerscale=25 ,fontsize=21, loc=1, ncol=1)#,title=legend_title)
p.xticks(fontsize=22)
p.yticks(fontsize=22)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlim((4, 3000))
ax.set_ylim((7e36, 5e37))
fig.tight_layout()#(rect=[0, 0.03, 1, 0.99])
fig.savefig( figure_name )
fig.clf()

figure_name = os.path.join(root_directory, 'DR9_Ebins_s4', 'fig_all_new','profile_all_M31.pdf')
fig, ax = p.subplots(figsize=[9, 9])
p.title("M31, subtract BG", fontsize=22)
p.tick_params(which='both', length=6, width=1.2)

plot_profile_line_sub(profile_path_GAL[3],colors_2D_gal,0)
plot_profile_line_sub(profile_path_GAL_full[3],colors_2D_gal,1)

ax.axhline(1, lw=3, color=colors_2D_gal[1], ls='dashed', alpha=0.7, label='1% background')
ax.axhline(1, lw=2, zorder=3, ls=':', color=colors_2D_gal[1], label='empirical PSF - background')
ax.set_xlabel(r'$R_p$ [kpc]', fontsize=22)
ax.set_ylabel(r'$S_X$ $[erg\; s^{-1}\; kpc^{-2}]$', fontsize=22)  # ; keV^{-1}]$')# / (\Delta keV=40/1000)$ ')
ax.legend(frameon=True, markerscale=25, fontsize=20, loc=1, ncol=1)  # ,title=legend_title)
# p.xlim((E_min-0.01, E_max+0.01))
ax.set_yscale('log')
ax.set_xscale('log')
p.xticks(fontsize=22)
p.yticks(fontsize=22)
ax.set_xlim((4, 3000))
ax.set_ylim((2e34, 5e37))
fig.tight_layout()  # (rect=[0, 0.03, 1, 0.95])
fig.savefig(figure_name)
fig.clf()

profile_path = glob.glob(os.path.join(root_directory, 'DR9_Ebins_s4','data','profile','profile_*_nomask_2D_GAL.txt'))
sort=n.argsort(n.array([os.path.basename(i).split('_')[1] for i in profile_path]).astype(float))
profile_path_GAL_nomask=n.array(profile_path)[sort]
profile_path = glob.glob(os.path.join(root_directory, 'DR9_Ebins_s4','data','profile','profile_*_nomask_2D_SF.txt'))
sort=n.argsort(n.array([os.path.basename(i).split('_')[1] for i in profile_path]).astype(float))
profile_path_SF_nomask=n.array(profile_path)[sort]
profile_path = glob.glob(os.path.join(root_directory, 'DR9_Ebins_s4','data','profile','profile_*_nomask_2D_QU.txt'))
sort=n.argsort(n.array([os.path.basename(i).split('_')[1] for i in profile_path]).astype(float))
profile_path_QU_nomask=n.array(profile_path)[sort]
profile_path = glob.glob(os.path.join(root_directory, 'DR9_Ebins_s4','data','profile','profile_*_nomask_full_GAL.txt'))
sort=n.argsort(n.array([os.path.basename(i).split('_')[1] for i in profile_path]).astype(float))
profile_path_GAL_full_nomask=n.array(profile_path)[sort]
profile_path = glob.glob(os.path.join(root_directory, 'DR9_Ebins_s4','data','profile','profile_*_nomask_full_SF.txt'))
sort=n.argsort(n.array([os.path.basename(i).split('_')[1] for i in profile_path]).astype(float))
profile_path_SF_full_nomask=n.array(profile_path)[sort]
profile_path = glob.glob(os.path.join(root_directory, 'DR9_Ebins_s4','data','profile','profile_*_nomask_full_QU.txt'))
sort=n.argsort(n.array([os.path.basename(i).split('_')[1] for i in profile_path]).astype(float))
profile_path_QU_full_nomask=n.array(profile_path)[sort]

figure_name = os.path.join( root_directory, 'DR9_Ebins_s4','fig_all_new','profile_nomask_MW_withBG.pdf' )
fig, ax = p.subplots(figsize=[9,9])
p.title("MW,no mask", fontsize=22)
p.tick_params(which='both',length=6, width=1.2)

plot_profile_line(profile_path_GAL_nomask[2],colors_2D_gal,0)
plot_profile_line(profile_path_GAL_full_nomask[2],colors_2D_gal,1)

ax.axhline(1, lw=3, color=colors_2D_gal[1], ls='dashed', label='background')
ax.axhline(1, lw=2, zorder=3, ls=':', color=colors_2D_gal[1], label='empirical PSF')

ax.set_xlabel(r'$R_p$ [kpc]',fontsize=22)
ax.set_ylabel(r'$S_X$ $[erg\; s^{-1}\; kpc^{-2}]$',fontsize=22) #; keV^{-1}]$')# / (\Delta keV=40/1000)$ ')
ax.legend(frameon=True, markerscale=25 ,fontsize=21, loc=1, ncol=1)#,title=legend_title)
p.xticks(fontsize=22)
p.yticks(fontsize=22)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlim((4, 3000))
ax.set_ylim((7e36, 5e37))
fig.tight_layout()#(rect=[0, 0.03, 1, 0.99])
fig.savefig( figure_name )
fig.clf()

figure_name = os.path.join(root_directory, 'DR9_Ebins_s4', 'fig_all_new','profile_nomask_MW.pdf')
fig, ax = p.subplots(figsize=[9, 9])
p.title("MW, no mask, subtract BG", fontsize=22)
p.tick_params(which='both', length=6, width=1.2)

plot_profile_line_sub(profile_path_GAL_nomask[2],colors_2D_gal,0)
plot_profile_line_sub(profile_path_GAL_full_nomask[2],colors_2D_gal,1)

ax.axhline(1, lw=3, color=colors_2D_gal[1], ls='dashed', alpha=0.7, label='1% background')
ax.axhline(1, lw=2, zorder=3, ls=':', color=colors_2D_gal[1], label='empirical PSF - background')
ax.set_xlabel(r'$R_p$ [kpc]', fontsize=22)
ax.set_ylabel(r'$S_X$ $[erg\; s^{-1}\; kpc^{-2}]$', fontsize=22)  # ; keV^{-1}]$')# / (\Delta keV=40/1000)$ ')
ax.legend(frameon=True, markerscale=25, fontsize=20, loc=1, ncol=1)  # ,title=legend_title)
# p.xlim((E_min-0.01, E_max+0.01))
ax.set_yscale('log')
ax.set_xscale('log')
p.xticks(fontsize=22)
p.yticks(fontsize=22)
ax.set_xlim((4, 3000))
ax.set_ylim((2e34, 5e37))
fig.tight_layout()  # (rect=[0, 0.03, 1, 0.95])
fig.savefig(figure_name)
fig.clf()
