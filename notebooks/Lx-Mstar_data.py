import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as plt

#from xsingle_cubes_lib import *
import numpy as n
import os, sys, glob

import astropy.io.fits as fits
from astropy.table import Table, Column
from numpy import ndarray
from sklearn.neighbors import BallTree
from astropy.coordinates import SkyCoord
import astropy.units as u
import time
from astropy.wcs import WCS
import matplotlib.pyplot as p
t0 = time.time()
from scipy.interpolate import interp1d
from scipy.signal import convolve
from astropy.cosmology import FlatLambdaCDM
from scipy.interpolate import RegularGridInterpolator
#from lmfit import Minimizer, Parameters,report_fit,Model
#import corner
from scipy.optimize import minimize
#import emcee
from IPython.display import display, Math
from colossus.cosmology import cosmology
cosmology.setCosmology('planck18')
from colossus.halo import mass_defs
from colossus.halo import concentration

def M200toM500(M200m,z):
    c200m = concentration.concentration(M200m, '200m', z, model='ishiyama21')
    M500c, R500c, c500c = mass_defs.changeMassDefinition(M200m, c200m, z, '200m', '500c')
    return M500c#, R500c, c500c
cosmoUNIT = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
h = 0.6774
L_box = 1000.0 / h
cosmo = cosmoUNIT
deg_to_rad = n.pi/180.

colors_3lev=['#1b9e77', '#d95f02', '#7570b3']
colors_prb=['#984ea3', '#e41a1c', '#377eb8']
colors_fcs=['#ff7f00','#984ea3','gray']

colors_pp=['#6a3d9a','#b15928']
colors_bb=['#1f78b4','#8dd3c7']
colors_rr=['#e31a1c','#fdc086']

####CEN CGM-M####
Mminlist=10**n.array([10.0,10.5,11.0,11.25])
Mmaxlist=10**n.array([10.5,11.0,11.25,11.5])
Mmedlist=n.array([1.8e10,5.5e10,1.3e11,2.2e11])
Mmederr=[Mmedlist-Mminlist,Mmaxlist-Mmedlist]
Mmederrcen=[Mmedlist[1:]-Mminlist[1:],Mmaxlist[1:]-Mmedlist[1:]]

cgm=n.array([2.55119274e+39, 1.14164999e+40, 6.05290342e+40, 2.28216597e+41])
cgm_err=n.array([8.23905307e+38, 1.97541431e+39, 5.27424379e+39, 1.13620498e+40])


###Anderson CGM-M####
an_M=10**n.array([11.95197127020064, 11.852509953938636, 11.751643375751424, 11.652036328623103, 11.552637468446665, 11.45264527545736, 11.351809925312931, 11.25094334712572, 11.153324485387879, 11.05090691440914, 10.949426184713873, 10.85129726494392, 10.74895255939834, 10.651396153746063, 10.553298462018892, 10.454815624430738, 10.350628464360996, 10.250594633981315, 10.152580217034897])
an_M_err=10**n.array([[11.898956462903687, 11.800681812267417, 11.704135113331771, 11.603934733390584, 11.502277044786219, 11.403846253936035, 11.303937335727483, 11.207536367658157, 11.102152132615089, 11.003388242641892, 10.90118926796263, 10.80143648996799, 10.702370728914564, 10.60262836026752, 10.502074062508132, 10.401249121711297, 10.30072605199469, 10.19815234080204, 10.100242017331563],[12.000645379550836, 11.898956462903687, 11.800681812267417, 11.704135113331771, 11.603934733390584, 11.502277044786219, 11.403846253936035, 11.303937335727483, 11.207536367658157, 11.102152132615089, 11.003388242641892, 10.90118926796263, 10.80143648996799, 10.702370728914564, 10.60262836026752,10.502074062508132,10.401249121711297,10.30072605199469,10.19815234080204]])
an_Lx = 10**n.array([43.50882052647933, 43.25142273715985, 43.181604297065746, 42.76599635992365, 42.469676388334, 42.06212988857815, 41.57141208327784, 41.24408931504417, 40.985555111643805, 40.56341279353664, 40.28467172726062, 0.1, 39.28061437386248, 39.92016691081813, 0.1,0.1,39.62789541439162, 39.06885071247836, 39.74991787632619])
an_Lx_err=10**n.array([[43.3397434190083, 43.12879655524482, 43.08405025081014, 42.71336618280286, 42.41402761131087, 42.00264571403205, 41.49949837972211, 41.18265192879655, 40.89410929107294, 40.031464464864385, 39.812562702534734, 38.525183113597016, 38.51353486926799, 39.113881120433255, 38.528840946419855, 38.51662449505039, 38.8290318284725, 38.524508367736495, 38.93074088871133],[43.69462422870333, 43.3639987570471, 43.28320681848449, 42.82249744750744, 42.53470058152439, 42.130989479291514, 41.67940693390154, 41.35062813512673, 41.117024015625695, 41.083322235539576, 40.76182358947041, 38, 38.50398188840059, 40.71224752519199, 38.71776978736627, 38.69575176454921, 40.44810227726728, 39.92805078350424, 40.617392462378476]])
an_M_err=n.abs(n.array(an_M_err)-an_M)
an_Lx_err=n.abs(n.array(an_Lx_err)-an_Lx)

###CEN+Anderson>11###
Mmed_addan=n.concatenate((n.log10(Mmedlist),n.log10(an_M[:9])))
sort_Mid=n.argsort(Mmed_addan)
Mmed_addan=Mmed_addan[sort_Mid]
cgm_addan=n.concatenate((n.log10(cgm),n.log10(an_Lx[:9])))
cgm_addan=cgm_addan[sort_Mid]
cgmerr_addan=n.concatenate((n.log10(cgm_err+cgm)-n.log10(cgm),n.log10(an_Lx[:9]+an_Lx_err[0,:9])-n.log10(an_Lx[:9])))
cgmerr_addan=cgmerr_addan[sort_Mid]


fig, ax = p.subplots(figsize=[8,8])
p.title(r"", fontsize=22)
p.tick_params(which='both', length=10, width=2)
frame_thickness = 1.5  # Adjust the frame thickness as needed
for side in ['top', 'right', 'bottom', 'left']:
    ax.spines[side].set_linewidth(frame_thickness)

#ax.errorbar(x_data, y_data, yerr=y_derr, color='r',ls='',lw=4)
ax.errorbar(Mmedlist,cgm,xerr=Mmederr,yerr=cgm_err,lw=3, ls='',color=colors_3lev[2],fmt='o',markersize=8,label=r"$L_{\mathrm{X, CGM}}$, this work")
ax.errorbar(an_M,an_Lx,xerr=an_M_err,yerr=an_Lx_err,lw=2,ls='',marker='*', markersize=5,alpha=0.7, color=colors_pp[1], label='$L_{\mathrm{X, CGM}}$, Anderson et al. 2015')

ax.set_xlabel(r'$M_*$ [M$_\odot$]', fontsize=22)
ax.set_ylabel(r'$L_{\mathrm{X}}$ [erg/s]', fontsize=22)  # ; keV^{-1}]$')# / (\Delta keV=40/1000)$ ')
ax.legend(fontsize=19, loc=2, ncol=1,labelspacing=0.3,borderaxespad=0.1,edgecolor='white')  # ,title=legend_title)
# p.xlim((E_min-0.01, E_max+0.01))
ax.set_yscale('log')
ax.set_xscale('log')
p.xticks(fontsize=22)
p.yticks(fontsize=22)
p.grid()
ax.set_xlim((1e10, 1e12))
ax.set_ylim((1e38, 1e44))
fig.tight_layout()  # (rect=[0, 0.03, 1, 0.95])
#p.show()
p.savefig('test.png')
fig.clf()



###Anderson CGM-M####
an_M=n.array([11.95197127020064, 11.852509953938636, 11.751643375751424, 11.652036328623103, 11.552637468446665, 11.45264527545736, 11.351809925312931, 11.25094334712572, 11.153324485387879, 11.05090691440914, 10.949426184713873, 10.85129726494392, 10.74895255939834, 10.651396153746063, 10.553298462018892, 10.454815624430738, 10.350628464360996, 10.250594633981315, 10.152580217034897])
an_M_err=n.array([[11.898956462903687, 11.800681812267417, 11.704135113331771, 11.603934733390584, 11.502277044786219, 11.403846253936035, 11.303937335727483, 11.207536367658157, 11.102152132615089, 11.003388242641892, 10.90118926796263, 10.80143648996799, 10.702370728914564, 10.60262836026752, 10.502074062508132, 10.401249121711297, 10.30072605199469, 10.19815234080204, 10.100242017331563],[12.000645379550836, 11.898956462903687, 11.800681812267417, 11.704135113331771, 11.603934733390584, 11.502277044786219, 11.403846253936035, 11.303937335727483, 11.207536367658157, 11.102152132615089, 11.003388242641892, 10.90118926796263, 10.80143648996799, 10.702370728914564, 10.60262836026752,10.502074062508132,10.401249121711297,10.30072605199469,10.19815234080204]])
an_Lx = n.array([43.50882052647933, 43.25142273715985, 43.181604297065746, 42.76599635992365, 42.469676388334, 42.06212988857815, 41.57141208327784, 41.24408931504417, 40.985555111643805, 40.56341279353664, 40.28467172726062, 0.1, 39.28061437386248, 39.92016691081813, 0.1,0.1,39.62789541439162, 39.06885071247836, 39.74991787632619])
an_Lx_err=n.array([[43.3397434190083, 43.12879655524482, 43.08405025081014, 42.71336618280286, 42.41402761131087, 42.00264571403205, 41.49949837972211, 41.18265192879655, 40.89410929107294, 40.031464464864385, 39.812562702534734, 38.525183113597016, 38.51353486926799, 39.113881120433255, 38.528840946419855, 38.51662449505039, 38.8290318284725, 38.524508367736495, 38.93074088871133],[43.69462422870333, 43.3639987570471, 43.28320681848449, 42.82249744750744, 42.53470058152439, 42.130989479291514, 41.67940693390154, 41.35062813512673, 41.117024015625695, 41.083322235539576, 40.76182358947041, 38, 38.50398188840059, 40.71224752519199, 38.71776978736627, 38.69575176454921, 40.44810227726728, 39.92805078350424, 40.617392462378476]])

Mminlist=n.array([10.0,10.5,11.0,11.25])
Mmaxlist=n.array([10.5,11.0,11.25,11.5])
Mmedlist=n.array([1.8e10,5.5e10,1.3e11,2.2e11])
Mmederr=Mmaxlist-Mminlist
Mmederrcen=[Mmedlist[1:]-Mminlist[1:],Mmaxlist[1:]-Mmedlist[1:]]

cgm=n.array([2.55119274e+39, 1.14164999e+40, 6.05290342e+40, 2.28216597e+41])
cgm_err=n.array([8.23905307e+38, 1.97541431e+39, 5.27424379e+39, 1.13620498e+40])

#an_dy = abs((an_Lx_err[1]-an_Lx_err[0])/an_Lx)#/2.
s_an = (an_M>11.0)
x_data       =  n.hstack(( n.log10(Mmedlist)        , an_M         [s_an]))
x_data_err   =  n.hstack(( Mmederr        , an_M_err[1][s_an]-an_M_err[0][s_an]))
y_data    =  n.hstack(( n.log10(cgm)             , an_Lx        [s_an]))
y_data_up =  n.hstack(( n.log10(cgm + cgm_err)   , an_Lx_err[1] [s_an]  ))
y_data_lo =  n.hstack(( n.log10(cgm - cgm_err)   , an_Lx_err[0] [s_an]  ))
y_derr=abs(y_data_up-y_data_lo)/2.
s1 = (y_data>38) & (x_data<11.8)




fig, ax = p.subplots(figsize=[8,8])
p.title(r"", fontsize=22)
p.tick_params(which='both', length=10, width=2)
frame_thickness = 1.5  # Adjust the frame thickness as needed
for side in ['top', 'right', 'bottom', 'left']:
    ax.spines[side].set_linewidth(frame_thickness)

#ax.errorbar(x_data, y_data, xerr=x_data_err, yerr=y_derr, color='r',ls='',lw=4)
ax.errorbar(x_data[s1], y_data[s1], xerr=x_data_err[s1], yerr=y_derr[s1], color='b',ls='',lw=2)
#ax.errorbar(Mmedlist,cgm,xerr=Mmederr,yerr=cgm_err,lw=3, ls='',color=colors_3lev[2],fmt='o',markersize=8,label=r"$L_{\mathrm{X, CGM}}$, this work")
#ax.errorbar(an_M,an_Lx,xerr=an_M_err,yerr=an_Lx_err,lw=2,ls='',marker='*', markersize=5,alpha=0.7, color=colors_pp[1], label='$L_{\mathrm{X, CGM}}$, Anderson et al. 2015')

ax.set_xlabel(r'$M_*$ [M$_\odot$]', fontsize=22)
ax.set_ylabel(r'$L_{\mathrm{X}}$ [erg/s]', fontsize=22)  # ; keV^{-1}]$')# / (\Delta keV=40/1000)$ ')
ax.legend(fontsize=19, loc=2, ncol=1,labelspacing=0.3,borderaxespad=0.1,edgecolor='white')  # ,title=legend_title)
# p.xlim((E_min-0.01, E_max+0.01))
#ax.set_yscale('log')
#ax.set_xscale('log')
p.xticks(fontsize=22)
p.yticks(fontsize=22)
p.grid()
#ax.set_xlim((1e10, 1e12))
#ax.set_ylim((38, 45))
fig.tight_layout()  # (rect=[0, 0.03, 1, 0.95])
#p.show()
p.savefig('test2.png')
fig.clf()


n.savetxt('LxMs.ascii', n.transpose([x_data[s1], y_data[s1], x_data_err[s1], y_derr[s1]]))
