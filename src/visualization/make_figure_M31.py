"""
Makes a few figures illustrating the results for M31

python make_figure_M31.py

"""
import sys, os, glob
sys.path.append( os.path.join(os.environ['GIT_STMOD'], 'src') )
from ReadObservations import Profile as PR

import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as p


fig_dir = os.path.join( os.environ['GIT_STMOD'], 'results', 'stacked_galaxy_profile', 'Mstar_selection' )
os.system( 'mkdir -p ' + fig_dir )
p2_file_fit_out = os.path.join(fig_dir, 'MstarSel_ALL_SATfrac_params.fits')

# specZ stacks, SDSS - Tinker catalog

prefix = 'Ti20_SDSS_kdgroups'
root_directory = os.environ['DATA_S4'] # /ptmp/joco/mpecl/comparat/data_s4/
mergedCube_dir = os.path.join(root_directory, 'mergedCubes', prefix)
galaxy_directory = os.path.join(root_directory, 'galaxy_catalogues', prefix )
simulated_directory = os.path.join(root_directory, 'simulated_catalogues', prefix )

Z_SEL = 'Z3P' # sys.argv[1]
Z_SEL = 'Z2P' # sys.argv[1]
M_val, M_up = 11.0, 11.25
ALL_ANY = PR.get_profile_SDSS_Ti21 (glob.glob(os.path.join(mergedCube_dir,'ALL-ANY-0.01_'+Z_SEL+'_*-'+str(M_val)+'_Mstar_'+str(M_up)+'_STACKEDprofiles.fits' ))[0] )
CEN_ANY = PR.get_profile_SDSS_Ti21 (glob.glob(os.path.join(mergedCube_dir,'CEN-ANY-0.01_'+Z_SEL+'_*-'+str(M_val)+'_Mstar_'+str(M_up)+'_STACKEDprofiles.fits' ))[0] )
#CEN_RS = PR.get_profile_SDSS_Ti21 (glob.glob(os.path.join(mergedCube_dir,'CEN-RS-0.01_'+Z_SEL+'_*-'+str(M_val)+'_Mstar_'+str(M_up)+'_STACKEDprofiles.fits' ))[0] )
#CEN_BC = PR.get_profile_SDSS_Ti21 (glob.glob(os.path.join(mergedCube_dir,'CEN-BC-0.01_'+Z_SEL+'_*-'+str(M_val)+'_Mstar_'+str(M_up)+'_STACKEDprofiles.fits' ))[0] )
SAT_ANY = PR.get_profile_SDSS_Ti21 (glob.glob(os.path.join(mergedCube_dir,'SAT-ANY-0.01_'+Z_SEL+'_*-'+str(M_val)+'_Mstar_'+str(M_up)+'_STACKEDprofiles.fits' ))[0] )
#SAT_RS = PR.get_profile_SDSS_Ti21 (glob.glob(os.path.join(mergedCube_dir,'SAT-RS-0.01_'+Z_SEL+'_*-'+str(M_val)+'_Mstar_'+str(M_up)+'_STACKEDprofiles.fits' ))[0] )
#SAT_BC = PR.get_profile_SDSS_Ti21 (glob.glob(os.path.join(mergedCube_dir,'SAT-BC-0.01_'+Z_SEL+'_*-'+str(M_val)+'_Mstar_'+str(M_up)+'_STACKEDprofiles.fits' ))[0] )


colors_2D_gal = ['pink','purple',p.cm.tab20b(0)]
colors_2D_sf = [p.cm.tab20b(2),p.cm.tab20c(0),p.cm.tab20b(0)]
colors_2D_qu = ['pink','orange','darkred']

bins = np.hstack(( 0., 10., 30., 50., 75., 100.,150, 200.,250,300.,350,400.,500.,600.,700.,800.,900., 1000. ,1100,1200,1300,1400,1600,1800,2000,2200,2400,2600,2800,3000))#l_bins # np.hstack(( 0., 10**np.arange(1.2, np.log10(3000.), 0.2), 3000. ))

profile_path = glob.glob(os.path.join(root_directory, 'mergedCubes', 'DR9_Ebins_s4','data','profile','profile_*_nomask_2D_GAL.txt'))
sort=np.argsort(np.array([os.path.basename(i).split('_')[1] for i in profile_path]).astype(float))
profile_path_GAL_nomask=np.array(profile_path)[sort]

profile_path = glob.glob(os.path.join(root_directory, 'mergedCubes', 'DR9_Ebins_s4','data','profile','profile_*_nomask_full_GAL.txt'))
sort=np.argsort(np.array([os.path.basename(i).split('_')[1] for i in profile_path]).astype(float))
profile_path_GAL_full_nomask=np.array(profile_path)[sort]

def plot_profile_line(p_2_pro,colors_2D_gal,i):
    Profile = np.transpose(np.loadtxt(p_2_pro))
    xb = Profile[0]
    totalerr = Profile[2]
    BG_new_mean = Profile[3][0]
    BG_new_std = Profile[4][0]
    title_list = os.path.basename(p_2_pro).split('_')
    # [xb,Ix_mean,totalerr, np.ones_like(xb)*BG_new_mean, np.ones_like(xb)*BG_new_std,ps_profile_normed])
    p.errorbar(xb, Profile[1], yerr=totalerr, xerr=[xb - bins[:-1], bins[1:] - xb], lw=3, ls='',
                color=colors_2D_gal[i], zorder=20, label=title_list[-2])
    p.fill_between(xb, Profile[1] - totalerr, Profile[1] + totalerr, color=colors_2D_gal[i], alpha=0.3, zorder=5)
    p.axhline(BG_new_mean, lw=3, color=colors_2D_gal[i], ls='dashed')
    p.fill_between(xb, BG_new_mean - BG_new_std, BG_new_mean + BG_new_std, color=colors_2D_gal[i], alpha=0.3, zorder=3)
    p.step(bins[1:], Profile[5], lw=2, zorder=3, ls=':', color=colors_2D_gal[i])
    p.axvline(Profile[-1][0], color=colors_2D_gal[i], ls='--')

def plot_profile_line_sub(p_2_pro,colors_2D_gal,i):
    Profile = np.transpose(np.loadtxt(p_2_pro))
    xb = Profile[0]
    totalerr = Profile[2]
    BG_new_mean = Profile[3][0]
    BG_new_std = Profile[4][0]
    title_list = os.path.basename(p_2_pro).split('_')
    # [xb,Ix_mean,totalerr, np.ones_like(xb)*BG_new_mean, np.ones_like(xb)*BG_new_std,ps_profile_normed])
    sel_lim = (Profile[1] - BG_new_mean - totalerr > 0)
    p.errorbar(xb[sel_lim], Profile[1][sel_lim] - BG_new_mean, yerr=totalerr[sel_lim],
                xerr=[xb[sel_lim] - bins[:-1][sel_lim], bins[1:][sel_lim] - xb[sel_lim]], lw=2, ls='',
                color=colors_2D_gal[i], label=title_list[-2],
                zorder=5)
    for x, y, dy, wid, xup, xlow, y_ln, x_low in zip(xb, Profile[ 1] - BG_new_mean + totalerr, totalerr / 4,
                                                     bins[1:] - bins[:-1],
                                                     xb - bins[:-1], bins[1:] - xb, Profile[ 1] - BG_new_mean,
                                                     Profile[1] - BG_new_mean - totalerr):
        if x_low <= 0:
            p.errorbar(x, y, xerr=[[xup], [xlow]], yerr=dy, uplims=True, color=colors_2D_gal[i],
                        )
    p.step(bins[1:],  Profile[5] - BG_new_mean, lw=2, zorder=3, ls=':', color=colors_2D_gal[i])
    p.axhline(BG_new_mean / 100, lw=3, color=colors_2D_gal[i], alpha=0.7, ls='dashed')
    p.axvline( Profile[-1][0], color=colors_2D_gal[i], ls='--')

    kernel = ( Profile[5]-np.min( Profile[5]) ) / np.max( Profile[5]-np.min( Profile[5]))
    #xbk = xb[kernel>0.03]
    kernel = kernel[kernel>0.03]
    return kernel
#
#
# Measurement
#
#
p2_fig = os.path.join(fig_dir, 'M31_MstarSel_FULL_ISO_ALL_CEN_SAT_SBprofile.png')
p.figure(3, (5,4.5) )
# ALL
p.step(np.ones_like(ALL_ANY['dd_xb'])-1000, np.ones_like(ALL_ANY['ps_profile_normed'])-1000., lw=0.8 , color='k', where='mid', label='PSF profile' )
p.axvline(CEN_ANY['R500c'][0], ls='--', lw=2, c='k', label='$R_{500c}=$'+str(np.round(CEN_ANY['R500c'][0],1))+'kpc')
str_NG = ', N='+str(ALL_ANY['N_gal'][0])
str_MS = ', M$_{200m}$='+str(np.round(ALL_ANY['GAL_Mhalo_mean'][0],1))+'$\pm$'+str(np.round(ALL_ANY['GAL_Mhalo_std'][0],1))
N_all = ALL_ANY['N_gal'][0]
p.errorbar( ALL_ANY['dd_xb'], ALL_ANY['dd_profile'],
        yerr = ALL_ANY['dd_profile_err'],
        xerr = [ ALL_ANY['dd_xb'] - ALL_ANY['x_lo'], ALL_ANY['x_up'] - ALL_ANY['dd_xb']],
        lw=2, ls='', color=PR.cs['ALL_ANY'], label='ALL')# + str_NG + str_MS )
p.step(ALL_ANY['dd_xb'], ALL_ANY['ps_profile_normed'], lw=0.8 , color=PR.cs['ALL_ANY'], where='mid' )

#CEN
str_NG = ', N='+str(CEN_ANY['N_gal'][0])
str_MS = ', M$_{200m}$='+str(np.round(CEN_ANY['GAL_Mhalo_mean'][0],1))+'$\pm$'+str(np.round(CEN_ANY['GAL_Mhalo_std'][0],1))
N_cen = CEN_ANY['N_gal'][0]
#frac_signal = N_cen*1./N_all
p.errorbar( CEN_ANY['dd_xb'], CEN_ANY['dd_profile'],
        yerr = CEN_ANY['dd_profile_err'],
        xerr = [ CEN_ANY['dd_xb'] - CEN_ANY['x_lo'], CEN_ANY['x_up'] - CEN_ANY['dd_xb']],
        lw=2, ls='', color=PR.cs['CEN_ANY'], label='Cen')# + str_NG + str_MS )
p.step(CEN_ANY['dd_xb'], CEN_ANY['ps_profile_normed'], lw=0.8 , color=PR.cs['CEN_ANY'], where='mid' )

#SAT
N_sat = SAT_ANY['N_gal'][0]
#frac_signal = N_sat*1./N_all
str_NG = ', N='+str(SAT_ANY['N_gal'][0])
str_MS = ', M$_{200m}$='+str(np.round(SAT_ANY['GAL_Mhalo_mean'][0],1))+'$\pm$'+str(np.round(SAT_ANY['GAL_Mhalo_std'][0],1))
p.errorbar( SAT_ANY['dd_xb'], SAT_ANY['dd_profile'],
        yerr = SAT_ANY['dd_profile_err'],
        xerr = [ SAT_ANY['dd_xb'] - SAT_ANY['x_lo'], SAT_ANY['x_up'] - SAT_ANY['dd_xb']],
        lw=2, ls='', color=PR.cs['SAT_ANY'], label='Sat')# + str_NG + str_MS )
p.step(SAT_ANY['dd_xb'], SAT_ANY['ps_profile_normed'], lw=0.8 , color=PR.cs['SAT_ANY'], where='mid' )

# FULL & ISOLATED
plot_profile_line(profile_path_GAL_nomask[3],colors_2D_gal,0)
plot_profile_line(profile_path_GAL_full_nomask[3],colors_2D_gal,1)

p.axhline(1, lw=3, color=colors_2D_gal[1], ls='dashed', label='background')
p.axhline(1, lw=2, zorder=3, ls=':', color=colors_2D_gal[1], label='empirical PSF')

p.xlabel(r'$R_p$ [kpc]')
p.ylabel(r'$S_X$ $[erg\; s^{-1}\; kpc^{-2}]$')
p.legend(fontsize=10, loc=0, ncol=1)
p.xlim((4, 2000))
p.ylim((2e36, 6e37))
p.yscale('log')
p.xscale('log')
p.title(str(M_val)+'$<\log_{10}(M^*/M_\odot)<$'+str(M_up), fontsize=12)
p.tight_layout()
p.savefig( p2_fig )
p.clf()
print(p2_fig, 'written')


#
#
# Measurement BG subtracted
#
#
p2_fig = os.path.join(fig_dir, 'M31_MstarSel_FULL_ISO_ALL_CEN_SAT_BGsubSBprofile.png')
p.figure(4, (6.5,6.5) )
# ALL
#p.step(np.ones_like(ALL_ANY['dd_xb'])-1000, np.ones_like(ALL_ANY['ps_profile_normed'])-1000., lw=0.8 , color='k', where='mid', label='PSF profile' )
p.axhline(ALL_ANY['bg'][0]/300., ls='--', lw=2, c='grey', label='background/300')
p.axvline(CEN_ANY['R500c'][0], ls='--', lw=2, c='k', label='$R_{500c}=$'+str(np.round(CEN_ANY['R500c'][0],1))+'kpc')
str_NG = ', N='+str(ALL_ANY['N_gal'][0])
str_MS = ', M$_{200m}$='+str(np.round(ALL_ANY['GAL_Mhalo_mean'][0],1))+'$\pm$'+str(np.round(ALL_ANY['GAL_Mhalo_std'][0],1))
N_all = ALL_ANY['N_gal'][0]
p.errorbar( ALL_ANY['dd_xb'], ALL_ANY['dd_profile_BGSUB'],
        yerr = [ ALL_ANY['dd_profile_BGSUB'] - ALL_ANY['dd_profile_lo_BGSUB'], ALL_ANY['dd_profile_up_BGSUB'] - ALL_ANY['dd_profile_BGSUB'] ],
        xerr = [ ALL_ANY['dd_xb'] - ALL_ANY['x_lo'], ALL_ANY['x_up'] - ALL_ANY['dd_xb']],
        lw=2, ls='', color=PR.cs['ALL_ANY'], label='ALL')# + str_NG + str_MS )
p.step(ALL_ANY['dd_xb'], ALL_ANY['ps_profile_normed_BGSUB'], lw=0.8 , color=PR.cs['ALL_ANY'], where='mid')#, label='PSF profile' )

# CEN ANY
str_NG = ', N='+str(CEN_ANY['N_gal'][0])
str_MS = ', M$_{200m}$='+str(np.round(CEN_ANY['GAL_Mhalo_mean'][0],1))+'$\pm$'+str(np.round(CEN_ANY['GAL_Mhalo_std'][0],1))
N_cen = CEN_ANY['N_gal'][0]
frac_signal = N_cen*1./N_all
p.errorbar( CEN_ANY['dd_xb'], CEN_ANY['dd_profile_BGSUB'],
        yerr = [ CEN_ANY['dd_profile_BGSUB'] - CEN_ANY['dd_profile_lo_BGSUB'], CEN_ANY['dd_profile_up_BGSUB'] - CEN_ANY['dd_profile_BGSUB'] ],
        xerr = [ CEN_ANY['dd_xb'] - CEN_ANY['x_lo'], CEN_ANY['x_up'] - CEN_ANY['dd_xb']],
        lw=2, ls='', color=PR.cs['CEN_RS'], label='CEN intrinsic')# + str_NG + str_MS )
#p.step(CEN_ANY['dd_xb'], CEN_ANY['ps_profile_normed_BGSUB'], lw=0.8 , color=PR.cs['CEN_ANY'], where='mid', label='PSF profile' )

# CEN ANY
str_NG = ', N='+str(CEN_ANY['N_gal'][0])
str_MS = ', M$_{200m}$='+str(np.round(CEN_ANY['GAL_Mhalo_mean'][0],1))+'$\pm$'+str(np.round(CEN_ANY['GAL_Mhalo_std'][0],1))
N_cen = CEN_ANY['N_gal'][0]
frac_signal = N_cen*1./N_all
p.errorbar( CEN_ANY['dd_xb'], CEN_ANY['dd_profile_BGSUB']*frac_signal,
        yerr = [ CEN_ANY['dd_profile_BGSUB']*frac_signal - CEN_ANY['dd_profile_lo_BGSUB']*frac_signal, CEN_ANY['dd_profile_up_BGSUB']*frac_signal - CEN_ANY['dd_profile_BGSUB']*frac_signal ],
        xerr = [ CEN_ANY['dd_xb'] - CEN_ANY['x_lo'], CEN_ANY['x_up'] - CEN_ANY['dd_xb']],
        lw=2, ls='', color=PR.cs['CEN_ANY'], label='CEN frac')# + str_NG + str_MS )
#p.step(CEN_ANY['dd_xb'], CEN_ANY['ps_profile_normed_BGSUB'], lw=0.8 , color=PR.cs['CEN_ANY'], where='mid', label='PSF profile' )

# SAT ANY
N_sat = SAT_ANY['N_gal'][0]
frac_signal = N_sat*1./N_all
str_NG = ', N='+str(SAT_ANY['N_gal'][0])
str_MS = ', M$_{200m}$='+str(np.round(SAT_ANY['GAL_Mhalo_mean'][0],1))+'$\pm$'+str(np.round(SAT_ANY['GAL_Mhalo_std'][0],1))
p.errorbar( SAT_ANY['dd_xb'], SAT_ANY['dd_profile_BGSUB']*frac_signal,
        yerr = [ SAT_ANY['dd_profile_BGSUB']*frac_signal - SAT_ANY['dd_profile_lo_BGSUB']*frac_signal, SAT_ANY['dd_profile_up_BGSUB']*frac_signal - SAT_ANY['dd_profile_BGSUB']*frac_signal ],
        xerr = [ SAT_ANY['dd_xb'] - SAT_ANY['x_lo'], SAT_ANY['x_up'] - SAT_ANY['dd_xb']],
        lw=2, ls='', color=PR.cs['SAT_ANY'], label='SAT frac')# + str_NG + str_MS )
#p.step(SAT_ANY['dd_xb'], SAT_ANY['ps_profile_normed_BGSUB'], lw=0.8 , color=PR.cs['SAT_ANY'], where='mid', label='PSF profile' )

## FULL
#kernel = plot_profile_line_sub(profile_path_GAL_full_nomask[3],colors_2D_gal,1)


#y_value = np.convolve(ALL_ANY['dd_profile_BGSUB'], kernel/kernel.sum(), mode='same')
#y_value_up = np.convolve(ALL_ANY['dd_profile_up_BGSUB'] - ALL_ANY['dd_profile_BGSUB'], kernel/kernel.sum(), mode='same')
#y_value_lo = np.convolve(ALL_ANY['dd_profile_BGSUB'] - ALL_ANY['dd_profile_lo_BGSUB'], kernel/kernel.sum(), mode='same')
#p.errorbar( ALL_ANY['dd_xb'], y_value,
        #yerr = [ y_value_lo, y_value_up ],
        #xerr = [ ALL_ANY['dd_xb'] - ALL_ANY['x_lo'], ALL_ANY['x_up'] - ALL_ANY['dd_xb']],
        #lw=2, ls='', color=PR.cs['CEN_Mhalo'], label='ALL conv Z2P')# + str_NG + str_MS )

# ISOLATED
kernel = plot_profile_line_sub(profile_path_GAL_nomask[3],colors_2D_gal,0)

p.xlabel(r'$R_p$ [kpc]')
p.ylabel(r'$S_X$ $[erg\; s^{-1}\; kpc^{-2}]$')
p.legend(fontsize=10, loc=3, ncol=1)
p.xlim((4, 2000))
p.ylim((ALL_ANY['bg'][0]/10000., np.max(ALL_ANY['dd_profile_BGSUB'])*1.5))
p.yscale('log')
p.xscale('log')
p.title(str(M_val)+'$<\log_{10}(M^*/M_\odot)<$'+str(M_up), fontsize=12)
p.tight_layout()
p.savefig( p2_fig )
p.clf()
print(p2_fig, 'written')
