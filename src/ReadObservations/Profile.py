"""
Library to read stacked profiles



"""
from astropy.table import Table, Column, hstack, vstack
import time
t0 = time.time()

import numpy as n
import numpy as np
import os, sys, glob
import astropy.io.fits as fits
from astropy.table import Table, Column
from scipy.stats import scoreatpercentile

import astropy.units as u
import astropy.constants as cc
from astropy.cosmology import FlatLambdaCDM
cosmoUNIT = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
h = 0.6774
L_box = 1000.0 / h
cosmo = cosmoUNIT
def Delta_virial_crit( z ):
	OZ = cosmo.Om(0)*(1.+z)**3 / ( cosmo.Om(z)*(1.+z)**3 + cosmo.Ode(z) ) ** 2
	x = OZ - 1.
	return (18*n.pi*n.pi+82.*x-39.*x*x)/(1.+x)

from scipy.optimize import curve_fit
from scipy.special import erf
def f_sat_mod(log10radius, log10r0, sigma, min_val, max_val) :
        return (0.5 + 0.5 * erf((log10radius - log10r0) / sigma))  * (max_val - min_val) + min_val

# color scheme
cs = {}
cs ['ALL_ANY'] = 'black'
cs ['CEN_ANY'] = 'darkgreen'
cs ['PS_CEN_ANY'] = 'red'
cs ['CEN_RS']  = 'firebrick'
cs ['CEN_BC']  = 'navy'
cs ['CEN_RS_hotGAS']  = 'goldenrod'
cs ['CEN_BC_hotGAS']  = 'aqua'
cs ['SAT_ANY'] = 'aquamarine'
cs ['SAT_BC']  = 'cornflowerblue'
cs ['SAT_RS']  = 'lightcoral'

cs ['ALL_Mstar']  = 'deeppink'
cs ['CEN_Mstar']  = 'fuchsia'
cs ['SAT_Mstar']  = 'plum'
cs ['ALL_Mhalo']  = 'forestgreen'
cs ['CEN_Mhalo']  = 'chartreuse'
cs ['SAT_Mhalo']  = 'turquoise'

cs['AGN'] = 'darkmagenta'
cs['XRB'] = 'hotpink'
cs['hotGAS'] = 'darkorange'
cs ['CEN_ANY_MS'] = 'lime'
cs['hotGAS_MS'] = 'moccasin'

cs ['FULL_2D'] = 'dodgerblue'
cs ['ISO_2D'] = 'fuchsia'


# = ['#FFC20A', '#5D3A9B', 'darkgreen', '#0C7BDC', 'silver']

from scipy.interpolate import interp1d
z_array = n.arange(0.01, 0.33, 0.0001)
dl2_4pi_itp = interp1d (z_array, (4*n.pi*cosmo.luminosity_distance(z_array).to(u.cm)**2).value )

def get_profile_SDSS_Ti21( p2_prof ) :
    """
    Function to read a profile from the SDSS stacks using the Tinker 2021 catalog
    input:
     - p2_prof: path to the profile
    """
    prefix = 'Ti20_SDSS_kdgroups'
    root_directory = os.environ['DATA_S4'] # /ptmp/joco/mpecl/comparat/data_s4/
    mergedCube_dir = os.path.join(root_directory, 'mergedCubes', prefix)
    #mergedCube_500E1000_dir  = os.path.join(root_directory, 'mergedCubes_500_1000', prefix)
    #mergedCube_1000E1200_dir = os.path.join(root_directory, 'mergedCubes_1000_2000', prefix)
    galaxy_directory = os.path.join(root_directory, 'galaxy_catalogues', prefix )
    simulated_directory = os.path.join(root_directory, 'simulated_catalogues', prefix )

    # simulated catalog
    s_cat = os.path.join(simulated_directory, os.path.basename(p2_prof)[:-21]+'-SIM.fits' )
    # galaxy catalog
    g_cat = os.path.join(galaxy_directory, os.path.basename(s_cat)[:-9]+'.fits')

    print(g_cat, s_cat, p2_prof)
    print(os.path.isfile(g_cat), os.path.isfile(s_cat), os.path.isfile(p2_prof) )
    if os.path.isfile(g_cat) and os.path.isfile(s_cat) and os.path.isfile(p2_prof):
        GAL = Table.read(g_cat)
        bn = os.path.basename(p2_prof).split('-')
        sim = Table.read(s_cat)
        t_out = Table.read(p2_prof)
        t_out['z_min'] = float(bn[2].split('_')[0])
        t_out['z_SEL'] = bn[2].split('_')[1]
        t_out['z_max'] = float(bn[2].split('_')[2])
        t_out['M_min'] = float(bn[3].split('_')[0])
        t_out['M_max'] = float(bn[3].split('_')[2])
        t_out['SFR'] = bn[1]
        t_out['CoS'] = bn[0]
        #
        t_out['redshift_mean'] =  n.mean(GAL['Z'])
        t_out['redshift_std'] =  n.std(GAL['Z'])
        # DD
        area_shell = ( t_out['x_up']**2 -  t_out['x_lo']**2 ) * n.pi
        t_out['bg'         ] = n.mean(t_out['dd_profile'    ][-2:])
        PROF_GAL_BG1 = t_out['bg'][0]
        sub_P = (t_out['dd_profile'    ])-PROF_GAL_BG1
        LXX = sub_P * area_shell
        sub_P_up = ((t_out['dd_profile'    ]+t_out['dd_profile_err']))-PROF_GAL_BG1
        LXX_up = sub_P_up * area_shell
        sub_P_low = ((t_out['dd_profile'    ]-t_out['dd_profile_err']))-PROF_GAL_BG1
        LXX_low = sub_P_low * area_shell

        t_out['BG'] = PROF_GAL_BG1
        t_out['area_shell'] = area_shell
        t_out['dd_profile_BGSUB'] = sub_P
        t_out['dd_profile_up_BGSUB'] = sub_P_up
        t_out['dd_profile_lo_BGSUB'] = sub_P_low
        # PS
        t_out['ps_profile_err']=t_out['ps_profile'    ]*t_out['ps_profile_NN']**-0.5
        sub_P = (t_out['ps_profile'    ])-PROF_GAL_BG1
        LXX = sub_P * area_shell
        sub_P_up = ((t_out['ps_profile'    ]+t_out['ps_profile_err']))-PROF_GAL_BG1
        LXX_up = sub_P_up * area_shell
        sub_P_low = ((t_out['ps_profile'    ]-t_out['ps_profile_err']))-PROF_GAL_BG1
        LXX_low = sub_P_low * area_shell
        #
        t_out['ps_profile_BGSUB'] = sub_P
        t_out['ps_profile_up_BGSUB'] = sub_P_up
        t_out['ps_profile_lo_BGSUB'] = sub_P_low
        #
        ps_profile = t_out['ps_profile'    ]
        profile = t_out['dd_profile'    ]
        PROF_GAL_BG1 = t_out['bg'][0]
        ps_profile_normed  = ((ps_profile - n.min(ps_profile)) / n.max(ps_profile-n.min(ps_profile))) * (n.max(profile[:3])-PROF_GAL_BG1) + PROF_GAL_BG1
        sub_PS = (ps_profile_normed-PROF_GAL_BG1)
        t_out['ps_profile_normed'] = ps_profile_normed
        t_out['ps_profile_normed_BGSUB'] = sub_PS
        LXX_ps_normed = n.sum(ps_profile_normed * area_shell)
        LXX_ps_normed_BG = n.sum(sub_PS * area_shell)
        t_out['file_ID']=os.path.basename(p2_prof)[:-5]
        # AGN
        AGN_LX_vals = sim['AGN_LX_soft'][(sim['AGN_LX_soft']>0)&(sim['AGN_agn_type']>20)]
        AGN_frac = 0.10 # len(AGN_LX_vals)*1./len(sim)
        N_total = len(AGN_LX_vals)/AGN_frac # len(AGN_LX_vals)*1./len(sim)
        AGN_LX_mean = n.log10(n.sum(10**(AGN_LX_vals))/N_total)
        t_out['AGN_Co19_mean_LX_all'] = AGN_LX_mean
        t_out['AGN_Co19_frac_in_simulation'] = len(AGN_LX_vals)*1./len(sim)
        t_out['AGN_Co19_frac'] = AGN_frac
        t_out['AGN_profile_normed'] = ps_profile_normed / LXX_ps_normed * 10**AGN_LX_mean
        t_out['AGN_profile_normed_BGSUB'] = sub_PS / LXX_ps_normed_BG * 10**AGN_LX_mean
        t_out['XRB_profile_normed'] = ps_profile_normed / LXX_ps_normed * 10**t_out['XRB_Ai17_med_LX_all'][0]
        t_out['XRB_profile_normed_BGSUB'] = sub_PS / LXX_ps_normed_BG * 10**t_out['XRB_Ai17_med_LX_all'][0]
        t_out['R200b']=n.mean(sim['R200b']) / (1+t_out['redshift_mean'])
        t_out['R200c']=n.mean(sim['R200c']) / (1+t_out['redshift_mean'])
        t_out['R500c']=n.mean(sim['R500c']) / (1+t_out['redshift_mean'])
        sim['delta_Mvir_crit'] = Delta_virial_crit(sim['redshift_R'])
        rho_crit = cosmo.critical_density(sim['redshift_R'])
        rho_bg = rho_crit * cosmo.Om(sim['redshift_R'])
        sim['Rvir'] = ( ( 3 * sim['Mvir'] * u.Msun ) / ( 4. * n.pi * sim['delta_Mvir_crit'] * rho_bg.to(u.Msun/u.kpc**3) ) ) ** ( 1. / 3. )
        t_out['Rvir']=n.mean(sim['Rvir'])
        #
        t_out['Mvir_mean']=n.mean(sim['Mvir'])
        t_out['M200b_mean']=n.mean(sim['M200b'])
        t_out['M200c_mean']=n.mean(sim['M200c'])
        t_out['M500c_mean']=n.mean(sim['M500c'])
        #
        t_out ['Mvir_std']=n.std(sim['Mvir'])
        t_out['M200b_std']=n.std(sim['M200b'])
        t_out['M200c_std']=n.std(sim['M200c'])
        t_out['M500c_std']=n.std(sim['M500c'])
        #
        t_out ['Mvir_min']=n.min(sim['Mvir'])
        t_out['M200b_min']=n.min(sim['M200b'])
        t_out['M200c_min']=n.min(sim['M200c'])
        t_out['M500c_min']=n.min(sim['M500c'])
        #
        t_out ['Mvir_max']=n.max(sim['Mvir'])
        t_out['M200b_max']=n.max(sim['M200b'])
        t_out['M200c_max']=n.max(sim['M200c'])
        t_out['M500c_max']=n.max(sim['M500c'])
        #
        t_out['GAL_stellarMass_mean'] =  n.mean(GAL['log10M_star'])
        t_out['GAL_stellarMass_std'] =  n.std(GAL['log10M_star'])
        t_out['GAL_Lgal_mean'] =  n.mean(GAL['L_gal'])
        t_out['GAL_Lgal_std']  =  n.std(GAL['L_gal'])
        t_out['GAL_Mhalo_mean'] =  n.mean(n.log10(GAL['M_halo']))
        t_out['GAL_Mhalo_std']   =  n.std(n.log10(GAL['M_halo']))
        t_out['GAL_Ltot_mean'] =  n.mean(GAL['L_tot'])
        t_out['GAL_Ltot_std'] =  n.std(GAL['L_tot'])
        t_out['GAL_Magg_mean'] =  n.mean(GAL['Mag_g'])
        t_out['GAL_Magg_std']  =   n.std(GAL['Mag_g'])
        t_out['GAL_Magr_mean'] =  n.mean(GAL['Mag_r'])
        t_out['GAL_Magr_std']  =   n.std(GAL['Mag_r'])
        t_out['GAL_Dn4000_mean'] =  n.mean(GAL['Dn4000'])
        t_out['GAL_Dn4000_std']  =   n.std(GAL['Dn4000'])
        t_out['GAL_concentration_mean'] =  n.mean(GAL['concentration'])
        t_out['GAL_concentration_std']  =   n.std(GAL['concentration'])
        # proper distance conversion
        t_out['dd_xb'] = t_out['dd_xb'] / (1+t_out['redshift_mean'])
        t_out['x_lo']  = t_out['x_lo']  / (1+t_out['redshift_mean'])
        t_out['x_up']  = t_out['x_up']  / (1+t_out['redshift_mean'])
        return t_out
    else:
        print('files missing')
