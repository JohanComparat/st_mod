"""
Validation of the AGN model
Figures :
 * stellar mass function
 * luminosity functions
 * logN-logS
"""
import time
t0 = time.time()
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as plt
import sys, os, glob
from scipy.integrate import quad

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

#import extinction

from scipy.interpolate import interp1d

from astropy.table import Table

import numpy as np

print('Tabulates logNlogS X-ray')
print('------------------------------------------------')
print('------------------------------------------------')


cosmoUCHUU = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
h = 0.6774
cosmo = cosmoUCHUU
zs = np.arange(0.0000001, 7.1, 0.001)
dl_itp = interp1d(zs, cosmo.luminosity_distance(zs).to(u.cm).value)

nl = lambda sel : len(sel.nonzero()[0])
zs = np.arange(0.0000001, 7.1, 0.001)
dm_itp = interp1d(zs, cosmo.distmod(zs).value)

#
# AGN MODEL
#
sys.path.append( os.path.join(os.environ['GIT_STMOD'], 'src') )
from models import AGN as GG

z_dir = sys.argv[1]
LC_dir = sys.argv[2] # 'FullSky'
C_AGN = GG.AGN(z_dir, LC_dir=LC_dir)
#LC_dir='FullSky'

#str_scatter_0 = '0.8'
#str_fsat = '8'


##

validation_dir           = os.path.join(os.environ['GIT_STMOD'], 'data', 'validation','validation_AGN')
validation_dir_Stell     = os.path.join(validation_dir, 'StellarMassFunction', z_dir)
validation_dir_XLF_hard  = os.path.join(validation_dir, 'XrayhardLuminosityFunction', z_dir)
validation_dir_XLF_soft  = os.path.join(validation_dir, 'XraySoftLuminosityFunction', z_dir)
validation_dir_lNlS_hard = os.path.join(validation_dir, 'XrayhardLogNlogS', z_dir)
validation_dir_lNlS_soft = os.path.join(validation_dir, 'XraySoftLogNlogS', z_dir)
validation_dir_Rband     = os.path.join(validation_dir, 'RbandLuminosityFunction', z_dir)
validation_dir_LSAR     = os.path.join(validation_dir, 'XrayLSAR', z_dir)
validation_dir_DC     = os.path.join(validation_dir, 'DutyCycle', z_dir)
#validation_dir_Kband     = os.path.join(validation_dir, 'KbandLuminosityFunction')

os.system('mkdir -p ' + validation_dir            )
os.system('mkdir -p ' + validation_dir_Stell      )
os.system('mkdir -p ' + validation_dir_XLF_hard   )
os.system('mkdir -p ' + validation_dir_XLF_soft   )
os.system('mkdir -p ' + validation_dir_lNlS_hard  )
os.system('mkdir -p ' + validation_dir_lNlS_soft  )
os.system('mkdir -p ' + validation_dir_Rband      )
os.system('mkdir -p ' + validation_dir_LSAR      )
os.system('mkdir -p ' + validation_dir_DC      )

##

#f_duty = C_AGN.f_duty


##
#
#
# tabulate smfs predicted by simulation
#
#
##
fdex = 0.05
fbins = np.arange(-20, -8, fdex)

##enough_area = (C_AGN.LC_MetaData['area_DC_max']>=0.5*np.max(C_AGN.LC_MetaData['area_DC_max'])) & (C_AGN.LC_MetaData['area_DC_max']>0)
##C_AGN.LC_MetaData['mean_area'] = (C_AGN.LC_MetaData['area_DC_min'] + C_AGN.LC_MetaData['area_DC_max']) / 2.
##small_difference_minmax_1 = ( C_AGN.LC_MetaData['area_DC_min'] / C_AGN.LC_MetaData['mean_area'] >= 0.8 ) & ( C_AGN.LC_MetaData['area_DC_min'] / C_AGN.LC_MetaData['mean_area'] <= 1.2 )
##small_difference_minmax_2 = ( C_AGN.LC_MetaData['area_DC_max'] / C_AGN.LC_MetaData['mean_area'] >= 0.8 ) & ( C_AGN.LC_MetaData['area_DC_max'] / C_AGN.LC_MetaData['mean_area'] <= 1.2 )

for meta in C_AGN.LC_MetaData:#[(enough_area)&(small_difference_minmax_1)&(small_difference_minmax_2)]:
    #
    print(meta)
    # retrieve the resulting catalogues and meta data
    p_2_catalogue = os.path.join(os.environ['UCHUU'], LC_dir, z_dir, 'replication_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']), 'glist.fits')
    #p_2_catal_MAG = os.path.join(os.environ['UCHUU'], LC_dir, z_dir, 'replication_'+str(meta['jx'])+'_'+str(meta['jy'])+'_'+str(meta['jz']), 'Kmatch_mags.fits')
    p_2_catal_AGNs = np.array( glob.glob( os.path.join( os.path.dirname(p_2_catalogue), 'AGN_list_sigma_*_fsat_*.fits' ) ) )
    #print(not os.path.isfile(p_2_catalogue))
    #print(not os.path.isfile(p_2_catal_MAG))
    #print(len(p_2_catal_AGNs)==0)
    print(p_2_catalogue)
    print(p_2_catal_AGNs)
    if len(p_2_catal_AGNs)==0 :
        continue
    #
    sky_frac_Dmax = meta['area_DC_max'] / C_AGN.LC_MetaData['area_DC_max'].sum()
    C_AGN.sky_frac = sky_frac_Dmax * meta['N_frac_'+LC_dir] # ( sky_frac_Dmin + sky_frac_Dmax ) / 2.
    #
    # create a dictionnary of catalogs
    AGN_cat_names = {}
    for p_2_catal_AGN in p_2_catal_AGNs:
        bn = os.path.basename(p_2_catal_AGN)[9:-5]
        AGN_cat_names[p_2_catal_AGN] = bn
    # loop over
    #p_2_catal_AGN = os.path.join( os.path.dirname(p_2_catalogue), 'AGN_list_sigma_'+str_scatter_0+'_fsat_'+str_fsat+'.fits' )
    AGNs = {}
    for p_2_catal_AGN in p_2_catal_AGNs:
        AGNs[AGN_cat_names[p_2_catal_AGN]] = Table.read(p_2_catal_AGN)

    #GAL = Table.read(p_2_catalogue)
    #MAG = Table.read(p_2_catal_MAG)
    #z_min, z_max = np.min(GAL['redshift_S']), np.max(GAL['redshift_S'])
    #print('z_min, z_max=', z_min, z_max)
    area_mock =  C_AGN.sky_frac * 129600 / np.pi
    #z_mean = np.mean(GAL['redshift_S'])
    # galaxy data
    #logm_gal = np.log10(GAL['sm'])
    #z_gal = GAL['redshift_S']
    #rmag_gal = MAG['rmag']
    #kmagABS_gal = MAG['K_mag_abs']

    for p_2_catal_AGN in p_2_catal_AGNs:
        logNlogS_out = os.path.join(os.path.dirname(p_2_catal_AGN), 'logNlogS_'+os.path.basename(p_2_catal_AGN))
        AGN = AGNs[AGN_cat_names[p_2_catal_AGN]]
        z = AGN['redshift_S']
        #dl2_cm = np.log10(dl_itp(z)) * 2 + np.log10(4*np.pi)
        fx_hard = AGN['FX_hard'] #- dl2_cm
        fx_soft = AGN['FX_soft'] #- dl2_cm
        t_out = Table()
        t_out['FX_lo'] = fbins[:-1]
        t_out['FX_hi'] = fbins[1:]
        t_out['N_hard'] = np.histogram(fx_hard, fbins)[0]
        t_out['N_soft'] = np.histogram(fx_soft, fbins)[0]
        t_out['area'] = area_mock
        print(t_out)
        t_out.write(logNlogS_out, overwrite=True)
        print(logNlogS_out, 'written')

