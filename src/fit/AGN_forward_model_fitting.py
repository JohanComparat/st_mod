"""
Model from Powell et al. 2022
https://ui.adsabs.harvard.edu/abs/2022ApJ...938...77P

Ananna et al. 2022
https://ui.adsabs.harvard.edu/abs/2022ApJS..261....9A/abstract

"""

import time
t0 = time.time()
import sys, os, glob
from astropy.table import Table
import numpy as np
from scipy.stats import norm

sys.path.append( os.path.join(os.environ['GIT_STMOD'], 'src') )
from models import AGN as GG

## option 1, Mhalo
#slope_Mhalo = 1.62
#norm_Mhalo = -12.38
#scatter_MBH_Mhalo = 0.4
#log10_MBH_halo = lambda log10_Mhalo, slope_Mhalo, norm_Mhalo : slope_Mhalo * log10_Mhalo + norm_Mhalo

def SMBH_SAR_Mstar_model(sim_data, params):
    """
    inputs:
    * sim_data : simulated data set predicted by Uchuu + UM
    * params : set of parameters
    """
    #slope_MBH, norm_MBH, scatter_MBH, lEs, delta1, epsilonL = params

    slope_MBH       = 0.51
    norm_MBH        = 7.55
    scatter_MBH = 0.43
    lEs    = 10**(-1.338)
    delta1 = 0.38
    epsilonL = 2.260 # need to be positive

    N_objects = len(sim_data)

    # option 2, Mstar
    #slope_Mstar       = 0.51
    #norm_Mstar        = 7.55
    #scatter_MBH_Mstar = 0.43
    log10_MBH_star = lambda log10_Mstar, slope_MBH, norm_MBH : slope_MBH * log10_Mstar + norm_MBH
    def scatter_MBH(NN): return norm.rvs(loc=0.0, scale=scatter_MBH, size=NN)

    log10_MBH_values = log10_MBH_star( np.log10(sim_data['obs_sm']), slope_MBH, norm_MBH ) + scatter_MBH(N_objects)

    # distribution of eddington ratios
    # draw a set of values from
    #lEs    = 10**(-1.338)
    #delta1 = 0.38
    #epsilonL = 2.260 # need to be positive
    lambda dNdlE = lE, lEs, delta1, delta2 : 1./( (lE/lEs)**delta1 + (lE/lEs)**(delta1+epsilonL))
    # normalize the function
    log10_lEDD_values = DRAWS from dNdlE

    LXhard = np.log10( 1.26 * 1e38 ) + log10_MBH_values + log10_lEDD_values - np.log10( 8 ) # s-1 erg

    # predict XLF in each redshift bin
    # compare to DATA

    # predict logNlogS by concatenating all bins
    # compare to DATA



z_dir = sys.argv[1]
#LC_dir = sys.argv[2] # 'FullSky'
LC_dir = 'LC0020'
scatter_0 = float(sys.argv[3])
f_sat = 0.0

C_AGN = GG.AGN(z_dir, LC_dir=LC_dir, scatter_0 = scatter_0, f_sat = f_sat )

print(len(C_AGN.p_2_catalogues), 'catalogues to loop over', ', Dt=', time.time()-t0, 'seconds')
is_cat = np.array([os.path.isfile(el) for el in C_AGN.p_2_catalogues])
for p_2_catalogue in C_AGN.p_2_catalogues[is_cat]:
    print('='*100)
    ##
    ## retrieve the area subtended by the catalog, skip step if area=0
    ##
    replication_dir = p_2_catalogue.split('/')[7]
    ix = float(replication_dir.split('_')[1])
    iy = float(replication_dir.split('_')[2])
    iz = float(replication_dir.split('_')[3])
    #print(ix, iy, iz)
    selection_area = ( C_AGN.LC_MetaData['jx'] == ix ) & ( C_AGN.LC_MetaData['jy'] == iy ) & ( C_AGN.LC_MetaData['jz'] == iz )
    #print(C_AGN.LC_MetaData[selection_area])
    sky_frac_Dmin = C_AGN.LC_MetaData['area_DC_min'][selection_area].sum() / C_AGN.LC_MetaData['area_DC_min'].sum()
    sky_frac_Dmax = C_AGN.LC_MetaData['area_DC_max'][selection_area].sum() / C_AGN.LC_MetaData['area_DC_max'].sum()
    C_AGN.sky_frac = sky_frac_Dmax # ( sky_frac_Dmin + sky_frac_Dmax ) / 2.
    print('sky fraction', C_AGN.sky_frac, sky_frac_Dmin, sky_frac_Dmax)
    if C_AGN.sky_frac==0.:
        print('area = 0 deg2')
        continue
    print('area = ', C_AGN.sky_frac, ' deg2')

    ##
    ## OPENS Galaxy light cone
    ##
    p_2_catalogue_out = os.path.join( os.path.dirname(p_2_catalogue), 'AGN_list_sigma_'+C_AGN.str_scatter_0+'_fsat_'+C_AGN.str_fsat+'.fits')
    print('opens', p_2_catalogue)
    t_sim_gal = Table.read(p_2_catalogue)
    t_sim_gal['ID_glist'] = np.arange(len(t_sim_gal))
    ZZ = t_sim_gal['redshift_R']
    C_AGN.z_min = np.min(ZZ)
    C_AGN.z_max = np.max(ZZ)
    C_AGN.z_mean = np.mean(ZZ)
    print(C_AGN.z_min,'<z<',C_AGN.z_max)
    #print(t_sim_gal.info())
    t_sim_gal.remove_columns([
                'x',
                'y',
                'z',
                'Mvir',
                'icl',
                'id',
                #'obs_sfr',
                #'obs_sm',
                'obs_uv',
                'sfr',
                'sm',
                #'upid',
                'A_UV',
                'Mpeak',
                'Vmax_Mpeak',
                'desc_id',
                'vmax',
                'vx',
                'vy',
                'vz',
                'RA',
                'DEC',
                'g_lat',
                'g_lon',
                'ecl_lat',
                'ecl_lon',
                'redshift_R',
                #'redshift_S',
                #'dL',
                #'nH',
                #'ebv'
                ])
    #print(t_sim_gal.info())

    ##
    ## retrieves tabulated AGN FILE
    ##
    C_AGN.get_tabulated_AGN()

    ##
    ## applies f_sat and duty cycle
    ##
    C_AGN.get_z_mass_id(t_sim_gal)
    ## scatter parameter is used here :
    C_AGN.abundance_matching()
    #
    C_AGN.get_obscured_fractions()
    C_AGN.compute_fluxes()
    C_AGN.AGN['agn_type'] = C_AGN.compute_agn_type(C_AGN.AGN['redshift_S'], C_AGN.AGN['LX_hard'], C_AGN.AGN['logNH'])
    C_AGN.compute_r_mag()

    t_out = C_AGN.AGN
    #print(t_out.info())
    t_out.remove_columns([
                'dL',
                'nH',
                'ebv'
                ])

    #print(t_out.info())

    t_out.write(p_2_catalogue_out, overwrite = True)
    print(p_2_catalogue_out, 'written', time.time()-t0)
