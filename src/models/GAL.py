"""
Model to paint a galaxy light cone with observed broad band photometry spectral energy distributions (SED) on the Uchuu UniverseMachine lightcone

Developped initially to create 4MOST mock catalogues with the collaboration of N. Taylor, C. Haines, A. Variu, J. Richard, A. Finoguenov.

input :
 - z_dir : name of the redshift slice (to retrieve the list of simulated galaxies, glist.fits files). Accessible via the UCHUU environment variable.
 - GAMA+KIDS+VIKING+LSDR9+RedMapper observed catalogue
 - COSMOS observed catalogue

output :
 - Kmatch_mags.fits : file containing SEDs corresponding to its glist.fits equivalent in the light cone (galaxy list)

"""
import time
t0 = time.time()

from sklearn.neighbors import BallTree
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import os, glob
from scipy.interpolate import interp1d
from scipy.stats import norm
from astropy.table import Table
import numpy as np

cosmoUCHUU = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
h = 0.6774
cosmo = cosmoUCHUU
zs = np.arange(0.0000001, 7.1, 0.001)
dm_itp = interp1d(zs, cosmo.distmod(zs).value)

class GAL:
    def __init__(self, z_dir, LC_dir='FullSky'):
        print('Initiates a GAL class model')
        print('------------------------------------------------')
        print('------------------------------------------------')
        self.z_dir = z_dir
        self.mean_z = int(z_dir[1])+int(z_dir[3:])/100.
        self.LC_dir = LC_dir
        print('directory:',self.z_dir, ', mean redshift=',self.mean_z)
        # get galaxy catalogues in the redshift slice of the light cone
        self.p_2_catalogues = np.array( glob.glob( os.path.join(os.environ['UCHUU'], self.LC_dir, z_dir, 'replication_*_*_*', 'glist.fits') ) )
        self.p_2_catalogues.sort()
        # pathes to observed catalogs
        self.path_2_SDSS = os.path.join(os.environ['DATA'], 'mpecl/comparat/data_s4/galaxy_catalogues/Ti20_SDSS_kdgroups','Ti20_full_MERGED.fits')
        self.path_2_KIDS = os.path.join(os.environ['DATA'], 'GAMA', 'G09.GAMADR4+LegacyDR9.galreference+RM.fits')
        self.path_2_GAMA = os.path.join(os.environ['DATA'], 'GAMA', 'forJohan.fits')
        self.path_2_COSMOS = os.path.join(os.environ['DATA'], 'COSMOS', 'photoz_vers2.0_010312.fits')
        self.LC_MetaData = Table.read( os.path.join(os.environ['UCHUU'], 'area_per_replica_'+self.z_dir+'.fits') )

    def construct_OBS_kdTree(self):
        """
        Constructs the kDTree with the observed redshift and stellar mass to match to.
        It creates the following attributes to the class :

        * self.min_Z : rescaling min/max for the variables going in the tree
        * self.max_Z : rescaling min/max for the variables going in the tree
        * self.min_K : rescaling min/max for the variables going in the tree
        * self.max_K : rescaling min/max for the variables going in the tree
        * self.z_01  : rescaled values in the Tree
        * self.k_01  : rescaled values in the Tree
        * self.t_ref : catalogue used for matching
        * self.Tree_Obs : tree created for the matching

        Next dev : split the populations into star foming and quiescent galaxies to assign SED.
        """
        if self.mean_z<0.8 :
            print('uses KIDS as reference catalog')
            t_ref = Table.read(self.path_2_KIDS)
            # quality cuts
            keep_kids = (8.9 - 2.5*np.log10(t_ref['flux_Kt'])>0) & (t_ref['z_peak']>0.01) & (t_ref['z_peak']<0.85)
            #
            t_ref = t_ref[keep_kids]
            t_ref['kmag'] = 8.9 - 2.5*np.log10(t_ref['flux_Kt'])
            t_ref['zmag'] = 8.9 - 2.5*np.log10(t_ref['flux_Zt'])
            t_ref['rmag'] = 8.9 - 2.5*np.log10(t_ref['flux_rt'])
            t_ref['imag'] = 8.9 - 2.5*np.log10(t_ref['flux_it'])
            t_ref['gmag'] = 8.9 - 2.5*np.log10(t_ref['flux_gt'])
            t_ref['umag'] = 8.9 - 2.5*np.log10(t_ref['flux_ut'])
            #
            # computing K-corrected K band absolute magnitudes
            # file below contains :
            # redshift, distmod, bandpass, kcorr_median, kcorr_16pc, kcorr_84pc = KCORR_DATA
            KCORR_DATA = np.loadtxt( os.path.join( os.environ['GIT_STMOD'], 'data', 'models','model_GAL', 'VISTA_Ks_kcorrections.txt'), unpack = True)
            kcorr_itp = interp1d(KCORR_DATA[0], KCORR_DATA[3])
            # kmag = ab_smag + distmod + bandpass + kcorrection
            # ab_smag = appmag - distmod - kcorr
            # is properly zero-centred; or at least zero-ish; the median +/- NMAD is 0.05 +/- 0.12.
            dm_values = dm_itp(t_ref['z_peak'].data.data)
            kmag_abs = t_ref['kmag'] - ( dm_values + kcorr_itp(t_ref['z_peak'].data.data) )
            t_ref['Kmag_abs'] = kmag_abs
            # rescale redshift variable
            self.min_Z = 0.85
            self.max_Z = 0.0
            # rescale MAG variables
            self.min_K = -28. #np.min(t_ref['Kmag_abs'])
            self.max_K = -15. #np.max(t_ref['Kmag_abs'])
            self.z_01 = (t_ref['z_peak'] - self.min_Z ) / ( self.max_Z - self.min_Z )
            self.k_01 = (t_ref['Kmag_abs'] - self.min_K ) / ( self.max_K - self.min_K )
            self.t_ref = t_ref

        if self.mean_z>=0.8 :
            print('uses COSMOS as reference catalog')

            t = Table.read(self.path_2_COSMOS)
            # quality cuts
            good = (t['photoz']>0.75 )&( t['photoz']< 6. ) & ( t['K']>0 )& ( t['MK']<0 )&( t['MK']>-40 )&( t['mass_med']<14 )&( t['mass_med']>6 )
            t_ref = t[good]
            t_ref['Kmag_abs'] = t_ref['MK']
            t_ref['umag'] = t_ref['U']
            t_ref['gmag'] = t_ref['G']
            t_ref['rmag'] = t_ref['R']
            t_ref['imag'] = t_ref['I']
            t_ref['zmag'] = t_ref['Z']
            t_ref['kmag'] = t_ref['K']
            # rescale redshift variable
            self.min_Z = 0.75
            self.max_Z = 6.0
            # rescale MAG variables
            self.min_K = -28. #np.min(t_ref['Kmag_abs'])
            self.max_K = -15. #np.max(t_ref['Kmag_abs'])
            # defines
            self.z_01 = (t_ref['photoz'] - self.min_Z ) / ( self.max_Z - self.min_Z )
            self.k_01 = (t_ref['Kmag_abs'] - self.min_K ) / ( self.max_K - self.min_K )
            self.t_ref = t_ref

        print(len(self.t_ref),'lines in reference catalogue')
        self.Tree_Obs = BallTree(np.transpose([self.z_01, self.k_01]))
        # placeholder for next DEV
        #self.Tree_Obs_SF = BallTree(np.transpose([self.z_01[SF], self.k_01[QU]]))
        #self.Tree_Obs_QU = BallTree(np.transpose([self.z_01[SF], self.k_01[QU]]))
        print('kd tree is constructed, t=', time.time()-t0, 'seconds')

    def add_kmag(self, t_sim_gal, str_redshift = 'redshift_S', str_stellar_mass = 'obs_sm'):
        """
        Adds absolute magnitude in the K-band to a catalogue containing redshift and stellar mass.

        input :

        * t_sim_gal : table containing a list of galaxies. Requires at least two columns redshift and stellar mass

        return :

        * the same table with a new column ('K_mag_abs') populated with absolute K-band magnitude

        """
        t_sim_gal ['K_mag_abs'] = -99.
        z0, z1, slope, OO, zpt, scatter = np.loadtxt(os.path.join( os.environ['GIT_STMOD'], 'data', 'models','model_GAL', 'logMs-MK-look-up-table-2parameters_2023Aug22.txt'), unpack=True)
        fun = lambda x, a, b : a * x + b
        for zmin, zmax, slope_i, OO_i, zpt_i, scatter_i in zip(z0, z1, slope, OO, zpt, scatter):
            #print(zmin, zmax, p_b)
            scatter_i = 0.15
            s_gal  = (t_sim_gal[str_redshift]>=zmin)  & (t_sim_gal[str_redshift]<=zmax)
            mag = lambda x : fun( x, slope_i, OO_i)
            if len(t_sim_gal ['K_mag_abs'][s_gal])>0:
                t_sim_gal ['K_mag_abs'][s_gal]  = mag(np.log10(t_sim_gal [str_stellar_mass][s_gal]))
                t_sim_gal ['K_mag_abs'][s_gal]+=norm.rvs(loc=0, scale=scatter_i, size=len(t_sim_gal[str_redshift][s_gal]))
        return t_sim_gal

    def add_magnitudes_direct_match_K(self, t_sim_gal, str_redshift = 'redshift_S', str_stellar_mass = 'obs_sm'):
        """
        Direct match using (redshift, K_mag_abs) between the observed catalogue and the simulated catalogue.

        input :

        * t_sim_gal : table containing a list of galaxies. Requires at least two columns redshift and stellar mass

        return :

        * a table with : line number of the input table (ID_glist), redshift, log10(stellar mass), absolute magnitude in the K-band, six observed magnitude columns (umag, gmag, rmag, imag, zmag, kmag), the distance in the matching procedure 'distance_match'.

        Next dev : split the populations into star foming and quiescent galaxies to assign SED.

        """
        print('starting direct match')
        # simulated quantities
        sim_redshift = t_sim_gal[str_redshift]
        sim_k_mag    = t_sim_gal['K_mag_abs']
        # rescale variables
        SIM_z_01      = ( sim_redshift - self.min_Z ) / ( self.max_Z - self.min_Z )
        SIM_k_01      = ( sim_k_mag - self.min_K ) / ( self.max_K - self.min_K )
        SIM_DATA = np.transpose([SIM_z_01, SIM_k_01])
        # search nearest neighbour in tree
        dist_out, ids_out = self.Tree_Obs.query(SIM_DATA, k=1, return_distance = True)
        kids_ID = np.arange(len(self.t_ref))
        ids = np.hstack((ids_out))
        dist_map = np.hstack(( dist_out ))
        id_to_map = kids_ID[ids]
        columns_to_add = np.array([ 'umag',
                                    'gmag',
                                    'rmag',
                                    'imag',
                                    'zmag',
                                    'kmag' ])

        print('outputing results')
        t_out = Table()
        t_out['ID_glist']   = np.arange(len(t_sim_gal))
        t_out[str_redshift] = sim_redshift
        t_out['log10_sm']   = np.log10(t_sim_gal [str_stellar_mass])
        t_out['K_mag_abs']  = sim_k_mag
        for el in columns_to_add :
            t_out[el] = self.t_ref[el][id_to_map]
        t_out['distance_match'] = dist_map
        return t_out


