"""
Model to paint a dark matter halo light cone with X-ray properties related to the virialized hot gas

Follows the cluster model of Comparat et al. 2019 used and improved on by Liu et al. 2022 (eFEDs simulations), Seppi et al. 2022 (eRASS1 simulation going down to group scale)

input :
 - z_dir : name of the redshift slice (to retrieve the list of simulated galaxies, glist.fits files). Accessible via the UCHUU environment variable.
 -

output :
 - Clusters.fits : file containing cluster properties

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

from colossus.cosmology import cosmology
cosmology.setCosmology('planck18')
from colossus.halo import mass_so
from colossus.halo import mass_defs
from colossus.halo import concentration

class GAS:
    def __init__(self, z_dir, b_HS=0.8, logM500c_min=11., logFX_min=-18, LC_dir='FullSky'):
        print('Initiates a GAS class model')
        print('------------------------------------------------')
        print('------------------------------------------------')
        self.z_dir = z_dir
        self.mean_z = int(z_dir[1])+int(z_dir[3:])/100.
        self.LC_dir = LC_dir
        print('directory:',self.z_dir, ', mean redshift=',self.mean_z)
        # get halo + galaxy catalogues in the redshift slice of the light cone
        self.p_2_catalogues = np.array( glob.glob( os.path.join(os.environ['UCHUU'], self.LC_dir, z_dir, 'replication_*_*_*', 'glist.fits') ) )
        self.p_2_catalogues.sort()
        # meta data
        self.LC_MetaData = Table.read( os.path.join(os.environ['UCHUU'], 'area_per_replica_'+self.z_dir+'.fits') )
        self.b_HS = b_HS
        self.logM500c_min = logM500c_min
        self.logFX_min = logFX_min

    def prepare_halo_cat(self, t1):
        """
        Prepares the halo catalogue

        input:
            -   t1 : catalog of distince haloes

        output stored in the class object :
            -   self.CAT : catalogue of distinct haloes for which X-ray quantities will be computed
                - adds converted masses to 500c, 200c and 200b using mass-concentration relation from Ishiyama 2021 (for consistency)
            -   self.N_cat : number of haloes considered
            -   self.zmin, self.zmax : minimum and maximum redshifts
        """
        Mvir = t1['Mvir']
        zzz = t1['redshift_R']
        cvir = concentration.concentration(Mvir, 'vir', zzz, model = 'ishiyama21')
        t1['M200b'], t1['R200b'], t1['c200b'] = mass_defs.changeMassDefinition(Mvir, cvir, zzz.mean(), 'vir', '200m')
        t1['M200c'], t1['R200c'], t1['c200c'] = mass_defs.changeMassDefinition(Mvir, cvir, zzz.mean(), 'vir', '200c')
        t1['M500c'], t1['R500c'], t1['c500c'] = mass_defs.changeMassDefinition(Mvir, cvir, zzz.mean(), 'vir', '500c')
        self.CAT = t1[(t1['M500c']>10**self.logM500c_min)]
        self.CAT.remove_columns([
                'x',
                'y',
                'z',
                #'Mvir',
                'icl',
                'id',
                #'obs_sfr',
                #'obs_sm',
                'obs_uv',
                'sfr',
                'sm',
                'upid',
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
                'dL',
                'nH',
                'ebv'])
        t1=0
        self.N_cat = len(self.CAT )
        self.zmin = np.min(self.CAT['redshift_S'] )
        self.zmax = np.max(self.CAT['redshift_S'] )


    def calc_lx(self, prof,kt,m5,z,fraction=1):
        """
        Compute the X-ray luminosity in the profile to be extended to 3x r500c.

        .. math::
            r_{500c} = \left(\\frac{3 M_{500c}}{ 4. \pi 500 \\rho_c(z)  }\\right)^{1/3} [ cm ]


        * profile\_emission = profile x rescale_factor
        * rescale\_factor = $\sqrt(kT/10.0) E^3(z)$
        * CF(kT) = cooling function, show the curve
        * L$_X(r)$ = $\Sigma_{<r}$( profile_emission $r_{500c}^2 2 \pi x CF(kT)$ Mpc=3.0856776e+24 dx )
        * L$_{500c}$ = L$_X$(1)

        """
        Mpc=3.0856776e+24
        msun=1.98892e33
        ez2 = cosmo.efunc(z)**2
        rhoc = cosmo.critical_density(z).value
        r500 = np.power(m5*msun/4.*3./np.pi/500./rhoc,1./3.)
        resfact = np.sqrt(kt/10.0)*np.power(ez2,3./2.)
        prof_em = prof * resfact # emission integral
        tlambda = np.interp(kt,self.coolfunc[:,0],self.coolfunc[:,1]) # cooling function
        dx = np.empty(len(self.xgrid_ext))
        dx[0]=self.xgrid_ext[0]
        dx[1:len(self.xgrid_ext)]=(np.roll(self.xgrid_ext,-1)-self.xgrid_ext)[:len(self.xgrid_ext)-1]
        #print(prof_em*self.xgrid_ext*r500**2*2.*np.pi*tlambda*Mpc*dx)
        lxcum = np.cumsum(prof_em*self.xgrid_ext*r500**2*2.*np.pi*tlambda*Mpc*dx) # riemann integral
        lx_500=np.interp(fraction,self.xgrid_ext,lxcum) # evaluated at R500
        return lx_500


    def draw_simulated_values(self, nsim = 100, zmin=0., zmax=0.1):
        """
        Draws simulated values from the covariance matrix of Comparat, Eckert et al. 2019

        Adds as attribute to the class the covariance matrix and its dependencies
            - self.path_2_cbp: path to the tabulated data
            - self.covor     : covariance matric
            - self.xgrid_ext : radial grid
            - self.mean_log  : mean parameters
            - self.coolfunc  : cooling function

        Then draws nsim simulated vectors and adds them as attribute to the class
            - self.allz      : redshifts
            - self.allkt     : temperatures
            - self.allm5     : M500c
            - self.profiles  : profiles
            - self.nsim2     : number of simulated vectros in the redshift range of interest
            - self.frac_r200c = R200c/R500c
            - self.frac_rVIR  = Rvir/R500c
            - self.alllx       : integrated LX up to R500c
            - self.alllx_R200c : integrated LX up to R200c
            - self.alllx_Rvir  : integrated LX up to Rvir

        """
        self.path_2_cbp= os.path.join(os.environ['GIT_STMOD'], 'data/models/model_GAS')
        self.covor     = np.loadtxt(os.path.join(self.path_2_cbp, 'covmat_xxl_hiflugcs_xcop.txt'))
        self.xgrid_ext = np.loadtxt(os.path.join(self.path_2_cbp, 'radial_binning.txt'))
        self.mean_log  = np.loadtxt(os.path.join(self.path_2_cbp, 'mean_pars.txt'))
        self.coolfunc  = np.loadtxt(os.path.join(self.path_2_cbp, 'coolfunc.dat'))

        profs=np.exp(np.random.multivariate_normal(self.mean_log,self.covor,size=nsim))

        allz_i  = profs[:,len(self.mean_log)-3]
        allkt_i = profs[:,len(self.mean_log)-1]
        allm5_i = profs[:,len(self.mean_log)-2]

        profiles_i = profs[:,:len(self.xgrid_ext)]
        # filters output to cover the relevant redshift range
        in_zbin = (allz_i<zmax+0.01)&(allz_i>zmin-0.01)

        self.allz     = allz_i    [in_zbin]
        self.allkt    = allkt_i   [in_zbin]
        self.allm5    = allm5_i   [in_zbin]
        self.profiles = profiles_i[in_zbin]
        self.nsim2    = len(self.allz)
        c500c = concentration.concentration(self.allm5, '500c', self.allz, model = 'ishiyama21')
        Mvir , Rvir , cvir  = mass_defs.changeMassDefinition(self.allm5, c500c, np.mean(self.allz), '500c', 'vir')
        M200c, R200c, c200c = mass_defs.changeMassDefinition(self.allm5, c500c, np.mean(self.allz), '500c', '200c')
        M500c, R500c, c500cBis = mass_defs.changeMassDefinition(self.allm5, c500c, np.mean(self.allz), '500c', '500c')
        self.frac_r200c = R200c/R500c
        self.frac_rVIR  = Rvir/R500c

        # integrate to get LX within R500c, R200C, Rvir
        self.alllx       = np.empty(self.nsim2)
        self.alllx_R200c = np.empty(self.nsim2)
        self.alllx_Rvir  = np.empty(self.nsim2)
        for i in range(self.nsim2):
            tprof = self.profiles[i]
            self.alllx[i]         = self.calc_lx(tprof, self.allkt[i], self.allm5[i], self.allz[i], 1.0)
            self.alllx_R200c[i]   = self.calc_lx(tprof, self.allkt[i], self.allm5[i], self.allz[i], self.frac_r200c[i] )
            self.alllx_Rvir [i]   = self.calc_lx(tprof, self.allkt[i], self.allm5[i], self.allz[i], self.frac_rVIR[i] )


        ##if self.nsim2<nsim/100.:
            ##profs=np.exp(np.random.multivariate_normal(self.mean_log,self.covor,size=nsim*10))

            ##allz_i  = profs[:,len(self.mean_log)-3]
            ##allkt_i = profs[:,len(self.mean_log)-1]
            ##allm5_i = profs[:,len(self.mean_log)-2]

            ##profiles_i = profs[:,:len(self.xgrid_ext)]
            ### filters output to cover the relevant redshift range
            ##in_zbin = (allz_i<zmax+0.01)&(allz_i>zmin-0.01)

            ##self.allz     = allz_i    [in_zbin]
            ##self.allkt    = allkt_i   [in_zbin]
            ##self.allm5    = allm5_i   [in_zbin]
            ##self.profiles = profiles_i[in_zbin]
            ##self.nsim2 = len(self.allz)

    def populate_cat_Seppi2022(self):
        """
        Assigns to each catalog entry a vector of quantities drawn by the draw_simulated_values function
            - finds the nearest neighbour in redshift and M500c

        Compared to previous incarnations of the algorithm, it ignores the halo dynamical state (Xoff). Indeed Xoff is not available in the Uchuu light cone.

        """
        Tree_profiles = BallTree(np.transpose([self.allz, np.log10(self.allm5)]))
        DATA = np.transpose([self.CAT['redshift_S'], np.log10(self.CAT['M500c'] * self.b_HS)])
        ids_out = Tree_profiles.query(DATA, k=1, return_distance = False)
        ids = np.hstack((ids_out))

        # LX_out 0.5-2 rest frame erg/s
        # to convert to 0.5-2 observed frame
        LX_out        = self.alllx[ids]
        LX_out_R200c  = self.alllx_R200c[ids]
        LX_out_Rvir   = self.alllx_Rvir [ids]
        KT_OUT        = self.allkt[ids]
        CBP_redshift  = self.allz[ids]
        CBP_M500c     = np.log10(self.allm5)[ids]

        # maximum kT for the k-correction
        KT_OUT[KT_OUT>19.4] = 19.4

        self.CAT['CBP_CLUSTER_LX_soft_RF_R500c'] = np.log10( LX_out       )
        self.CAT['CBP_CLUSTER_LX_soft_RF_R200c'] = np.log10( LX_out_R200c )
        self.CAT['CBP_CLUSTER_LX_soft_RF_Rvir' ] = np.log10( LX_out_Rvir  )
        self.CAT['CBP_redshift' ] = CBP_redshift
        self.CAT['CBP_M500c' ] = CBP_M500c
        self.CAT['CBP_kT' ] = KT_OUT


        # Correction to match Anderson et al. 2015
        # scaling in Anderson 2015
        # logLX_A15 = lambda logMs : 2.6 * logMs + 12.5
        # logLX_A15 = lambda logMs : 2.7 * logMs + 11.35
        logLX_A15 = lambda logMs : 3 * logMs + 7.8

        MS = np.log10(self.CAT['obs_sm'])
        dlogM = 0.1
        MS_min = np.min(MS)-0.01
        MS_max = np.max(MS)
        MS_bins = np.arange( MS_min , MS_max + dlogM, dlogM )
        mean_LX_values = np.zeros_like(MS_bins)
        std_LX_values = np.zeros_like(MS_bins)
        for jj, (m0, m1) in enumerate( zip( MS_bins, MS_bins + dlogM ) ):
            s0 = ( MS >= m0 ) & ( MS < m1 )
            if len(np.log10(LX_out)[s0])>0:
                mean_LX_values[jj] = np.mean(np.log10(LX_out)[s0])
                std_LX_values[jj] = np.std(np.log10(LX_out)[s0])
            elif len(np.log10(LX_out)[s0])==0 and m0>11.0:
                mean_LX_values[jj] = mean_LX_values[jj-1]
                std_LX_values[jj] = std_LX_values[jj-1]
            else : # len(np.log10(LX_out)[s0])==0 :
                mean_LX_values[jj] = logLX_A15(m0 + dlogM/2.)
                std_LX_values[jj] = 1.0

        goal_mean = logLX_A15(MS_bins + dlogM/2.)
        correction = goal_mean - mean_LX_values

        itp = interp1d(
            np.hstack(( MS_bins[0] - dlogM, MS_bins + dlogM/2., MS_bins[-1] +dlogM )),
            np.hstack(( correction[0], correction, correction[-1] ))
            )

        LX_GroupCorrection = 10**itp(MS)

        Corr_LX_out       = LX_GroupCorrection * LX_out
        Corr_LX_out_R200c = LX_GroupCorrection * LX_out_R200c
        Corr_LX_out_Rvir  = LX_GroupCorrection * LX_out_Rvir

        self.CAT['CLUSTER_LX_soft_RF_R500c'] = np.log10( Corr_LX_out       )
        self.CAT['CLUSTER_LX_soft_RF_R200c'] = np.log10( Corr_LX_out_R200c )
        self.CAT['CLUSTER_LX_soft_RF_Rvir' ] = np.log10( Corr_LX_out_Rvir  )

        # attenuation grid should cover down to 0.05 keV
        #itp_logNH, itp_z, itp_kt, itp_frac_obs = np.loadtxt( os.path.join( os.environ['GIT_AGN_MOCK'], "data", "xray_k_correction", "fraction_observed_clusters.txt"), unpack=True )


    def populate_cat_forcing_SR(self):
        """
        Assigns to each catalog entry a vector of quantities drawn by the draw_simulated_values function
            - finds the nearest neighbour in redshift and M500c

        Compared to previous incarnations of the algorithm, it ignores the halo dynamical state (Xoff). Indeed Xoff is not available in the Uchuu light cone.

        """
        Tree_profiles = BallTree(np.transpose([self.allz, np.log10(self.allm5)]))
        DATA = np.transpose([self.CAT['redshift_S'], np.log10(self.CAT['M500c'] * self.b_HS)])
        ids_out = Tree_profiles.query(DATA, k=1, return_distance = False)
        ids = np.hstack((ids_out))

        # LX_out 0.5-2 rest frame erg/s
        # to convert to 0.5-2 observed frame
        LX_out        = self.alllx[ids]
        LX_out_R200c  = self.alllx_R200c[ids]
        LX_out_Rvir   = self.alllx_Rvir [ids]
        KT_OUT        = self.allkt[ids]
        CBP_redshift  = self.allz[ids]
        CBP_M500c     = np.log10(self.allm5)[ids]

        # maximum kT for the k-correction
        KT_OUT[KT_OUT>19.4] = 19.4

        self.CAT['CBP_CLUSTER_LX_soft_RF_R500c'] = np.log10( LX_out       )
        self.CAT['CBP_CLUSTER_LX_soft_RF_R200c'] = np.log10( LX_out_R200c )
        self.CAT['CBP_CLUSTER_LX_soft_RF_Rvir' ] = np.log10( LX_out_Rvir  )
        self.CAT['CBP_redshift' ] = CBP_redshift
        self.CAT['CBP_M500c' ] = CBP_M500c
        self.CAT['CBP_kT' ] = KT_OUT


        # Correction to match log10(LX_05_20/EZ) = 1.5 log10(M500c EZ) + 22
        fun_M_L =  lambda x : 1.5 * x + 22.
        EZ = cosmo.efunc(self.CAT['redshift_S'])
        M500cEZ = self.CAT['CBP_M500c'] + np.log10( EZ)
        LX_oEZ = self.CAT['CBP_CLUSTER_LX_soft_RF_R500c'] - np.log10(EZ)
        dlogM = 0.1
        M500cEZ_min = np.min(M500cEZ)-0.01
        M500cEZ_max = np.max(M500cEZ)
        M500cEZ_bins = np.arange( M500cEZ_min , M500cEZ_max + dlogM, dlogM )
        mean_LX_values = np.zeros_like(M500cEZ_bins)
        std_LX_values = np.zeros_like(M500cEZ_bins)
        for jj, (m0, m1) in enumerate( zip( M500cEZ_bins, M500cEZ_bins + dlogM ) ):
            s0 = ( M500cEZ >= m0 ) & ( M500cEZ < m1 )
            if len(LX_oEZ[s0])>0:
                mean_LX_values[jj] = np.mean(LX_oEZ[s0])
                std_LX_values[jj] = np.std(LX_oEZ[s0])
            elif len(LX_oEZ[s0])==0 and m0>11.0:
                mean_LX_values[jj] = mean_LX_values[jj-1]
                std_LX_values[jj] = std_LX_values[jj-1]
            else : # len(LX_oEZ[s0])==0 :
                mean_LX_values[jj] = fun_M_L(m0 + dlogM/2.)
                std_LX_values[jj] = 1.0

        goal_mean = fun_M_L(M500cEZ_bins + dlogM/2.)
        correction = goal_mean - mean_LX_values

        itp = interp1d(
            np.hstack(( M500cEZ_bins[0] - dlogM, M500cEZ_bins + dlogM/2., M500cEZ_bins[-1] +dlogM )),
            np.hstack(( correction[0], correction, correction[-1] ))
            )

        LX_GroupCorrection = 10**itp(M500cEZ)

        self.CAT['LX_GroupCorrection'] = LX_GroupCorrection
        self.CAT['CLUSTER_LX_soft_RF_R500c'] = np.log10( LX_GroupCorrection * LX_out )
        self.CAT['CLUSTER_LX_soft_RF_R200c'] = np.log10( LX_GroupCorrection * LX_out_R200c )
        self.CAT['CLUSTER_LX_soft_RF_Rvir' ] = np.log10( LX_GroupCorrection * LX_out_Rvir )
        self.CAT['CLUSTER_kT'] = LX_GroupCorrection*self.CAT['CBP_kT' ]

        # Correction to match log10(kT/EZ^(2/3)) = 0.6 log10(M500c) - 8
        fun_M_T =  lambda x : 0.6 * x - 8.
        M500c = self.CAT['CBP_M500c']
        kT_oEZ23 = np.log10( self.CAT['CBP_kT' ] / EZ**(2./3.) )
        dlogM = 0.1
        M500c_min = np.min(M500c)-0.01
        M500c_max = np.max(M500c)
        M500c_bins = np.arange( M500c_min , M500c_max + dlogM, dlogM )
        mean_kT_values = np.zeros_like(M500c_bins)
        std_kT_values = np.zeros_like(M500c_bins)
        for jj, (m0, m1) in enumerate( zip( M500c_bins, M500c_bins + dlogM ) ):
            s0 = ( M500c >= m0 ) & ( M500c < m1 )
            if len(kT_oEZ23[s0])>0:
                mean_kT_values[jj] = np.mean(kT_oEZ23[s0])
                std_kT_values[jj] = np.std(kT_oEZ23[s0])
            elif len(kT_oEZ23[s0])==0 and m0>11.0:
                mean_kT_values[jj] = mean_kT_values[jj-1]
                std_kT_values[jj] = std_kT_values[jj-1]
            else : # len(kT_oEZ23[s0])==0 :
                mean_kT_values[jj] = fun_M_L(m0 + dlogM/2.)
                std_kT_values[jj] = 1.0

        goal_mean = fun_M_T(M500c_bins + dlogM/2.)
        correction = goal_mean - mean_kT_values

        itp = interp1d(
            np.hstack(( M500c_bins[0] - dlogM, M500c_bins + dlogM/2., M500c_bins[-1] +dlogM )),
            np.hstack(( correction[0], correction, correction[-1] ))
            )

        kT_GroupCorrection = 10**itp(M500c)

        self.CAT['kT_GroupCorrection'] = kT_GroupCorrection
        self.CAT['CLUSTER_LX_soft_RF_R500c_kTcorr'] = np.log10( kT_GroupCorrection * LX_out )
        self.CAT['CLUSTER_LX_soft_RF_R200c_kTcorr'] = np.log10( kT_GroupCorrection * LX_out_R200c )
        self.CAT['CLUSTER_LX_soft_RF_Rvir_kTcorr' ] = np.log10( kT_GroupCorrection * LX_out_Rvir )
        self.CAT['CLUSTER_kT_kTcorr'] = kT_GroupCorrection*self.CAT['CBP_kT' ]

        # attenuation grid should cover down to 0.05 keV
        #itp_logNH, itp_z, itp_kt, itp_frac_obs = np.loadtxt( os.path.join( os.environ['GIT_AGN_MOCK'], "data", "xray_k_correction", "fraction_observed_clusters.txt"), unpack=True )


    def populate_cat(self):
        """
        Assigns to each catalog entry a vector of quantities drawn by the draw_simulated_values function
            - finds the nearest neighbour in redshift and M500c

        Compared to previous incarnations of the algorithm, it ignores the halo dynamical state (Xoff). Indeed Xoff is not available in the Uchuu light cone.

        """
        Tree_profiles = BallTree(np.transpose([self.allz, np.log10(self.allm5)]))
        DATA = np.transpose([self.CAT['redshift_S'], np.log10(self.CAT['M500c'] * self.b_HS)])
        ids_out = Tree_profiles.query(DATA, k=1, return_distance = False)
        ids = np.hstack((ids_out))

        # LX_out 0.5-2 rest frame erg/s
        # to convert to 0.5-2 observed frame
        LX_out        = self.alllx[ids]
        LX_out_R200c  = self.alllx_R200c[ids]
        LX_out_Rvir   = self.alllx_Rvir [ids]
        KT_OUT        = self.allkt[ids]
        CBP_redshift  = self.allz[ids]
        CBP_M500c     = np.log10(self.allm5)[ids]

        # maximum kT for the k-correction
        KT_OUT[KT_OUT>19.4] = 19.4

        self.CAT['CBP_CLUSTER_LX_soft_RF_R500c'] = np.log10( LX_out       )
        self.CAT['CBP_CLUSTER_LX_soft_RF_R200c'] = np.log10( LX_out_R200c )
        self.CAT['CBP_CLUSTER_LX_soft_RF_Rvir' ] = np.log10( LX_out_Rvir  )
        self.CAT['CBP_redshift' ] = CBP_redshift
        self.CAT['CBP_M500c' ] = CBP_M500c
        self.CAT['CBP_kT' ] = KT_OUT

        # intrinsic scatter
        #sigma_kT_int = np.random.normal(loc=0.0, scale=0.15, size=len(self.CAT))
        #sigma_LXX_int = 0.25 / 0.15 * sigma_kT_int
        #sigma_var  = np.random.normal(loc=0.0, scale=0.02, size=len(self.CAT))
        #sigma_var2 = np.random.normal(loc=0.0, scale=0.02, size=len(self.CAT))
        # Matching log10(LX_05_20/EZ) = 1.5 log10(M500c EZ) + 22
        fun_M_L =  lambda log10M500c : 1.5 * log10M500c + 22.
        sigma_LX = 0.3
        # Matching log10(kT/EZ^(2/3)) = 0.6 log10(M500c) - 8
        fun_M_T =  lambda log10M500c : 0.6 * log10M500c - 8.
        sigma_kT = 0.2
        # covariance
        rho = 0.9
        cov_kT_LX = np.array([[sigma_kT**2, rho*sigma_kT*sigma_LX ],[rho*sigma_kT*sigma_LX, sigma_LX**2]])
        # generates means
        EZ = cosmo.efunc(self.CAT['redshift_S'])
        self.CAT['kT_Mean_oEzm23'] = ( fun_M_T(np.log10(self.CAT['M500c'])) ) #+ 2./3. * np.log10(EZ)
        self.CAT['LX_Mean_eZm1'] = fun_M_L( np.log10(self.CAT['M500c']) ) #+ 2*np.log10(EZ)
        # generate values with correlated scatter
        corr_scat = np.random.multivariate_normal([0,0], cov_kT_LX, size=len(self.CAT))
        corr_scat_kT = corr_scat.T[0]
        corr_scat_LX = corr_scat.T[1]
        self.CAT['CLUSTER_kT'] = self.CAT['kT_Mean_oEzm23'] + corr_scat_kT + 2./3. * np.log10(EZ)
        self.CAT['CLUSTER_LX_soft_RF_R500c'] = self.CAT['LX_Mean_eZm1'] + corr_scat_LX + 2*np.log10(EZ)

    def make_photons(self):
        """
        creates a list of photons on the sky"""




