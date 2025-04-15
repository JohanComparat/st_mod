"""
Model to paint a dark matter halo light cone with X-ray properties related to the Active Galactic Nuclei

Follows the cluster model of Comparat et al. 2019 used and improved on by Liu et al. 2022 (eFEDs simulations, more accurate K-corrections).

input :
 - z_dir : name of the redshift slice (to retrieve the list of simulated galaxies, glist.fits files). Accessible via the UCHUU environment variable.
 - LC_dir : directory of the light cone
 - str_scatter_0 : scatter between hard LX and stellar mass
 - str_fsat : satellite fraction

output :
 - AGN\_list\_sigma\_'+str_scatter_0+'_fsat_'+str_fsat+'.fits : file containing AGN properties

"""
import time
t0 = time.time()

import os, glob, sys
from scipy.interpolate import interp1d
from scipy.special import erf
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
from scipy.stats import norm
from astropy.table import Table, vstack
import numpy as np
import extinction

class AGN:
    def __init__(self, z_dir, LC_dir = 'FullSky', scatter_0 = 0.8, f_sat = 10 ):
        print('Initiates AGN class model')
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
        self.LC_MetaData = Table.read( os.path.join(os.environ['UCHUU'], self.LC_dir, 'area_per_replica_'+self.z_dir+'.fits') )
        # AGN HAM parameters
        self.scatter_0 = scatter_0
        self.str_scatter_0 = str(self.scatter_0)
        self.str_fsat = str(f_sat)
        self.f_sat = f_sat/100.

        #
        # Duty cycle
        #
        #f_duty = interp1d(np.array([0., 0.75, 2., 3.5, 10.1]), np.array([0.1, 0.2, 0.3, 0.3, 0.3]))
        f_duty_1 = interp1d(np.array([0.   , 0.25,  0.75, 1.75, 10.1]) ,
                        10**np.array([-1.0, -0.90, -0.75, -0.6, -0.40]) )
        S_trans  = 1.5
        LX_LIM_DC = 40.
        DC_MAX = f_duty_1(self.mean_z)
        x_M = np.arange(6, 13.5, 0.01)
        # add exponential cut-ff below 1e9.5
        #def f_DC_mod(MS):
            #return DC_MAX * (0.5 + 0.5 * erf((MS - (LX_LIM_DC-(30.4+self.mean_z/5.))) / S_trans)) #* np.e**(2*(MS-9))
        def f_DC_mod(MS):
            return DC_MAX *np.ones_like(MS) # * 2 * ( 0.5 + 0.5 * erf((MS - 9.25) / 0.15) )
            #return DC_MAX * (0.5 + 0.5 * erf((MS - (LX_LIM_DC-(30.4+self.mean_z/5.))) / S_trans)) * ( 0.5 + 0.5 * erf((MS - 9.25) / 0.15) )
            #return DC_MAX * (0.5 + 0.5 * erf((MS - (LX_LIM_DC-(30.4+self.mean_z/5.))) / S_trans)) * ( 0.5 + 0.5 * erf((MS - 9.5) / 0.05) )
        self.f_duty = interp1d( x_M, f_DC_mod(x_M) )
        #f_duty_E = interp1d( x_M, f_DC_mod_E(x_M) )
        #print(f_duty([8, 8.5, 9, 9.5, 10]))
        print('duty cycle MAX', DC_MAX)
        print('duty cycle vs Ms: 8,12.5 by 0.5', self.f_duty([8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5]))


    def get_tabulated_AGN(self, z_bins_val, dCx_min, dCx_max, dCy_min, dCy_max, dCz_min, dCz_max):
        """
        Concatenates all tabulated AGN in the redshift shell of interest in a single table, saved as attribute self.AGN as well as the total number of tabulated AGN self.N_AGN_tabulated
        """
        # z bins in which AGN are tabulated
        #DZ=0.01
        #z_bins = np.arange(0.0, 6., DZ)
        #z_ups = z_bins+DZ
        #jj = np.argmin( abs(z_bins - z_mean))
        #print(jj)
        #j0 = np.searchsorted(z_bins, self.z_min, side='right')-1
        #j1 = np.searchsorted(z_bins, self.z_max, side='right')
        #print(self.z_min, self.z_max)
        #AGN_ALL = []
        #for jj in np.arange(j0, j1+1, 1):
        p_2_AGN = os.path.join( os.environ["UCHUU"], 'AGN_LX_tables', 'LX_table_' +str(np.round(z_bins_val,2))+ '.fits')
        AGN = Table.read(p_2_AGN)
        AGN.remove_columns(['FX_hard', 'logNH', 'agn_type' ])
        p_2_xyz = os.path.join( os.environ["UCHUU"], 'AGN_LX_tables', 'XYZ_position_' +str(np.round(z_bins_val,2))+ '.fits' )
        XYZ = Table.read(p_2_xyz)
        N_AGN_tabulated = len(AGN)
        # only works for the first replication d_comoving 2 Gpc/h !
        dCx_sel = (XYZ['X']>=dCx_min) & (XYZ['X'] < dCx_max)
        dCy_sel = (XYZ['Y']>=dCy_min) & (XYZ['Y'] < dCy_max)
        dCz_sel = (XYZ['Z']>=dCz_min) & (XYZ['Z'] < dCz_max)
        z_sel=(AGN['z']>=self.z_min)&(AGN['z']<=self.z_max)&(dCx_sel)&(dCy_sel)&(dCz_sel)
        #& (np.random.random(size=N_AGN_tabulated)<self.sky_frac)
        #AGN_ALL.append( AGN[z_sel] )
        AGN = AGN[z_sel]
        #print(p_2_AGN, len(AGN), len(AGN[z_sel]))

        #AGN = vstack(( AGN_ALL ))
        self.AGN_tab = AGN[np.argsort(AGN['LX_hard'])[::-1]]
        self.N_AGN_tabulated = len(AGN)
        print(self.N_AGN_tabulated, 'AGN tabulated')
        #AGN_ALL = 0
        AGN = 0
        XYZ = 0

    def get_z_mass_id(self, t_sim_gal):
        """
        Downsamples to the f_sat parameter and applies the duty cycle

        Saves the following attributes :
            - self.IDS the indexes of the galaxies with an AGN
            - self.N_active : number of AGN
            - self.AGN : glist catalog restricted to the active IDS
        """
        ZZ = t_sim_gal['redshift_S']
        MM = np.log10(t_sim_gal['obs_sm'])
        HALO_pid = t_sim_gal['upid']
        cen = (HALO_pid==-1)
        sat = (cen==False)
        N_galaxies = len(HALO_pid)
        N_galaxies_cen = len(HALO_pid[cen])
        N_galaxies_sat = len(HALO_pid[sat])
        native_f_sat = N_galaxies_sat*1./N_galaxies_cen
        print('native N, cen, sat, f_sat', N_galaxies, N_galaxies_cen, N_galaxies_sat, native_f_sat)
        rds = np.random.random(N_galaxies)
        if self.f_sat==0:
            keep = (cen)
            N_kept1 = len(keep.nonzero()[0])
        else:
            if self.f_sat > native_f_sat:
                # downsample centrals
                print('downsamples centrals')
                N_cen_goal = N_galaxies_sat / self.f_sat
                print('N cen goal', N_cen_goal)
                sel_cen = (rds < N_cen_goal / N_galaxies_cen)
                all_cen = (cen)&(sel_cen)
                keep = (sat)|(all_cen)
                N_kept1 = len(keep.nonzero()[0])
                #print(N_kept1)
            if self.f_sat <= native_f_sat:
                # downsample sat
                print('downsamples sat')
                N_sat_goal = N_galaxies_cen * self.f_sat
                print('N sat goal', N_sat_goal)
                sel_sat = (rds<N_sat_goal/N_galaxies_sat)
                all_sat = (sat)&(sel_sat)
                keep = (cen)|(all_sat)
                N_kept1 = len(keep.nonzero()[0])
                #print(N_kept1)
        # APPLY DUTY CYCLE
        print('renormalizing DC by', N_galaxies/N_kept1)
        f_duty_realization = self.f_duty(MM)*N_galaxies/N_kept1
        active = (np.random.random(size=len(ZZ)) <= f_duty_realization)
        z_range = (ZZ >= self.z_min) & (ZZ < self.z_max) & (t_sim_gal['nH']>0) & (t_sim_gal['ebv']>=0)
        self.IDS = np.arange(len(ZZ))    [active & keep & z_range]
        self.AGN = t_sim_gal[self.IDS]
        self.N_active = len(self.AGN)
        print('N AGN', self.N_active)

    def abundance_matching(self):
        """
        Matches LX hard and stellar mass after applying the scatter

        Adds the log_{10}(L_X^{2-10 keV}) 'LX_hard' and the realization of the scatter parameter (scatter_LX_Ms) to the self.AGN table

        """
        # match the exact number of tabulated AGN and active galaxies
        if self.N_AGN_tabulated < self.N_active  :
            # too many active galaxies (box resolution is higher than AGN tabulated
            # just populate the higher stellar mass ones in the HAM
            mass_sort_idx = np.argsort(np.log10(self.AGN['obs_sm']))
            id_selection = mass_sort_idx [ - self.N_AGN_tabulated : ]
            self.AGN = self.AGN[id_selection]
            self.N_active = len(self.AGN)
        #
        #
        # HAM procedure
        #
        #
        self.AGN['scatter_LX_Ms'] = norm.rvs(loc=0, scale=self.scatter_0, size=self.N_active )
        M_scatt = np.log10(self.AGN['obs_sm']) + self.AGN['scatter_LX_Ms']
        ids_M_scatt = np.argsort(M_scatt)[::-1]
        # output numbers
        lx = np.zeros_like(M_scatt)
        #print(self.AGN_tab['LX_hard'])
        print( np.transpose([M_scatt[ids_M_scatt], self.AGN_tab['LX_hard'][:len(lx)]]) )
        lx[ids_M_scatt] = self.AGN_tab['LX_hard'][:len(lx)]
        self.AGN['LX_hard'] = lx
        #lsar = np.zeros_like(lx)
        #lsar[ids_M_scatt] = self.AGN_tab['LX_hard'] - np.log10(self.AGN['obs_sm'])[ids_M_scatt]

    def get_obscured_fractions(self):
        """
        Computes the obscured fraction model (equations 4-11, 12-15 of Comparat et al. 2019)

        Adds $log_{10}(n_H)$ values ('logNH') to the self.AGN table
        """

        # too many CTK at high luminosity
        # Eq. 4
        #def f_thick(LXhard, z): return 0.30
        def thick_LL(z, lx0 = 41.5): return lx0 + np.arctan(z*5)*1.5
        def f_thick(LXhard, z): return 0.30 * (0.5 + 0.5 * erf((thick_LL(z) - LXhard) / 0.25))

        # too many absorbed ones
        # Eq. 7
        def f_2(LXhard, z): return 0.9 * (41 / LXhard)**0.5

        # fiducial
        # Eq. 8
        def f_1(LXhard, z): return f_thick(LXhard, z) + 0.01 + erf(z / 4.) * 0.3

        # Eq. 10
        def LL(z, lx0 = 43.2): return lx0 + erf(z) * 1.2

        # Eq. 5,6
        def fraction_ricci(LXhard, z, width = 0.6): return f_1(LXhard,z) + (f_2(LXhard, z) - f_1(LXhard,z)) * (0.5 + 0.5 * erf((LL(z) - LXhard) / width))

        # initializes logNH
        logNH = np.zeros(self.N_active)

        # obscuration, after the equations above
        randomNH = np.random.rand(self.N_active)

        # unobscured 20-22
        frac_thin = fraction_ricci(self.AGN['LX_hard'], self.AGN['redshift_S'])
        thinest = (randomNH >= frac_thin)

        # thick obscuration, 24-26
        thick = (randomNH < f_thick(self.AGN['LX_hard'], self.AGN['redshift_S']))

        # obscured 22-24
        obscured = (thinest == False) & (thick == False)

        # assigns logNH values randomly :
        logNH[thick]    = np.random.uniform(24, 26, len(logNH[thick]))
        logNH[obscured] = np.random.uniform(22, 24, len(logNH[obscured]))
        logNH[thinest]  = np.random.uniform(20, 22, len(logNH[thinest]))
        self.AGN['logNH'] = logNH


    def compute_fluxes(self):
        """
        Computes the fluxes in observed frame (K-corrections with obscuration in the 0.5-2 and 2-10 keV bands)

        Adds FX_hard, LX_soft, FX_soft, FX_soft_attenuated to the self.AGN table
        """
        NHS = np.arange(20, 26 + 0.05, 0.4)
        path_2_RF_obs_hard = os.path.join(
            os.environ['GIT_STMOD_DATA'],
            "data", "models", "model_AGN",
            "xray_k_correction",
            "v3_fraction_observed_A15_RF_hard_Obs_hard_fscat_002.txt")

        obscuration_z_grid, obscuration_nh_grid, obscuration_fraction_obs_erosita = np.loadtxt(
            path_2_RF_obs_hard, unpack=True)

        xy = np.c_[obscuration_z_grid, obscuration_nh_grid]
        obscuration_itp_H_H = LinearNDInterpolator(xy, obscuration_fraction_obs_erosita)

        percent_observed_itp = interp1d(
            np.hstack((20 - 0.1, NHS, 26 + 0.1)),
            np.hstack((
                obscuration_itp_H_H(self.z_mean, 20.),
                np.array([obscuration_itp_H_H(z_i, logNH_i) for z_i, logNH_i in zip(self.z_mean * np.ones_like(NHS), NHS)]),
                obscuration_itp_H_H(self.z_mean, 26.))))

        percent_observed_H_H = percent_observed_itp(self.AGN['logNH'])
        lx_obs_frame_2_10_attenuated = self.AGN['LX_hard'] + np.log10(percent_observed_H_H)

        dl_cm = self.AGN['dL']
        #print(self.AGN['dL'])
        self.AGN['FX_hard'] = lx_obs_frame_2_10_attenuated - np.log10 (4 * np.pi ) - 2. * np.log10(dl_cm) # / h**3

        # obs X-ray 2-10 keV ==>> obs 0.5-2
        # v3_fraction_observed_A15_RF_hard_Obs_soft_fscat_
        # path_2_hard_RF_obs_soft
        # link to X-ray K-correction and attenuation curves
        path_2_hard_RF_obs_soft = os.path.join(
            os.environ['GIT_STMOD_DATA'],
            "data", "models", "model_AGN",
            "xray_k_correction",
            "v3_fraction_observed_A15_RF_hard_Obs_soft_fscat_002.txt")

        obscuration_z_grid, obscuration_nh_grid, obscuration_fraction_obs_erosita = np.loadtxt(path_2_hard_RF_obs_soft, unpack=True)

        xy = np.c_[obscuration_z_grid, obscuration_nh_grid]
        obscuration_itp_H_S = LinearNDInterpolator(xy, obscuration_fraction_obs_erosita)

        percent_observed_itp = interp1d(
            np.hstack((20 - 0.1, NHS, 26 + 0.1)),
            np.hstack((
                obscuration_itp_H_S(self.z_mean, 20.),
                np.array([obscuration_itp_H_S(z_i, logNH_i) for z_i, logNH_i in zip(self.z_mean * np.ones_like(NHS), NHS)]),
                obscuration_itp_H_S(self.z_mean, 26.))))

        percent_observed_H_S = percent_observed_itp(self.AGN['logNH'])

        path_2_hard_RF_obs_soft_3D = np.array( glob.glob( os.path.join(
            os.environ['GIT_STMOD_DATA'],
            "data", "models", "model_AGN",
            "xray_k_correction",
            "v3_fraction_observed_A15_RF_hard_Obs_soft_fscat_002_GALnH_*.txt") ) )
        path_2_hard_RF_obs_soft_3D.sort()
        # create the 3D interpolation
        grid_z, grid_nh, grid_galnh, transmission = [], [], [], []
        for el in path_2_hard_RF_obs_soft_3D:
            grid_z_i, grid_nh_i, grid_galnh_i, transmission_i = np.loadtxt(el, unpack=True)
            grid_z      .append(grid_z_i      )
            grid_nh     .append(grid_nh_i     )
            grid_galnh  .append(grid_galnh_i  )
            transmission.append(transmission_i)

        points = np.transpose([
            np.hstack( grid_z ),
            np.hstack( grid_nh ),
            np.hstack( grid_galnh )])

        values = np.hstack( transmission)

        ITP3d = NearestNDInterpolator(points, values)
        # redshift, AGN nH, galactic nH
        percent_observed_H_S_galNH = ITP3d(self.AGN['redshift_S'], self.AGN['logNH'], np.log10(self.AGN['nH']))

        self.AGN['LX_soft'] = self.AGN['LX_hard'] + np.log10(percent_observed_H_S)
        self.AGN['FX_soft'] = self.AGN['LX_soft'] - np.log10(4 * np.pi) - 2*np.log10(dl_cm)
        self.AGN['LX_soft_MWattenuated'] = self.AGN['LX_hard'] + np.log10(percent_observed_H_S_galNH)
        self.AGN['FX_soft_MWattenuated'] = self.AGN['LX_soft'] - np.log10(4 * np.pi) - 2*np.log10(dl_cm)
        print(self.AGN)

    def compute_agn_type(self, z, lx, logNH):
        """
        Assigns a type to an AGN population

        parameters:
        - z: redshift
        - lx: hard X-ray luminosity (log10)
        - logNH: nH value (log10)

        return: array of AGN types
        """

        # Adds type 11, 12, 21, 22
        # Follows Merloni et al. 2014
        # equation 16 of Comparat et al. 2019

        def fraction_22p21_merloni(lx): return (
            0.5 + 0.5 * erf((-lx + 44.) / 0.9)) * 0.69 + 0.26

        # LF in the mock, starting parameters
        dlogf = 0.05
        Lbin_min = 36
        fbins = np.arange(Lbin_min, 48, dlogf)
        xf = fbins[:-1] + dlogf / 2.

        # boundary between the 22 and the 21 populations
        limit = fraction_22p21_merloni((fbins[1:] + fbins[:-1]) * 0.5)
        # selection per obscuration intensity
        nh_21 = (logNH <= 22.)
        nh_23 = (logNH > 22.)  # &(logNH<=26.)
        # initiate columns to compute
        opt_type = np.zeros(self.N_active).astype('int')
        rd = np.random.rand(self.N_active)
        # compute histograms of LX for different obscurations
        nall = np.histogram(lx, fbins)[0]       # all
        nth = np.histogram(lx[nh_23], fbins)[0]  # thin
        nun = np.histogram(lx[nh_21], fbins)[0]  # unobscured
        fr_thk = nth * 1. / nall  # fraction of obscured
        fr_un = nun * 1. / nall  # fraction of unobscured
        # first get the type 12: NH absorption but optically unobscured
        # to be chosen in obscured population
        n_per_bin_12 = (fr_thk - limit) * nall
        sel_12 = (np.ones(len(z)) == 0)
        for bin_low, bin_high, num_needed, nn_un in zip(
                fbins[:-1], fbins[1:], n_per_bin_12.astype('int'), nth):
            if num_needed > 0 and nn_un > 0:
                frac_needed = num_needed * 1. / nn_un
                sel_12 = (sel_12) | (
                    (lx > bin_low) & (
                        lx < bin_high) & (nh_23) & (
                        rd < frac_needed))
        t_12 = (nh_23) & (sel_12)
        # second the types 21
        # to be chosen in nun
        n_per_bin_21 = (-fr_thk + limit) * nall
        sel_21 = (np.ones(len(z)) == 0)
        for bin_low, bin_high, num_needed, nn_un in zip(
                fbins[:-1], fbins[1:], n_per_bin_21.astype('int'), nun):
            if num_needed > 0 and nn_un > 0:
                frac_needed = num_needed * 1. / nn_un
                sel_21 = (sel_21) | (
                    (lx > bin_low) & (
                        lx < bin_high) & (nh_21) & (
                        rd < frac_needed))
        t_21 = (nh_21) & (sel_21)
        # finally the types 11 and 22
        t_11 = (nh_21) & (t_21 == False)
        t_22 = (nh_23) & (t_12 == False)
        opt_type[t_22] = 22
        opt_type[t_12] = 12
        opt_type[t_11] = 11
        opt_type[t_21] = 21
        return opt_type

    def compute_r_mag(self):
        """
        Computes the observed r-band magnitude from X-ray using the alpha OX relation.
        Adds : self.AGN['SDSS_r_AB'] and self.AGN['SDSS_r_AB_attenuated']
        """

        def r_mean(log_FX0520): return -2. * log_FX0520 - 7.

        def scatter_t1(NN): return norm.rvs(loc=0.0, scale=1.0, size=NN)

        self.AGN['SDSS_r_AB'] = r_mean(self.AGN['FX_soft']) + scatter_t1(int(self.N_active))
        # optical extinction, Fitzpatrick 99
        ebv_values = np.hstack((np.arange(0., 5., 0.01), 10**np.arange(1, 4, 0.1)))
        ext_values = np.array([extinction.fitzpatrick99(
            np.array([6231.]), 3.1 * EBV, r_v=3.1, unit='aa')[0] for EBV in ebv_values])
        ext_interp = interp1d(ebv_values, ext_values)
        self.AGN['SDSS_r_AB_attenuated'] = self.AGN['SDSS_r_AB'] + ext_interp(self.AGN['ebv'])

#class SMBH_SAR_model:
#    def __init__(self, z_dir, LC_dir = 'FullSky', scatter_0 = 0.8, f_sat = 10 ):
