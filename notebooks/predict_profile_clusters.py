
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as plt
import os, sys
import glob
import numpy as np
from astropy.table import Table, vstack
from scipy.interpolate import interp1d
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoUNIT = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
h = 0.6774
L_box = 2000.0 / h
cosmo = cosmoUNIT
nl = lambda sel : len(sel.nonzero()[0])

from scipy.interpolate import interp1d
from scipy.integrate import quad

MD40 = Table.read('/home/comparat/sf_Shared/data/MultiDark/MD_4.0Gpc/MD40_eRO_CLU_b8_CM_0_pixS_10.0_M500c_13.0_FX_-14.5_25Apr2022_Profiles.fits')

fig_dir       = os.path.join(os.environ['GIT_STMOD'], 'data', 'validation','validation_GAS', 'Xray_profiles')

z_mins = np.array ([0.0, 0.1,  0.2     ])
z_maxs = np.array ([0.1, 0.2,  0.3     ])
z_strs = np.array (['z0001', 'z0102',  'z0203' ])

for z_min, z_max, z_str in zip(z_mins, z_maxs, z_strs):

    z_MD40 = (MD40['redshift_R']>z_min) & (MD40['redshift_R']<z_max)
    SM_arr_MD40 = MD40['galaxy_SMHMR_mass']
    prof = MD40['profiles']
    FX = MD40['CLUSTER_FX_soft']
    r500c = MD40['R500c_kpc']

    path_2_cbp= os.path.join(os.environ['GIT_STMOD'], 'data/models/model_GAS')
    covor     = np.loadtxt(os.path.join(path_2_cbp, 'covmat_xxl_hiflugcs_xcop.txt'))
    xgrid_ext = np.loadtxt(os.path.join(path_2_cbp, 'radial_binning.txt'))
    mean_log  = np.loadtxt(os.path.join(path_2_cbp, 'mean_pars.txt'))
    coolfunc  = np.loadtxt(os.path.join(path_2_cbp, 'coolfunc.dat'))
    xgrid_ext_0 = np.hstack((0., xgrid_ext[:-1] ))
    xgrid_ext_1 = xgrid_ext
    dx = np.empty(len(xgrid_ext))
    dx[0]=xgrid_ext[0]
    dx[1:len(xgrid_ext)]=(np.roll(xgrid_ext,-1)-xgrid_ext)[:len(xgrid_ext)-1]

    ##def calc_lx_w_Mgas(prof,r500,fraction=1):
        ##"""
        ##Compute the X-ray luminosity in the profile
        ##to be extended to 3x r500c

        ##Compute r_{500c} :
            ##r_{500c} = \left(\frac{3 M_{500c}}{ 4. \pi 500 \rho_c(z)  }\right)^{1/3} [ cm ].

        ##profile_emission = profile x rescale_factor
        ##rescale_factor = \sqrt(kT/10.0) E^3(z)
        ##CF(kT) = cooling function, show the curve
        ##L_X(r) = \Sigma_{<r}( profile_emission r_{500c}^2 2 \pi x CF(kT) Mpc=3.0856776e+24 dx )
        ##L_{500c} = L_X(1)

        ##"""

    def calc_lx(prof,kt,m5,z,fraction=1):
        """
        Compute the X-ray luminosity in the profile
        to be extended to 3x r500c

        Compute r_{500c} :
            r_{500c} = \left(\frac{3 M_{500c}}{ 4. \pi 500 \rho_c(z)  }\right)^{1/3} [ cm ].

        profile_emission = profile x rescale_factor
        rescale_factor = \sqrt(kT/10.0) E^3(z)
        CF(kT) = cooling function, show the curve
        L_X(r) = \Sigma_{<r}( profile_emission r_{500c}^2 2 \pi x CF(kT) Mpc=3.0856776e+24 dx )
        L_{500c} = L_X(1)

        """
        ez2 = cosmo.efunc(z)**2
        rhoc = cosmo.critical_density(z).value
        r500 = n.power(m5*msun/4.*3./n.pi/500./rhoc,1./3.)
        resfact = n.sqrt(kt/10.0)*n.power(ez2,3./2.)
        prof_em = prof * resfact # emission integral
        tlambda = n.interp(kt,coolfunc[:,0],coolfunc[:,1]) # cooling function
        #print(prof_em*xgrid_ext*r500**2*2.*n.pi*tlambda*Mpc*dx)
        lxcum = n.cumsum(prof_em*xgrid_ext*r500**2*2.*n.pi*tlambda*Mpc*dx) # riemann integral
        lx_500=n.interp(fraction,xgrid_ext,lxcum) # evaluated at R500
        return lx_500


    #profs=n.exp(n.random.multivariate_normal(mean_log,covor,size=nsim))

    #allz_i  = profs[:,len(mean_log)-3]
    #allkt_i = profs[:,len(mean_log)-1]
    #allm5_i = profs[:,len(mean_log)-2]
    #profiles = profs[:,:len(xgrid_ext)]

    #ez2 = cosmo.efunc(CBP_redshift)**2
    #tlambda = n.interp(KT_OUT,coolfunc[:,0],coolfunc[:,1]) # cooling function
    #rescale_F = n.sqrt(KT_OUT/10.0)*n.power(ez2,3./2.) # emission integral
    #profiles_out = n.multiply( profiles[ids].T,  tlambda * rescale_F * GroupCorrection).T

    #dx = np.empty(len(xgrid_ext))
    #dx[0]=xgrid_ext[0]
    #dx[1:len(xgrid_ext)]=(np.roll(xgrid_ext,-1)-xgrid_ext)[:len(xgrid_ext)-1]
    conversion = 3.0856776e+24

    def get_lum(profile, r500crit):
        fraction = 1
        lxcum = np.cumsum(profile * (xgrid_ext * r500crit*u.kpc.to(u.cm))**2 * 2.*np.pi * (conversion * dx)) # riemann integral
        lx_500=np.interp(fraction,xgrid_ext,lxcum) # evaluated at R500
        return lx_500

    #jj=0
    #profile = prof[jj]
    #r500crit = r500c[jj]
    def get_SX_from3D(profile, r500crit):
        """
        integrate
        SX(rp) = int_ rp ^ inf  [ profile(r) (r**2-rp**2)**-0.5 r ] dr
        """
        x_mid = ( xgrid_ext_0*r500crit + xgrid_ext_1*r500crit ) / 2.
        prof_itp = interp1d(x_mid, profile)
        to_int = lambda r, rp : prof_itp(r) * (r**2-rp**2)**(-0.5) * r
        rps = xgrid_ext[:-3]*r500crit
        prof_2d = []
        for rp in rps:
            prof_2d.append(quad(to_int, rp, r500crit*2, args=(rp))[0])
        return rps, 2*np.array(prof_2d)


    #def get_lum_2(profile, r500crit):
        #fraction = 1
        #dx = np.empty(len(xgrid_ext))
        #dx[0]=xgrid_ext[0]
        #dx[1:len(xgrid_ext)]=(np.roll(xgrid_ext,-1)-xgrid_ext)[:len(xgrid_ext)-1]
        #Mpc = 3.0856776e+24
        #lxcum = np.cumsum(profile * xgrid_ext * ( r500crit*u.kpc.to(u.cm))**2 * 2.*np.pi * Mpc * dx) # riemann integral
        #lx_500=np.interp(fraction,xgrid_ext,lxcum) # evaluated at R500
        #return lx_500

    #print(get_lum_2(prof[jj], r500c[jj]), get_lum(prof[jj], r500c[jj]))

    t_out = Table()
    t_out['x_min_R500c'] = np.hstack(( 0., xgrid_ext[:-1] ))
    t_out['x_max_R500c'] = xgrid_ext

    for M0 in np.arange(13, 15.5, 0.1):
        M1 = M0+0.1
        selection = (z_MD40) & (MD40['HALO_M500c']>=10**M0) & (MD40['HALO_M500c']<10**M1)
        #print(M0, nl(selection))
        M_str = 'M500c_'+str(int(np.round(M0*10)))+'_'
        r500crit = np.mean(r500c[selection])
        rps, sx2 = get_SX_from3D(np.mean(prof[selection] , axis=0), r500crit)
        Eq11_2dProj_mean = np.hstack((sx2, 0., 0., 0. ))* (u.kpc.to(u.cm) ) **3 #*u.kpc.to(u.cm)/1000
        x_mid = ( xgrid_ext_0*r500crit + xgrid_ext_1*r500crit ) / 2.
        LUM_Eq11_2dProj_mean = np.cumsum( Eq11_2dProj_mean * 2 * np.pi * x_mid **2 )[15]
        RS_fact = 10**np.mean(MD40['CLUSTER_LX_soft_RF'][selection]) / LUM_Eq11_2dProj_mean
        print('='*100)
        print(M0, nl(selection), 10**np.mean(MD40['CLUSTER_LX_soft_RF'][selection]), LUM_Eq11_2dProj_mean, )
        print('='*100)
        rps, sx2 = get_SX_from3D(RS_fact*np.mean(prof[selection] , axis=0), r500crit)
        Eq11_2dProj_mean = np.hstack((sx2, 0., 0., 0. ))* (u.kpc.to(u.cm) ) **3 #*u.kpc.to(u.cm)/1000
        x_mid = ( xgrid_ext_0*r500crit + xgrid_ext_1*r500crit ) / 2.
        LUM_Eq11_2dProj_mean = np.cumsum( Eq11_2dProj_mean * 2 * np.pi * x_mid **2 )[15]
        print('='*100)
        print(M0, nl(selection), 10**np.mean(MD40['CLUSTER_LX_soft_RF'][selection]), LUM_Eq11_2dProj_mean, )
        print('='*100)
        #RS_fact = 10**np.mean(MD40['CLUSTER_LX_soft_RF'][selection])/get_lum(np.mean(prof[selection] , axis=0), np.mean(r500c[selection]))
        #pr_out = get_SX(np.mean(prof[selection] , axis=0), np.mean(r500c[selection])) * u.kpc.to(u.cm) * u.kpc.to(u.cm) # get_profiles( selection )
        t_out[M_str + 'x_min'] = np.hstack(( 0., xgrid_ext[:-1] )) * np.mean(r500c[selection]) * u.kpc #.to(u.cm) * u.cm
        t_out[M_str + 'x_max'] = xgrid_ext * np.mean(r500c[selection])* u.kpc #.to(u.cm) * u.cm
        t_out[M_str + 'Eq11_mean']  = np.mean( RS_fact*prof[selection] , axis=0)*(u.erg/u.cm**2/u.s/u.Mpc)
        #t_out[M_str + 'Eq11_2dProj_mean']  = Eq11_2dProj_mean
        t_out[M_str + 'SB_X_mean']  = Eq11_2dProj_mean*(u.erg/u.kpc**2/u.s)
        t_out[M_str + 'SB_X_std_fraction'] = np.std( RS_fact*prof[selection] , axis=0)/np.mean( RS_fact*prof[selection] , axis=0)
        #t_out[M_str + 'LUM_Eq11_2dProj_mean'] =
        #t_out[M_str + 'LUM_Eq11_2dProj_mean']  = np.mean(pr_out , axis=0)*(u.erg/u.kpc**2/u.s)
        #t_out[M_str + 'SB_X_std']   = np.std (pr_out , axis=0) *(u.erg/u.kpc**2/u.s)
        t_out[M_str + 'Nhalo']   = nl(selection)
        t_out[M_str + 'MeanStellarMass'] = np.mean(10**MD40['galaxy_SMHMR_mass'][selection])*u.Msun
        t_out[M_str + 'StdStellarMass'] = np.std(10**MD40['galaxy_SMHMR_mass'][selection])*u.Msun
        t_out[M_str + 'MeanRedshift'] = np.mean(MD40['redshift_R'][selection])
        t_out[M_str + 'MeanSoftLuminosity'] = get_lum(RS_fact *np.mean(prof[selection] , axis=0), np.mean(r500c[selection]))*(u.erg/u.s)
        t_out[M_str + 'MeanHaloM500c'] = np.mean(MD40['HALO_M500c'][selection])*(u.Msun)
        t_out[M_str + 'StdHaloM500c'] = np.std(MD40['HALO_M500c'][selection])*(u.Msun)
        #print('='*100)
        #print(M0, nl(selection), 10**np.mean(MD40['CLUSTER_LX_soft_RF'][selection]), t_out[M_str + 'LUM_Eq11_2dProj_mean'][0], t_out[M_str + 'LUM_Eq11_2dProj_mean'][0]/10**np.mean(MD40['CLUSTER_LX_soft_RF'][selection]))
        #print('='*100)

    t_out.write(os.path.join(fig_dir, 'SimulatedProfiles_'+z_str+'_haloMass_QN.fits'), overwrite = True)


    t_out = Table()
    t_out['x_min_R500c'] = np.hstack(( 0., xgrid_ext[:-1] ))
    t_out['x_max_R500c'] = xgrid_ext

    for M0 in np.arange(10.5, 12.2, 0.1):
        M1 = M0+0.1
        selection = (z_MD40) & (SM_arr_MD40>=M0) & (SM_arr_MD40<M1)
        #print(M0, nl(selection))
        M_str = 'M_'+str(int(np.round(M0*10)))+'_'
        r500crit = np.mean(r500c[selection])
        rps, sx2 = get_SX_from3D(np.mean(prof[selection] , axis=0), r500crit)
        Eq11_2dProj_mean = np.hstack((sx2, 0., 0., 0. ))* (u.kpc.to(u.cm) ) **3 #*u.kpc.to(u.cm)/1000
        x_mid = ( xgrid_ext_0*r500crit + xgrid_ext_1*r500crit ) / 2.
        LUM_Eq11_2dProj_mean = np.cumsum( Eq11_2dProj_mean * 2 * np.pi * x_mid **2 )[15]
        RS_fact = 10**np.mean(MD40['CLUSTER_LX_soft_RF'][selection]) / LUM_Eq11_2dProj_mean
        print('='*100)
        print(M0, nl(selection), 10**np.mean(MD40['CLUSTER_LX_soft_RF'][selection]), LUM_Eq11_2dProj_mean, )
        print('='*100)
        rps, sx2 = get_SX_from3D(RS_fact*np.mean(prof[selection] , axis=0), r500crit)
        Eq11_2dProj_mean = np.hstack((sx2, 0., 0., 0. ))* (u.kpc.to(u.cm) ) **3 #*u.kpc.to(u.cm)/1000
        x_mid = ( xgrid_ext_0*r500crit + xgrid_ext_1*r500crit ) / 2.
        LUM_Eq11_2dProj_mean = np.cumsum( Eq11_2dProj_mean * 2 * np.pi * x_mid **2 )[15]
        print('='*100)
        print(M0, nl(selection), 10**np.mean(MD40['CLUSTER_LX_soft_RF'][selection]), LUM_Eq11_2dProj_mean, )
        print('='*100)
        #RS_fact = 10**np.mean(MD40['CLUSTER_LX_soft_RF'][selection])/get_lum(np.mean(prof[selection] , axis=0), np.mean(r500c[selection]))
        #pr_out = get_SX(np.mean(prof[selection] , axis=0), np.mean(r500c[selection])) * u.kpc.to(u.cm) * u.kpc.to(u.cm) # get_profiles( selection )
        t_out[M_str + 'x_min'] = np.hstack(( 0., xgrid_ext[:-1] )) * np.mean(r500c[selection]) * u.kpc #.to(u.cm) * u.cm
        t_out[M_str + 'x_max'] = xgrid_ext * np.mean(r500c[selection])* u.kpc #.to(u.cm) * u.cm
        t_out[M_str + 'Eq11_mean']  = np.mean( RS_fact*prof[selection] , axis=0)*(u.erg/u.cm**2/u.s/u.Mpc)
        #t_out[M_str + 'Eq11_2dProj_mean']  = Eq11_2dProj_mean
        t_out[M_str + 'SB_X_mean']  = Eq11_2dProj_mean*(u.erg/u.kpc**2/u.s)
        t_out[M_str + 'SB_X_std_fraction'] = np.std( RS_fact*prof[selection] , axis=0)/np.mean( RS_fact*prof[selection] , axis=0)
        #t_out[M_str + 'LUM_Eq11_2dProj_mean'] =
        #t_out[M_str + 'LUM_Eq11_2dProj_mean']  = np.mean(pr_out , axis=0)*(u.erg/u.kpc**2/u.s)
        #t_out[M_str + 'SB_X_std']   = np.std (pr_out , axis=0) *(u.erg/u.kpc**2/u.s)
        t_out[M_str + 'Nhalo']   = nl(selection)
        t_out[M_str + 'MeanStellarMass'] = np.mean(10**MD40['galaxy_SMHMR_mass'][selection])*u.Msun
        t_out[M_str + 'StdStellarMass'] = np.std(10**MD40['galaxy_SMHMR_mass'][selection])*u.Msun
        t_out[M_str + 'MeanRedshift'] = np.mean(MD40['redshift_R'][selection])
        t_out[M_str + 'MeanSoftLuminosity'] = get_lum(RS_fact *np.mean(prof[selection] , axis=0), np.mean(r500c[selection]))*(u.erg/u.s)
        t_out[M_str + 'MeanHaloM500c'] = np.mean(MD40['HALO_M500c'][selection])*(u.Msun)
        t_out[M_str + 'StdHaloM500c'] = np.std(MD40['HALO_M500c'][selection])*(u.Msun)
        #print('='*100)
        #print(M0, nl(selection), 10**np.mean(MD40['CLUSTER_LX_soft_RF'][selection]), t_out[M_str + 'LUM_Eq11_2dProj_mean'][0], t_out[M_str + 'LUM_Eq11_2dProj_mean'][0]/10**np.mean(MD40['CLUSTER_LX_soft_RF'][selection]))
        #print('='*100)

    t_out.write(os.path.join(fig_dir, 'SimulatedProfiles_'+z_str+'_QN.fits'), overwrite = True)


