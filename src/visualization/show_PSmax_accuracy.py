sys.path.append( os.path.join(os.environ['GIT_STMOD'], 'src') )
from io import *

Z_SEL = sys.argv[1]

prof_dir = os.path.join( os.environ['GIT_STACK'], 'figures', 'Ti20_SDSS_stacked_galaxy_profile', 'PSF_max' )
os.system( 'mkdir -p ' + prof_dir )

fit_dir = os.path.join( os.environ['GIT_STACK'], 'data', 'Ti20_SDSS_stacked_galaxy_profile', 'PSF_max' )
os.system( 'mkdir -p ' + fit_dir )

M_bins = n.arange(11.0, 14.5, 0.5)
t_out = Table()
t_out['Mh0'] = M_bins
t_out['Mh1'] = M_bins+0.5
t_out['log10r0']     = 0.0
t_out['sigma']       = 0.0
t_out['min_val']     = 0.0
t_out['max_val']     = 0.0
t_out['err_log10r0'] = 0.0
t_out['err_sigma']   = 0.0
t_out['err_min_val'] = 0.0
t_out['err_max_val'] = 0.0
for jjj, M_val in enumerate(M_bins[2:]):
        #
        #ALL_ANY = get_profile (glob.glob(os.path.join(mergedCube_dir,'ALL-ANY-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )
        #
        CEN_ANY = get_profile (glob.glob(os.path.join(mergedCube_dir,'CEN-ANY-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )
        #CEN_RS = get_profile (glob.glob(os.path.join(mergedCube_dir,'CEN-RS-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )
        #CEN_BC = get_profile (glob.glob(os.path.join(mergedCube_dir,'CEN-BC-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )
        ##
        #SAT_ANY = get_profile (glob.glob(os.path.join(mergedCube_dir,'SAT-ANY-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )
        #SAT_RS = get_profile (glob.glob(os.path.join(mergedCube_dir,'SAT-RS-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )
        #SAT_BC = get_profile (glob.glob(os.path.join(mergedCube_dir,'SAT-BC-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )
        #
        # Sim profile comparat, eckert et al. 2022
        #
        #log10_m500c_min = n.log10(ALL_ANY['M500c_mean']-ALL_ANY['M500c_std'])[0]
        #log10_m500c_max = n.log10(ALL_ANY['M500c_mean']+ALL_ANY['M500c_std'])[0]
        #in_zbin = (allz_i <= ALL_ANY['z_max'][0]) & (allz_i >= ALL_ANY['z_min'][0]) & (n.log10(allm5_i) >= log10_m500c_min) & (n.log10(allm5_i) <= log10_m500c_max)

        # files for energy dependence

        #Elo_ALL_ANY = get_profile (glob.glob(os.path.join(mergedCube_500E1000_dir,'ALL-ANY-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )
        #Elo_CEN_ANY = get_profile (glob.glob(os.path.join(mergedCube_500E1000_dir,'CEN-ANY-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )
        #Elo_CEN_RS  = get_profile (glob.glob(os.path.join(mergedCube_500E1000_dir,'CEN-RS-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )
        #Elo_CEN_BC  = get_profile (glob.glob(os.path.join(mergedCube_500E1000_dir,'CEN-BC-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )
        #Elo_SAT_ANY = get_profile (glob.glob(os.path.join(mergedCube_500E1000_dir,'SAT-ANY-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )
        #Elo_SAT_RS  = get_profile (glob.glob(os.path.join(mergedCube_500E1000_dir,'SAT-RS-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )
        #Elo_SAT_BC  = get_profile (glob.glob(os.path.join(mergedCube_500E1000_dir,'SAT-BC-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )


        #Eup_ALL_ANY = get_profile (glob.glob(os.path.join(mergedCube_1000E1200_dir,'ALL-ANY-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )
        #Eup_CEN_ANY = get_profile (glob.glob(os.path.join(mergedCube_1000E1200_dir,'CEN-ANY-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )
        #Eup_CEN_RS = get_profile (glob.glob(os.path.join(mergedCube_1000E1200_dir,'CEN-RS-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )
        #Eup_CEN_BC = get_profile (glob.glob(os.path.join(mergedCube_1000E1200_dir,'CEN-BC-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )
        #Eup_SAT_ANY = get_profile (glob.glob(os.path.join(mergedCube_1000E1200_dir,'SAT-ANY-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )
        #Eup_SAT_RS = get_profile (glob.glob(os.path.join(mergedCube_1000E1200_dir,'SAT-RS-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )
        #Eup_SAT_BC = get_profile (glob.glob(os.path.join(mergedCube_1000E1200_dir,'SAT-BC-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )

        #
        # Measurement
        #
        bn = Z_SEL+'_Mh_'+str(M_val)
        p2_fig = os.path.join(prof_dir, 'PS_PSFmax-'+bn+'.png')
        p.figure(2, (5.,4.5) )
        # CEN
        p.step(n.ones_like(CEN_ANY['dd_xb'])-1000, n.ones_like(CEN_ANY['ps_profile_normed'])-1000., lw=0.8 , color='k', where='mid', label='PSF profile' )
        p.axvline(CEN_ANY['R500c'][0], ls='--', lw=2, c='k', label='$R_{500c}=$'+str(n.round(CEN_ANY['R500c'][0],1))+'kpc')

        str_NG = ', N='+str(CEN_ANY['N_gal'][0])
        p.errorbar( CEN_ANY['dd_xb'], CEN_ANY['dd_profile'],
                yerr = CEN_ANY['dd_profile_err'],
                xerr = [ CEN_ANY['dd_xb'] - CEN_ANY['x_lo'], CEN_ANY['x_up'] - CEN_ANY['dd_xb']],
                lw=2, ls='', color=cs['CEN_ANY'], label='CEN 0.5-2' + str_NG )
        p.step(CEN_ANY['dd_xb'], CEN_ANY['ps_profile_normed'], lw=0.8 , color=cs['CEN_ANY'], where='mid' )

        p.errorbar( CEN_ANY['dd_xb'], CEN_ANY['ps_profile'],
                yerr = CEN_ANY['ps_profile_err'],
                xerr = [ CEN_ANY['dd_xb'] - CEN_ANY['x_lo'], CEN_ANY['x_up'] - CEN_ANY['dd_xb']],
                lw=2, ls='', color=cs['PS_CEN_ANY'], label='PS for CEN' )

        p.axhline(CEN_ANY['bg'][0], ls='--', lw=2, c='grey', label='background')

        p.xlabel(r'$R_p$ [kpc]')
        p.ylabel(r'$S_X$ $[erg\; s^{-1}\; kpc^{-2}]$')
        p.legend(fontsize=10, loc=1, ncol=1)
        p.xlim((7, 3000))
        #p.ylim((4e36, 5e38))
        p.yscale('log')
        p.xscale('log')
        p.tight_layout()
        p.savefig( p2_fig )
        p.clf()
        print(p2_fig, 'written')

        #
        # Measurement BG subtracted
        #
        p2_fig = os.path.join(prof_dir, 'PS_BGsubSBprofile_'+bn+'.png')
        p.figure(2, (5.,4.5) )
        # ALL
        p.axhline(CEN_ANY['bg'][0]/300., ls='--', lw=2, c='grey', label='background/300')
        p.axvline(CEN_ANY['R500c'][0], ls='--', lw=2, c='k', label='$R_{500c}=$'+str(n.round(CEN_ANY['R500c'][0],1))+'kpc')
        N_all = CEN_ANY['N_gal'][0]
        p.errorbar( CEN_ANY['dd_xb'], CEN_ANY['dd_profile_BGSUB'],
                yerr = [ CEN_ANY['dd_profile_BGSUB'] - CEN_ANY['dd_profile_lo_BGSUB'], CEN_ANY['dd_profile_up_BGSUB'] - CEN_ANY['dd_profile_BGSUB'] ],
                xerr = [ CEN_ANY['dd_xb'] - CEN_ANY['x_lo'], CEN_ANY['x_up'] - CEN_ANY['dd_xb']],
                lw=2, ls='', color=cs['CEN_ANY'], label='CEN 0.5-2' + str_NG)
        p.step(CEN_ANY['dd_xb'], CEN_ANY['ps_profile_normed_BGSUB'], lw=0.8 , color=cs['CEN_ANY'], where='mid')

        p.errorbar( CEN_ANY['dd_xb'], CEN_ANY['ps_profile_BGSUB'],
                yerr = [ CEN_ANY['ps_profile_BGSUB'] - CEN_ANY['ps_profile_lo_BGSUB'], CEN_ANY['ps_profile_up_BGSUB'] - CEN_ANY['ps_profile_BGSUB'] ],
                xerr = [ CEN_ANY['dd_xb'] - CEN_ANY['x_lo'], CEN_ANY['x_up'] - CEN_ANY['dd_xb']],
                lw=2, ls='', color=cs['PS_CEN_ANY'], label='PS for CEN 0.5-2')
        s1 = (CEN_ANY['dd_xb']<CEN_ANY['R500c'][0])
        print(((CEN_ANY['ps_profile_up_BGSUB']-CEN_ANY['ps_profile_lo_BGSUB'])*0.5/CEN_ANY['ps_profile_BGSUB'])[s1])
        p.xlabel(r'$R_p$ [kpc]')
        p.ylabel(r'$S_X$ $[erg\; s^{-1}\; kpc^{-2}]$')
        p.legend(fontsize=10, loc=3, ncol=1)
        p.xlim((7, 3000))
        p.ylim((CEN_ANY['bg'][0]/10000., n.max(CEN_ANY['ps_profile_BGSUB'])*1.5))
        p.yscale('log')
        p.xscale('log')
        p.tight_layout()
        p.savefig( p2_fig )
        p.clf()
        print(p2_fig, 'written')
