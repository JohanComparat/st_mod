"""
Creates figures with the decomposition of the stacked signal into centrals and satellites
"""
sys.path.append( os.path.join(os.environ['GIT_STMOD'], 'src') )
from io.Profile import *

Z_SEL = 'Z3P'

prof_dir = os.path.join( os.environ['GIT_STMOD'], 'figures', 'Ti20_SDSS_stacked_galaxy_profile', 'profiles_decomp_Mh' )
os.system( 'mkdir -p ' + prof_dir )

fit_dir = os.path.join( os.environ['GIT_STMOD_DATA'], 'data', 'results', 'Ti20_SDSS_stacked_galaxy_profile', 'profiles_decomp_Mh' )
os.system( 'mkdir -p ' + fit_dir )
p2_file_fit_out = os.path.join(fit_dir, 'MhaloSel_ALL_SATfrac_params.fits')

prefix = 'Ti20_SDSS_kdgroups'
root_directory = os.environ['DATA_S4'] # /ptmp/joco/mpecl/comparat/data_s4/
mergedCube_dir = os.path.join(root_directory, 'mergedCubes', prefix)
mergedCube_500E1000_dir  = os.path.join(root_directory, 'mergedCubes_500_1000', prefix)
mergedCube_1000E1200_dir = os.path.join(root_directory, 'mergedCubes_1000_2000', prefix)
galaxy_directory = os.path.join(root_directory, 'galaxy_catalogues', prefix )
simulated_directory = os.path.join(root_directory, 'simulated_catalogues', prefix )



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
        ALL_ANY = get_profile (glob.glob(os.path.join(mergedCube_dir,'ALL-ANY-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )
        #
        CEN_ANY = get_profile (glob.glob(os.path.join(mergedCube_dir,'CEN-ANY-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )
        CEN_RS = get_profile (glob.glob(os.path.join(mergedCube_dir,'CEN-RS-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )
        CEN_BC = get_profile (glob.glob(os.path.join(mergedCube_dir,'CEN-BC-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )
        #
        SAT_ANY = get_profile (glob.glob(os.path.join(mergedCube_dir,'SAT-ANY-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )
        SAT_RS = get_profile (glob.glob(os.path.join(mergedCube_dir,'SAT-RS-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )
        SAT_BC = get_profile (glob.glob(os.path.join(mergedCube_dir,'SAT-BC-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )
        #

        # files for energy dependence

        Elo_ALL_ANY = get_profile (glob.glob(os.path.join(mergedCube_500E1000_dir,'ALL-ANY-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )
        Elo_CEN_ANY = get_profile (glob.glob(os.path.join(mergedCube_500E1000_dir,'CEN-ANY-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )
        Elo_CEN_RS  = get_profile (glob.glob(os.path.join(mergedCube_500E1000_dir,'CEN-RS-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )
        Elo_CEN_BC  = get_profile (glob.glob(os.path.join(mergedCube_500E1000_dir,'CEN-BC-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )
        Elo_SAT_ANY = get_profile (glob.glob(os.path.join(mergedCube_500E1000_dir,'SAT-ANY-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )
        Elo_SAT_RS  = get_profile (glob.glob(os.path.join(mergedCube_500E1000_dir,'SAT-RS-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )
        Elo_SAT_BC  = get_profile (glob.glob(os.path.join(mergedCube_500E1000_dir,'SAT-BC-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )


        Eup_ALL_ANY = get_profile (glob.glob(os.path.join(mergedCube_1000E1200_dir,'ALL-ANY-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )
        Eup_CEN_ANY = get_profile (glob.glob(os.path.join(mergedCube_1000E1200_dir,'CEN-ANY-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )
        Eup_CEN_RS = get_profile (glob.glob(os.path.join(mergedCube_1000E1200_dir,'CEN-RS-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )
        Eup_CEN_BC = get_profile (glob.glob(os.path.join(mergedCube_1000E1200_dir,'CEN-BC-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )
        Eup_SAT_ANY = get_profile (glob.glob(os.path.join(mergedCube_1000E1200_dir,'SAT-ANY-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )
        Eup_SAT_RS = get_profile (glob.glob(os.path.join(mergedCube_1000E1200_dir,'SAT-RS-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )
        Eup_SAT_BC = get_profile (glob.glob(os.path.join(mergedCube_1000E1200_dir,'SAT-BC-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0] )

        #
        # Measurement
        #
        p2_fig = os.path.join(prof_dir, 'MhaloSel_ALL_CEN_SAT_RS_BC_'+Z_SEL+'_SBprofile_'+str(M_val)+'.png')
        p.figure(2, (5.,4.5) )
        # ALL
        p.step(n.ones_like(ALL_ANY['dd_xb'])-1000, n.ones_like(ALL_ANY['ps_profile_normed'])-1000., lw=0.8 , color='k', where='mid', label='PSF profile' )
        p.axvline(CEN_ANY['R500c'][0], ls='--', lw=2, c='k', label='$R_{500c}=$'+str(n.round(CEN_ANY['R500c'][0],1))+'kpc')

        str_NG = ', N='+str(ALL_ANY['N_gal'][0])
        str_MS = ', M$^*$='+str(n.round(ALL_ANY['GAL_stellarMass_mean'][0],1))+'$\pm$'+str(n.round(ALL_ANY['GAL_stellarMass_std'][0],1))
        N_all = ALL_ANY['N_gal'][0]
        p.errorbar( ALL_ANY['dd_xb'], ALL_ANY['dd_profile'],
                yerr = ALL_ANY['dd_profile_err'],
                xerr = [ ALL_ANY['dd_xb'] - ALL_ANY['x_lo'], ALL_ANY['x_up'] - ALL_ANY['dd_xb']],
                lw=2, ls='', color=cs['ALL_ANY'], label='ALL')# + str_NG + str_MS)
        p.step(ALL_ANY['dd_xb'], ALL_ANY['ps_profile_normed'], lw=0.8 , color=cs['ALL_ANY'], where='mid' )
        ## CEN
        #str_NG = ', N='+str(CEN_ANY['N_gal'][0])
        #str_MS = ', M$^*$='+str(n.round(CEN_ANY['GAL_stellarMass_mean'][0],1))+'$\pm$'+str(n.round(CEN_ANY['GAL_stellarMass_std'][0],1))
        #N_cen = CEN_ANY['N_gal'][0]
        #frac_signal = N_cen*1./N_all
        #p.errorbar( CEN_ANY['dd_xb'], CEN_ANY['dd_profile'],
                #yerr = CEN_ANY['dd_profile_err'],
                #xerr = [ CEN_ANY['dd_xb'] - CEN_ANY['x_lo'], CEN_ANY['x_up'] - CEN_ANY['dd_xb']],
                #lw=2, ls='', color=cs['CEN_ANY'], label='Cen')# + str_NG + str_MS)
        #p.step(CEN_ANY['dd_xb'], CEN_ANY['ps_profile_normed'], lw=0.8 , color=cs['CEN_ANY'], where='mid' )
        # CEN RS
        str_NG = ', N='+str(CEN_RS['N_gal'][0])
        str_MS = ', M$^*$='+str(n.round(CEN_RS['GAL_stellarMass_mean'][0],1))+'$\pm$'+str(n.round(CEN_RS['GAL_stellarMass_std'][0],1))
        N_cen = CEN_RS['N_gal'][0]
        frac_signal = N_cen*1./N_all
        p.errorbar( CEN_RS['dd_xb'], CEN_RS['dd_profile'],
                yerr = CEN_RS['dd_profile_err'],
                xerr = [ CEN_RS['dd_xb'] - CEN_RS['x_lo'], CEN_RS['x_up'] - CEN_RS['dd_xb']],
                lw=2, ls='', color=cs['CEN_RS'], label='Cen RS')# + str_NG + str_MS)
        p.step(CEN_RS['dd_xb'], CEN_RS['ps_profile_normed'], lw=0.8 , color=cs['CEN_RS'], where='mid' )
        # CEN BC
        str_NG = ', N='+str(CEN_BC['N_gal'][0])
        str_MS = ', M$^*$='+str(n.round(CEN_BC['GAL_stellarMass_mean'][0],1))+'$\pm$'+str(n.round(CEN_BC['GAL_stellarMass_std'][0],1))
        N_cen = CEN_BC['N_gal'][0]
        frac_signal = N_cen*1./N_all
        p.errorbar( CEN_BC['dd_xb'], CEN_BC['dd_profile'],
                yerr = CEN_BC['dd_profile_err'],
                xerr = [ CEN_BC['dd_xb'] - CEN_BC['x_lo'], CEN_BC['x_up'] - CEN_BC['dd_xb']],
                lw=2, ls='', color=cs['CEN_BC'], label='Cen BC')# + str_NG + str_MS)
        p.step(CEN_BC['dd_xb'], CEN_BC['ps_profile_normed'], lw=0.8 , color=cs['CEN_BC'], where='mid' )

        # SAT
        #N_sat = SAT_ANY['N_gal'][0]
        #frac_signal = N_sat*1./N_all
        #str_NG = ', N='+str(SAT_ANY['N_gal'][0])
        #str_MS = ', M$^*$='+str(n.round(SAT_ANY['GAL_stellarMass_mean'][0],1))+'$\pm$'+str(n.round(SAT_ANY['GAL_stellarMass_std'][0],1))
        #p.errorbar( SAT_ANY['dd_xb'], SAT_ANY['dd_profile'],
                #yerr = SAT_ANY['dd_profile_err'],
                #xerr = [ SAT_ANY['dd_xb'] - SAT_ANY['x_lo'], SAT_ANY['x_up'] - SAT_ANY['dd_xb']],
                #lw=2, ls='', color=cs['SAT_ANY'], label='Sat')# + str_NG + str_MS)
        #p.step(SAT_ANY['dd_xb'], SAT_ANY['ps_profile_normed'], lw=0.8 , color=cs['SAT_ANY'], where='mid' )

        # SAT
        N_sat = SAT_RS['N_gal'][0]
        frac_signal = N_sat*1./N_all
        str_NG = ', N='+str(SAT_RS['N_gal'][0])
        str_MS = ', M$^*$='+str(n.round(SAT_RS['GAL_stellarMass_mean'][0],1))+'$\pm$'+str(n.round(SAT_RS['GAL_stellarMass_std'][0],1))
        p.errorbar( SAT_RS['dd_xb'], SAT_RS['dd_profile'],
                yerr = SAT_RS['dd_profile_err'],
                xerr = [ SAT_RS['dd_xb'] - SAT_RS['x_lo'], SAT_RS['x_up'] - SAT_RS['dd_xb']],
                lw=2, ls='', color=cs['SAT_RS'], label='Sat RS')# + str_NG + str_MS)
        p.step(SAT_RS['dd_xb'], SAT_RS['ps_profile_normed'], lw=0.8 , color=cs['SAT_RS'], where='mid' )

        # SAT
        N_sat = SAT_BC['N_gal'][0]
        frac_signal = N_sat*1./N_all
        str_NG = ', N='+str(SAT_BC['N_gal'][0])
        str_MS = ', M$^*$='+str(n.round(SAT_BC['GAL_stellarMass_mean'][0],1))+'$\pm$'+str(n.round(SAT_BC['GAL_stellarMass_std'][0],1))
        p.errorbar( SAT_BC['dd_xb'], SAT_BC['dd_profile'],
                yerr = SAT_BC['dd_profile_err'],
                xerr = [ SAT_BC['dd_xb'] - SAT_BC['x_lo'], SAT_BC['x_up'] - SAT_BC['dd_xb']],
                lw=2, ls='', color=cs['SAT_BC'], label='Sat BC')# + str_NG + str_MS)
        p.step(SAT_BC['dd_xb'], SAT_BC['ps_profile_normed'], lw=0.8 , color=cs['SAT_BC'], where='mid' )

        p.xlabel(r'$R_p$ [kpc]')
        p.ylabel(r'$S_X$ $[erg\; s^{-1}\; kpc^{-2}]$')
        p.legend(fontsize=12, loc=0, ncol=2)
        p.xlim((4, 2000))
        #p.ylim((4e36, 5e38))
        p.yscale('log')
        p.xscale('log')
        p.title(str(M_val)+'$<\log_{10}(M_{200m}/M_\odot)<$'+str(M_val+0.5), fontsize=12)
        p.tight_layout()
        p.savefig( p2_fig )
        p.clf()
        print(p2_fig, 'written')

        #
        # fraction
        #
        p2_fig = os.path.join(prof_dir, 'MhaloSel_ALL_CEN_SAT_RS_BC_'+Z_SEL+'_FRACprofile_'+str(M_val)+'.png')
        p.figure(2, (5.,4.5) )
        # ALL
        str_NG = ', N='+str(ALL_ANY['N_gal'][0])
        str_MS = ', M$^*$='+str(n.round(ALL_ANY['GAL_stellarMass_mean'][0],1))+'$\pm$'+str(n.round(ALL_ANY['GAL_stellarMass_std'][0],1))
        N_all = ALL_ANY['N_gal'][0]
        #p.errorbar( ALL_ANY['dd_xb'], ALL_ANY['dd_profile_BGSUB']/ALL_ANY['dd_profile_BGSUB'],
                ##yerr = (ALL_ANY['dd_profile_up_BGSUB']-ALL_ANY['dd_profile_lo_BGSUB'])*0.5/ALL_ANY['dd_profile'],
                #xerr = [ ALL_ANY['dd_xb'] - ALL_ANY['x_lo'], ALL_ANY['x_up'] - ALL_ANY['dd_xb']],
                #lw=2, ls='', color=cs['ALL_ANY'], label='ALL')# + str_NG + str_MS)

        # CEN
        str_NG = ', N='+str(CEN_RS['N_gal'][0])
        str_MS = ', M$^*$='+str(n.round(CEN_RS['GAL_stellarMass_mean'][0],1))+'$\pm$'+str(n.round(CEN_RS['GAL_stellarMass_std'][0],1))
        N_cen = CEN_RS['N_gal'][0]
        frac_signal = N_cen*1./N_all
        p.errorbar( CEN_RS['dd_xb'], CEN_RS['dd_profile_BGSUB']*frac_signal/ALL_ANY['dd_profile_BGSUB'],
                yerr = (CEN_RS['dd_profile_up_BGSUB']-CEN_RS['dd_profile_lo_BGSUB'])*0.5/CEN_RS['dd_profile_BGSUB'],
                xerr = [ CEN_RS['dd_xb'] - CEN_RS['x_lo'], CEN_RS['x_up'] - CEN_RS['dd_xb']],
                lw=2.5, ls='', color=cs['CEN_RS'], label='CEN RS', fmt='o',markersize=8,markerfacecolor='none')# + str_NG + str_MS)

        # CEN
        str_NG = ', N='+str(CEN_BC['N_gal'][0])
        str_MS = ', M$^*$='+str(n.round(CEN_BC['GAL_stellarMass_mean'][0],1))+'$\pm$'+str(n.round(CEN_BC['GAL_stellarMass_std'][0],1))
        N_cen = CEN_BC['N_gal'][0]
        frac_signal = N_cen*1./N_all
        p.errorbar( CEN_BC['dd_xb'], CEN_BC['dd_profile_BGSUB']*frac_signal/ALL_ANY['dd_profile_BGSUB'],
                yerr = (CEN_BC['dd_profile_up_BGSUB']-CEN_BC['dd_profile_lo_BGSUB'])*0.5/CEN_BC['dd_profile_BGSUB'],
                xerr = [ CEN_BC['dd_xb'] - CEN_BC['x_lo'], CEN_BC['x_up'] - CEN_BC['dd_xb']],
                lw=2, ls='', color=cs['CEN_BC'], label='CEN BC', fmt='o',markersize=8,markerfacecolor='none')# + str_NG + str_MS)

        # SAT
        N_sat = SAT_RS['N_gal'][0]
        frac_signal = N_sat*1./N_all
        str_NG = ', N='+str(SAT_RS['N_gal'][0])
        str_MS = ', M$^*$='+str(n.round(SAT_RS['GAL_stellarMass_mean'][0],1))+'$\pm$'+str(n.round(SAT_RS['GAL_stellarMass_std'][0],1))
        p.errorbar( SAT_RS['dd_xb'], SAT_RS['dd_profile_BGSUB']*frac_signal/ALL_ANY['dd_profile_BGSUB'],
                yerr = (SAT_RS['dd_profile_up_BGSUB']-SAT_RS['dd_profile_lo_BGSUB'])*0.5/SAT_RS['dd_profile_BGSUB'],
                xerr = [ SAT_RS['dd_xb'] - SAT_RS['x_lo'], SAT_RS['x_up'] - SAT_RS['dd_xb']],
                lw=1.5, ls='', color=cs['SAT_RS'], label='SAT RS', fmt='s',markersize=8,markerfacecolor='none')# + str_NG + str_MS)

        # SAT
        N_sat = SAT_BC['N_gal'][0]
        frac_signal = N_sat*1./N_all
        str_NG = ', N='+str(SAT_BC['N_gal'][0])
        str_MS = ', M$^*$='+str(n.round(SAT_BC['GAL_stellarMass_mean'][0],1))+'$\pm$'+str(n.round(SAT_BC['GAL_stellarMass_std'][0],1))
        p.errorbar( SAT_BC['dd_xb'], SAT_BC['dd_profile_BGSUB']*frac_signal/ALL_ANY['dd_profile_BGSUB'],
                yerr = (SAT_BC['dd_profile_up_BGSUB']-SAT_BC['dd_profile_lo_BGSUB'])*0.5/SAT_BC['dd_profile_BGSUB'],
                xerr = [ SAT_BC['dd_xb'] - SAT_BC['x_lo'], SAT_BC['x_up'] - SAT_BC['dd_xb']],
                lw=1, ls='', color=cs['SAT_BC'], label='SAT BC', fmt='s',markersize=8,markerfacecolor='none')# + str_NG + str_MS)

        p.xlabel(r'$R_p$ [kpc]')
        p.ylabel(r'$S_X$ fraction')# $[erg\; s^{-1}\; kpc^{-2}]$')
        p.legend(fontsize=12, loc=3, ncol=2)
        p.xlim((4, 2000))
        p.ylim((-0.25, 1.01)) #4e36, 5e38))
        p.xscale('log')
        p.tight_layout()
        p.savefig( p2_fig )
        p.clf()
        print(p2_fig, 'written')


        #
        # fraction
        #
        p2_fig = os.path.join(prof_dir, 'MhaloSel_ALL_CEN_SAT_only_'+Z_SEL+'_FRACprofile_'+str(M_val)+'.png')
        p.figure(2, (5.,4.5) )
        # ALL
        str_NG = ', N='+str(ALL_ANY['N_gal'][0])
        str_MS = ', M$^*$='+str(n.round(ALL_ANY['GAL_stellarMass_mean'][0],1))+'$\pm$'+str(n.round(ALL_ANY['GAL_stellarMass_std'][0],1))
        N_all = ALL_ANY['N_gal'][0]
        p.errorbar( ALL_ANY['dd_xb'], ALL_ANY['dd_profile_BGSUB']/ALL_ANY['dd_profile_BGSUB'],
                #yerr = ALL_ANY['dd_profile_BGSUB_err']/ALL_ANY['dd_profile_BGSUB'],
                xerr = [ ALL_ANY['dd_xb'] - ALL_ANY['x_lo'], ALL_ANY['x_up'] - ALL_ANY['dd_xb']],
                lw=2, ls='', color=cs['ALL_ANY'])#, label='ALL')# + str_NG + str_MS)
        # CEN
        str_NG = ', N='+str(CEN_ANY['N_gal'][0])
        str_MS = ', M$^*$='+str(n.round(CEN_ANY['GAL_stellarMass_mean'][0],1))+'$\pm$'+str(n.round(CEN_ANY['GAL_stellarMass_std'][0],1))
        N_cen = CEN_ANY['N_gal'][0]
        frac_signal = N_cen*1./N_all
        p.errorbar( CEN_ANY['dd_xb'], CEN_ANY['dd_profile_BGSUB']*frac_signal/ALL_ANY['dd_profile_BGSUB'],
                yerr = (CEN_ANY['dd_profile_up_BGSUB']-CEN_ANY['dd_profile_lo_BGSUB'])*0.5/CEN_ANY['dd_profile_BGSUB'],
                xerr = [ CEN_ANY['dd_xb'] - CEN_ANY['x_lo'], CEN_ANY['x_up'] - CEN_ANY['dd_xb']],
                lw=2, ls='', color=cs['CEN_ANY'])#, label='Cen')# + str_NG + str_MS)
        inner_val_cen = CEN_ANY['dd_profile_BGSUB'][0]*frac_signal/CEN_ANY['dd_profile_BGSUB'][0]
        out_val_cen = n.mean(CEN_ANY['dd_profile_BGSUB'][-4:]*frac_signal/CEN_ANY['dd_profile_BGSUB'][-4:])#*frac_signal

        # SAT
        N_sat = SAT_ANY['N_gal'][0]
        frac_signal = N_sat*1./N_all
        str_NG = ', N='+str(SAT_ANY['N_gal'][0])
        str_MS = ', M$^*$='+str(n.round(SAT_ANY['GAL_stellarMass_mean'][0],1))+'$\pm$'+str(n.round(SAT_ANY['GAL_stellarMass_std'][0],1))
        p.errorbar( SAT_ANY['dd_xb'], SAT_ANY['dd_profile_BGSUB']*frac_signal/ALL_ANY['dd_profile_BGSUB'],
                yerr = (SAT_ANY['dd_profile_up_BGSUB']-SAT_ANY['dd_profile_lo_BGSUB'])*0.5/SAT_ANY['dd_profile_BGSUB'],
                xerr = [ SAT_ANY['dd_xb'] - SAT_ANY['x_lo'], SAT_ANY['x_up'] - SAT_ANY['dd_xb']],
                lw=1.5, ls='', color=cs['SAT_ANY'])#, label='Sat')# + str_NG + str_MS)


        x_data = n.log10(SAT_ANY['dd_xb'])
        y_data = SAT_ANY['dd_profile_BGSUB']*frac_signal/ALL_ANY['dd_profile_BGSUB']
        y_data_err = (SAT_ANY['dd_profile_up_BGSUB']-SAT_ANY['dd_profile_lo_BGSUB'])*0.5/SAT_ANY['dd_profile_BGSUB']
        inner_val_sat = SAT_ANY['dd_profile_BGSUB'][0]*frac_signal/ALL_ANY['dd_profile_BGSUB'][0]
        min_val = n.min([1-inner_val_cen, inner_val_sat])
        out_val_sat = n.mean(SAT_ANY['dd_profile_BGSUB'][-4:]*frac_signal/ALL_ANY['dd_profile_BGSUB'][-4:])#*frac_signal
        max_val = n.max([1-out_val_cen, out_val_sat])
        if max_val>=1:
                max_val=1
        if max_val<0.5:
                max_val=0.5
        if min_val<=0:
                min_val=0
        log10r0 = n.log10( 0.35*CEN_ANY['R500c'][0])
        sigma = 0.3
        def fun2fit(log10radius, log10r0, sigma, min_val, max_val):
                return f_sat_mod(log10radius, log10r0, sigma, min_val, max_val)
        #out_p, out_cv = curve_fit(fun2fit, x_data, y_data, sigma= y_data_err, p0=(log10r0, sigma) )
        out_p, out_cv = curve_fit(fun2fit, x_data, y_data, sigma= y_data_err, p0=(log10r0, sigma, min_val, max_val), maxfev = 1000000 )
        label_eq = r'$\left(\frac{1}{2} + \frac{1}{2} erf\left(\frac{\log_{10}(R_p)-\log_{10}(R_0)}{\sigma}\right)\right)\times y_\Delta + y_0$'
        label_fit = r'$\log_{10}(R_0)=$'+str(n.round(out_p[0],2))+r', $\sigma=$'+str(n.round(out_p[1],2))
        label_fit2 = r'$y_\Delta=$'+str(n.round(out_p[3] - out_p[2],2))+r', $y_0=$'+str(n.round(out_p[2],2))
        def best_fit(log10radius):
                return f_sat_mod(log10radius, out_p[0], out_p[1], out_p[2], out_p[3])
        log10radius = n.arange(0.7, 3.4, 0.01)
        p.plot(10**log10radius, best_fit(log10radius), 'b--', lw=0.5, label=label_eq+'\n'+label_fit+'\n'+label_fit2)
        t_out['log10r0']     [jjj]= out_p[0]
        t_out['sigma']       [jjj]= out_p[1]
        t_out['min_val']     [jjj]= out_p[2]
        t_out['max_val']     [jjj]= out_p[3]
        t_out['err_log10r0'] [jjj]= out_cv[0][0]**0.5
        t_out['err_sigma']   [jjj]= out_cv[1][1]**0.5
        t_out['err_min_val'] [jjj]= out_cv[2][2]**0.5
        t_out['err_max_val'] [jjj]= out_cv[3][3]**0.5

        p.xlabel(r'$R_p$ [kpc]')
        #p.xlabel(r'$R_p/(R_{500c}=$'+str(n.round(CEN_ANY['R500c'][0],1))+'kpc)')
        p.ylabel(r'$S_X$ fraction')# $[erg\; s^{-1}\; kpc^{-2}]$')
        p.legend(fontsize=12, loc=4, ncol=1)
        p.xlim((4, 2000))
        p.ylim((-0.25, 1.01)) #4e36, 5e38))
        #p.yscale('log')
        p.xscale('log')
        #p.title(str(M_val)+'$<\log_{10}(M_{200m}/M_\odot)<$'+str(M_val+0.5)+', R$_{vir}=$'+str(n.round(CEN_ANY['R500c'][0],1))+'kpc')
        p.tight_layout()
        p.savefig( p2_fig )
        p.clf()
        print(p2_fig, 'written')


        #
        # Measurement BG subtracted, decomposition
        #
        p2_fig = os.path.join(prof_dir, 'MhaloSel_ALL_CEN_SAT_RS_BC_'+Z_SEL+'_BGsubSBprofile_'+str(M_val)+'.png')
        p.figure(2, (5.,4.5) )
        # ALL
        #p.step(n.ones_like(ALL_ANY['dd_xb'])-1000, n.ones_like(ALL_ANY['ps_profile_normed'])-1000., lw=0.8 , color='k', where='mid', label='PSF profile' )
        p.axhline(ALL_ANY['bg'][0]/300., ls='--', lw=2, c='grey', label='background/300')
        p.axvline(CEN_ANY['R500c'][0], ls='--', lw=2, c='k', label='$R_{500c}=$'+str(n.round(CEN_ANY['R500c'][0],1))+'kpc')
        str_NG = ', N='+str(ALL_ANY['N_gal'][0])
        str_MS = ', M$^*$='+str(n.round(ALL_ANY['GAL_stellarMass_mean'][0],1))+'$\pm$'+str(n.round(ALL_ANY['GAL_stellarMass_std'][0],1))
        N_all = ALL_ANY['N_gal'][0]
        p.errorbar( ALL_ANY['dd_xb'], ALL_ANY['dd_profile_BGSUB'],
                yerr = [ ALL_ANY['dd_profile_BGSUB'] - ALL_ANY['dd_profile_lo_BGSUB'], ALL_ANY['dd_profile_up_BGSUB'] - ALL_ANY['dd_profile_BGSUB'] ],
                xerr = [ ALL_ANY['dd_xb'] - ALL_ANY['x_lo'], ALL_ANY['x_up'] - ALL_ANY['dd_xb']],
                lw=2, ls='', color=cs['ALL_ANY'], label='ALL')# + str_NG + str_MS)
        #p.step(ALL_ANY['dd_xb'], ALL_ANY['ps_profile_normed_BGSUB'], lw=0.8 , color=cs['ALL_ANY'], where='mid', label='PSF profile' )

        # CEN RS
        str_NG = ', N='+str(CEN_RS['N_gal'][0])
        str_MS = ', M$^*$='+str(n.round(CEN_RS['GAL_stellarMass_mean'][0],1))+'$\pm$'+str(n.round(CEN_RS['GAL_stellarMass_std'][0],1))
        N_cen = CEN_RS['N_gal'][0]
        frac_signal = N_cen*1./N_all
        p.errorbar( CEN_RS['dd_xb'], CEN_RS['dd_profile_BGSUB']*frac_signal,
                yerr = [ CEN_RS['dd_profile_BGSUB']*frac_signal - CEN_RS['dd_profile_lo_BGSUB']*frac_signal, CEN_RS['dd_profile_up_BGSUB']*frac_signal - CEN_RS['dd_profile_BGSUB']*frac_signal ],
                xerr = [ CEN_RS['dd_xb'] - CEN_RS['x_lo'], CEN_RS['x_up'] - CEN_RS['dd_xb']],
                lw=2, ls='', color=cs['CEN_RS'], label='CEN RS')# + str_NG + str_MS)
        #p.step(CEN_RS['dd_xb'], CEN_RS['ps_profile_normed_BGSUB'], lw=0.8 , color=cs['CEN_RS'], where='mid', label='PSF profile' )
        # CEN BC
        str_NG = ', N='+str(CEN_BC['N_gal'][0])
        str_MS = ', M$^*$='+str(n.round(CEN_BC['GAL_stellarMass_mean'][0],1))+'$\pm$'+str(n.round(CEN_BC['GAL_stellarMass_std'][0],1))
        N_cen = CEN_BC['N_gal'][0]
        frac_signal = N_cen*1./N_all
        p.errorbar( CEN_BC['dd_xb'], CEN_BC['dd_profile_BGSUB']*frac_signal,
                yerr = [ CEN_BC['dd_profile_BGSUB']*frac_signal - CEN_BC['dd_profile_lo_BGSUB']*frac_signal, CEN_BC['dd_profile_up_BGSUB']*frac_signal - CEN_BC['dd_profile_BGSUB']*frac_signal ],
                xerr = [ CEN_BC['dd_xb'] - CEN_BC['x_lo'], CEN_BC['x_up'] - CEN_BC['dd_xb']],
                lw=2, ls='', color=cs['CEN_BC'], label='CEN BC')# + str_NG + str_MS)
        #p.step(CEN_BC['dd_xb'], CEN_BC['ps_profile_normed_BGSUB'], lw=0.8 , color=cs['CEN_BC'], where='mid', label='PSF profile' )

        # SAT
        N_sat = SAT_RS['N_gal'][0]
        frac_signal = N_sat*1./N_all
        str_NG = ', N='+str(SAT_RS['N_gal'][0])
        str_MS = ', M$^*$='+str(n.round(SAT_RS['GAL_stellarMass_mean'][0],1))+'$\pm$'+str(n.round(SAT_RS['GAL_stellarMass_std'][0],1))
        p.errorbar( SAT_RS['dd_xb'], SAT_RS['dd_profile_BGSUB']*frac_signal,
                yerr = [ SAT_RS['dd_profile_BGSUB']*frac_signal - SAT_RS['dd_profile_lo_BGSUB']*frac_signal, SAT_RS['dd_profile_up_BGSUB']*frac_signal - SAT_RS['dd_profile_BGSUB']*frac_signal ],
                xerr = [ SAT_RS['dd_xb'] - SAT_RS['x_lo'], SAT_RS['x_up'] - SAT_RS['dd_xb']],
                lw=2, ls='', color=cs['SAT_RS'], label='SAT RS')# + str_NG + str_MS)
        #p.step(SAT_RS['dd_xb'], SAT_RS['ps_profile_normed_BGSUB'], lw=0.8 , color=cs['SAT_RS'], where='mid', label='PSF profile' )

        # SAT
        N_sat = SAT_BC['N_gal'][0]
        frac_signal = N_sat*1./N_all
        str_NG = ', N='+str(SAT_BC['N_gal'][0])
        str_MS = ', M$^*$='+str(n.round(SAT_BC['GAL_stellarMass_mean'][0],1))+'$\pm$'+str(n.round(SAT_BC['GAL_stellarMass_std'][0],1))
        p.errorbar( SAT_BC['dd_xb'], SAT_BC['dd_profile_BGSUB']*frac_signal,
                yerr = [ SAT_BC['dd_profile_BGSUB']*frac_signal - SAT_BC['dd_profile_lo_BGSUB']*frac_signal, SAT_BC['dd_profile_up_BGSUB']*frac_signal - SAT_BC['dd_profile_BGSUB']*frac_signal ],
                xerr = [ SAT_BC['dd_xb'] - SAT_BC['x_lo'], SAT_BC['x_up'] - SAT_BC['dd_xb']],
                lw=2, ls='', color=cs['SAT_BC'], label='SAT BC')# + str_NG + str_MS)
        #p.step(SAT_BC['dd_xb'], SAT_BC['ps_profile_normed_BGSUB'], lw=0.8 , color=cs['SAT_BC'], where='mid', label='PSF profile' )

        p.xlabel(r'$R_p$ [kpc]')
        p.ylabel(r'$S_X$ $[erg\; s^{-1}\; kpc^{-2}]$')
        p.legend(fontsize=12, loc=3, ncol=2)
        p.xlim((4, 2000))
        p.ylim((ALL_ANY['bg'][0]/10000., n.max(ALL_ANY['dd_profile_BGSUB'])*1.5))
        p.yscale('log')
        p.xscale('log')
        p.title(str(M_val)+'$<\log_{10}(M_{200m}/M_\odot)<$'+str(M_val+0.5), fontsize=12)
        p.tight_layout()
        p.savefig( p2_fig )
        p.clf()
        print(p2_fig, 'written')


        #
        # Measurement BG subtracted, decopmosition
        #
        p2_fig = os.path.join(prof_dir, 'MhaloSel_ALL_CEN_SAT_only_'+Z_SEL+'_BGsubSBprofile_'+str(M_val)+'.png')
        p2_tex = os.path.join(prof_dir, 'MhaloSel_ALL_CEN_SAT_only_'+Z_SEL+'_BGsubSBprofile_'+str(M_val)+'.latex')
        p.figure(2, (5.,4.5) )
        # ALL
        p.step(n.ones_like(ALL_ANY['dd_xb'])-1000, n.ones_like(ALL_ANY['ps_profile_normed'])-1000., lw=0.8 , color='k', where='mid', label='PSF profile' )
        p.axhline(ALL_ANY['bg'][0]/300., ls='--', lw=2, c='grey', label='background/300')
        p.axvline(CEN_ANY['R500c'][0], ls='--', lw=2, c='k', label='$R_{500c}=$'+str(n.round(CEN_ANY['R500c'][0],1))+'kpc')
        str_NG = ', N='+str(ALL_ANY['N_gal'][0])
        str_MS = ', M$^*$='+str(n.round(ALL_ANY['GAL_stellarMass_mean'][0],1))+'$\pm$'+str(n.round(ALL_ANY['GAL_stellarMass_std'][0],1))
        N_all = ALL_ANY['N_gal'][0]
        p.errorbar( ALL_ANY['dd_xb'], ALL_ANY['dd_profile_BGSUB'],
                yerr = [ ALL_ANY['dd_profile_BGSUB'] - ALL_ANY['dd_profile_lo_BGSUB'], ALL_ANY['dd_profile_up_BGSUB'] - ALL_ANY['dd_profile_BGSUB'] ],
                xerr = [ ALL_ANY['dd_xb'] - ALL_ANY['x_lo'], ALL_ANY['x_up'] - ALL_ANY['dd_xb']],
                lw=2, ls='', color=cs['ALL_ANY'], label='ALL')# + str_NG + str_MS)
        #p.step(ALL_ANY['dd_xb'], ALL_ANY['ps_profile_normed_BGSUB'], lw=0.8 , color=cs['ALL_ANY'], where='mid', label='PSF profile' )
        # CEN
        str_NG = ', N='+str(CEN_ANY['N_gal'][0])
        str_MS = ', M$^*$='+str(n.round(CEN_ANY['GAL_stellarMass_mean'][0],1))+'$\pm$'+str(n.round(CEN_ANY['GAL_stellarMass_std'][0],1))
        N_cen = CEN_ANY['N_gal'][0]
        frac_signal = N_cen*1./N_all
        p.errorbar( CEN_ANY['dd_xb'], CEN_ANY['dd_profile_BGSUB'],
                yerr = [ CEN_ANY['dd_profile_BGSUB'] - CEN_ANY['dd_profile_lo_BGSUB'], CEN_ANY['dd_profile_up_BGSUB'] - CEN_ANY['dd_profile_BGSUB'] ],
                xerr = [ CEN_ANY['dd_xb'] - CEN_ANY['x_lo'], CEN_ANY['x_up'] - CEN_ANY['dd_xb']],
                lw=2, ls='', color=cs['CEN_ANY'], label='CEN ANY')# + str_NG + str_MS)
        p.step(CEN_ANY['dd_xb'], CEN_ANY['ps_profile_normed_BGSUB'], lw=0.8 , color=cs['CEN_ANY'], where='mid')#, label='PSF profile' )
        sss = (CEN_ANY['dd_xb']<=CEN_ANY['R500c'])&(CEN_ANY['dd_profile_BGSUB']>CEN_ANY['bg'][0]/300.)&(CEN_ANY['dd_xb']>10.)
        difference    = CEN_ANY['dd_profile_BGSUB'][sss]-CEN_ANY['ps_profile_normed_BGSUB'][sss]
        psf_err = CEN_ANY['ps_profile_normed_BGSUB'][sss]  *CEN_ANY['ps_profile_NN'][sss]**-0.5
        data_err = (CEN_ANY['dd_profile_up_BGSUB'][sss] - CEN_ANY['dd_profile_lo_BGSUB'][sss])/2.
        diff_err = ( psf_err**2 + data_err**2 ) **0.5
        significance = n.max( difference  / diff_err)
        print('='*100)
        print(M_val)
        print('='*100)
        print('significance max', significance, 'sigmas')
        print(n.transpose([n.round(CEN_ANY['x_lo'][sss],1), n.round(CEN_ANY['x_up'][sss],1), n.round(difference  / diff_err,2)]))
        print()
        print('='*100)
        print('='*100)

        t_CEN_prof = Table()
        t_CEN_prof['x_min']       = n.round(CEN_ANY['x_lo'][sss],1)
        t_CEN_prof['x_max']       = n.round(CEN_ANY['x_up'][sss],1)
        t_CEN_prof['S_x']        = n.round(CEN_ANY['dd_profile_BGSUB'][sss]/1e35                                                    ,2)
        t_CEN_prof['S_x_err']    = n.round(( CEN_ANY['dd_profile_up_BGSUB'][sss] - CEN_ANY['dd_profile_lo_BGSUB'][sss] ) /2. / 1e35 ,2)
        t_CEN_prof['S_psf']      = n.round(CEN_ANY['ps_profile_normed_BGSUB'][sss] / 1e35                                           ,2)
        t_CEN_prof['S_psf_err']  = n.round(CEN_ANY['ps_profile_normed_BGSUB'][sss] / 1e35  *CEN_ANY['ps_profile_NN'][sss]**-0.5     ,2)
        t_CEN_prof['sigdevpsf']  = n.round(difference  / diff_err,2)
        t_CEN_prof.write(p2_tex, format='latex', overwrite = True )
        print('direct sum ',n.sum(t_CEN_prof['sigdevpsf']))
        print('squared sum',n.sum(t_CEN_prof['sigdevpsf']**2)**0.5)
        print(p2_tex, 'written')

        # SAT
        N_sat = SAT_ANY['N_gal'][0]
        frac_signal = N_sat*1./N_all
        str_NG = ', N='+str(SAT_ANY['N_gal'][0])
        str_MS = ', M$^*$='+str(n.round(SAT_ANY['GAL_stellarMass_mean'][0],1))+'$\pm$'+str(n.round(SAT_ANY['GAL_stellarMass_std'][0],1))
        p.errorbar( SAT_ANY['dd_xb'], SAT_ANY['dd_profile_BGSUB'],
                yerr = [ SAT_ANY['dd_profile_BGSUB'] - SAT_ANY['dd_profile_lo_BGSUB'], SAT_ANY['dd_profile_up_BGSUB'] - SAT_ANY['dd_profile_BGSUB'] ],
                xerr = [ SAT_ANY['dd_xb'] - SAT_ANY['x_lo'], SAT_ANY['x_up'] - SAT_ANY['dd_xb']],
                lw=2, ls='', color=cs['SAT_ANY'], label='SAT ANY')# + str_NG + str_MS)
        p.step(SAT_ANY['dd_xb'], SAT_ANY['ps_profile_normed_BGSUB'], lw=0.8 , color=cs['SAT_ANY'], where='mid')#, label='PSF profile' )

        p.xlabel(r'$R_p$ [kpc]')
        p.ylabel(r'$S_X$ $[erg\; s^{-1}\; kpc^{-2}]$')
        p.legend(fontsize=12, loc=3, ncol=2)
        p.xlim((4, 2000))
        p.ylim((ALL_ANY['bg'][0]/5000., n.max(ALL_ANY['dd_profile_BGSUB'])*1.5))
        p.yscale('log')
        p.xscale('log')
        p.title(str(M_val)+'$<\log_{10}(M_{200m}/M_\odot)<$'+str(M_val+0.5), fontsize=12)
        p.tight_layout()
        p.savefig( p2_fig )
        p.clf()
        print(p2_fig, 'written')



        #
        # Measurement E DEP BG subtracted, decomposition
        #
        p2_fig = os.path.join(prof_dir, 'MhaloSel_ALL_CEN_Edep_'+Z_SEL+'_BGsubSBprofile_'+str(M_val)+'.png')
        p.figure(2, (5.,4.5) )

        # CEN 0.5<e<2 keV
        p.errorbar( CEN_ANY['dd_xb'], CEN_ANY['dd_profile_BGSUB'],
                yerr = [ CEN_ANY['dd_profile_BGSUB'] - CEN_ANY['dd_profile_lo_BGSUB'], CEN_ANY['dd_profile_up_BGSUB'] - CEN_ANY['dd_profile_BGSUB'] ],
                xerr = [ CEN_ANY['dd_xb'] - CEN_ANY['x_lo'], CEN_ANY['x_up'] - CEN_ANY['dd_xb']],
                lw=2, ls='', color=cs['CEN_ANY'], label='0.5<E<2')
        p.step(CEN_ANY['dd_xb'], CEN_ANY['ps_profile_normed_BGSUB'], lw=0.8 , color=cs['CEN_ANY'], where='mid')
        # CEN 0.5<e<1 keV
        p.errorbar( Elo_CEN_ANY['dd_xb'], Elo_CEN_ANY['dd_profile_BGSUB'],
                yerr = [ Elo_CEN_ANY['dd_profile_BGSUB'] - Elo_CEN_ANY['dd_profile_lo_BGSUB'], Elo_CEN_ANY['dd_profile_up_BGSUB'] - Elo_CEN_ANY['dd_profile_BGSUB'] ],
                xerr = [ Elo_CEN_ANY['dd_xb'] - Elo_CEN_ANY['x_lo'], Elo_CEN_ANY['x_up'] - Elo_CEN_ANY['dd_xb']],
                lw=2, ls='', color=cs['CEN_BC'], label='0.5<E<1')
        p.step(Elo_CEN_ANY['dd_xb'], Elo_CEN_ANY['ps_profile_normed_BGSUB'], lw=0.8 , color=cs['CEN_BC'], where='mid')
        # CEN 1<e<2 kev
        p.errorbar( Eup_CEN_ANY['dd_xb'], Eup_CEN_ANY['dd_profile_BGSUB'],
                yerr = [ Eup_CEN_ANY['dd_profile_BGSUB'] - Eup_CEN_ANY['dd_profile_lo_BGSUB'], Eup_CEN_ANY['dd_profile_up_BGSUB'] - Eup_CEN_ANY['dd_profile_BGSUB'] ],
                xerr = [ Eup_CEN_ANY['dd_xb'] - Eup_CEN_ANY['x_lo'], Eup_CEN_ANY['x_up'] - Eup_CEN_ANY['dd_xb']],
                lw=2, ls='', color=cs['CEN_RS'], label='1<E<2')
        p.step(Eup_CEN_ANY['dd_xb'], Eup_CEN_ANY['ps_profile_normed_BGSUB'], lw=0.8 , color=cs['CEN_RS'], where='mid')

        p.xlabel(r'$R_p$ [kpc]')
        p.ylabel(r'$S_X$ $[erg\; s^{-1}\; kpc^{-2}]$')
        p.legend(fontsize=12, loc=3, ncol=2, title='CEN')
        p.xlim((4, 2000))
        p.ylim((CEN_ANY['bg'][0]/10000., n.max(CEN_ANY['dd_profile_BGSUB'])*1.5))
        p.yscale('log')
        p.xscale('log')
        p.title(str(M_val)+'$<\log_{10}(M_{200m}/M_\odot)<$'+str(M_val+0.5), fontsize=12)
        p.tight_layout()
        p.savefig( p2_fig )
        p.clf()
        print(p2_fig, 'written')



t_out.write(p2_file_fit_out, overwrite = True)
