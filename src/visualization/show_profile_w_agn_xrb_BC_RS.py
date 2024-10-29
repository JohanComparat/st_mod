"""
Stacks data around positions and redshifts of GAMA galaxies


python cube_plot_profile_and_spectrum.py  CEN_bins_-14.0_ssfr_-11.0_10.0_mass_11.0 M1
python cube_plot_profile_and_spectrum.py  CEN_bins_-11.0_ssfr_-8.0_10.4_mass_11.0 M1


"""
sys.path.append( os.path.join(os.environ['GIT_STMOD'], 'src') )
from io import *

Z_SEL = sys.argv[1]

prof_dir = os.path.join( os.environ['GIT_STACK'], 'figures', 'Ti20_SDSS_stacked_galaxy_profile', 'profiles_decomp_AXH_BCRS' )
os.system( 'mkdir -p ' + prof_dir )

#M_val = 14.0
LX_outs = []
ALLT = []
t_out = Table()
t_out['AGN_frac_QU'] = n.array([ 0.2, 0.2, 0.2, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09])
for f_AGN_QU, M_val in zip(t_out['AGN_frac_QU'],n.arange(11.0, 15.5, 0.5)): #= 11.0
        #arr = [M_val, M_val+0.5]
        print('='*100)
        print('='*100)
        print(M_val)
        p2_prof = glob.glob(os.path.join(mergedCube_dir,'CEN-ANY-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0]
        #s_cat = glob.glob(os.path.join(simulated_directory,'CEN-ANY-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'-SIM.fits' ))[0]
        CEN_ANY = get_profile_SDSS_Ti21(p2_prof)#, s_cat)

        p2_prof = glob.glob(os.path.join(mergedCube_dir,'CEN-RS-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0]
        #s_cat = glob.glob(os.path.join(simulated_directory,'CEN-RS-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'-SIM.fits' ))[0]
        CEN_RS = get_profile_SDSS_Ti21(p2_prof)#, s_cat)

        p2_prof = glob.glob(os.path.join(mergedCube_dir,'CEN-BC-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'_STACKEDprofiles.fits' ))[0]
        #s_cat = glob.glob(os.path.join(simulated_directory,'CEN-BC-0.01_'+Z_SEL+'_*-'+str(n.round(M_val, 1))+'_Mhalo_'+str(n.round(M_val+0.5, 1))+'-SIM.fits' ))[0]
        CEN_BC = get_profile_SDSS_Ti21(p2_prof)#, s_cat)

        p2_fig = os.path.join(prof_dir, 'CEN_BC_RS_AGN_'+Z_SEL+'_Mhalo_'+str(M_val)+'.png')
        p.figure(33, (5, 4.5) )
        # ALL
        str_NG = ', N='+str(CEN_ANY['N_gal'][0])
        str_MS = ', M$^*$='+str(n.round(CEN_ANY['MS_mean'][0],1))+'$\pm$'+str(n.round(CEN_ANY['MS_std'][0],1))
        N_all = CEN_ANY['N_gal'][0]
        p.axvline(CEN_ANY['R500c'][0], ls='--', lw=2, c='k')#, label=r'$\bar{R}_{500c}=$'+str(int(CEN_ANY['R500c'][0])) +'kpc')
        # ANY
        p.errorbar( CEN_ANY['dd_xb'], CEN_ANY['dd_profile_BGSUB'],
                yerr = [ CEN_ANY['dd_profile_BGSUB'] - CEN_ANY['dd_profile_lo_BGSUB'], CEN_ANY['dd_profile_up_BGSUB'] - CEN_ANY['dd_profile_BGSUB'] ],
                xerr = [ CEN_ANY['dd_xb'] - CEN_ANY['x_lo'], CEN_ANY['x_up'] - CEN_ANY['dd_xb']],
                lw=2, ls='', color=cs ['CEN_ANY'], label='CEN ANY')# + str_NG + str_MS )
        #
        # BC
        #
        p.errorbar( CEN_BC['dd_xb'], CEN_BC['dd_profile_BGSUB'],
                yerr = [ CEN_BC['dd_profile_BGSUB'] - CEN_BC['dd_profile_lo_BGSUB'], CEN_BC['dd_profile_up_BGSUB'] - CEN_BC['dd_profile_BGSUB'] ],
                xerr = [ CEN_BC['dd_xb'] - CEN_BC['x_lo'], CEN_BC['x_up'] - CEN_BC['dd_xb']],
                lw=1.5, ls='', color=cs ['CEN_BC'], label='CEN BC')# + str_NG + str_MS )
        # PSF
        p.step(CEN_BC['dd_xb'], CEN_BC['ps_profile_normed_BGSUB'], lw=0.8 , color='darkblue', where='mid')#, label='PSF profile' )

        Max_values = CEN_BC['dd_profile_up_BGSUB']
        Min_values = CEN_BC['dd_profile_lo_BGSUB']
        Mean_value = CEN_BC['dd_profile_BGSUB']
        Max_values[Max_values<=1]=1
        Min_values[Min_values<=1]=1
        Mean_value[Mean_value<=1]=1
        area_shell = ( CEN_BC['x_up']**2 -  CEN_BC['x_lo']**2 ) * n.pi
        LX_mean    = interp1d( CEN_BC['x_up'], n.cumsum(Mean_value*area_shell) ) ( CEN_BC['R500c'][0] )
        LX_mean_up = interp1d( CEN_BC['x_up'], n.cumsum(Max_values*area_shell) ) ( CEN_BC['R500c'][0] )
        LX_mean_lo = interp1d( CEN_BC['x_up'], n.cumsum(Min_values*area_shell) ) ( CEN_BC['R500c'][0] )
        CEN_BC['total_LX_R500c']    = LX_mean
        CEN_BC['total_LX_R500c_up'] = LX_mean_up
        CEN_BC['total_LX_R500c_lo'] = LX_mean_lo

        #
        # RS
        #
        p.errorbar( CEN_RS['dd_xb'], CEN_RS['dd_profile_BGSUB'],
                yerr = [ CEN_RS['dd_profile_BGSUB'] - CEN_RS['dd_profile_lo_BGSUB'], CEN_RS['dd_profile_up_BGSUB'] - CEN_RS['dd_profile_BGSUB'] ],
                xerr = [ CEN_RS['dd_xb'] - CEN_RS['x_lo'], CEN_RS['x_up'] - CEN_RS['dd_xb']],
                lw=1.5, ls='', color=cs ['CEN_RS'], label='CEN RS')# + str_NG + str_MS )
        # PSF
        p.step(CEN_RS['dd_xb'], CEN_RS['ps_profile_normed_BGSUB'], lw=0.8 , color='darkred', where='mid')#, label='PSF profile' )

        Max_values = CEN_RS['dd_profile_up_BGSUB']
        Min_values = CEN_RS['dd_profile_lo_BGSUB']
        Mean_value = CEN_RS['dd_profile_BGSUB']
        Max_values[Max_values<=1]=1
        Min_values[Min_values<=1]=1
        Mean_value[Mean_value<=1]=1
        area_shell = ( CEN_RS['x_up']**2 -  CEN_RS['x_lo']**2 ) * n.pi
        LX_mean    = interp1d( CEN_RS['x_up'], n.cumsum(Mean_value*area_shell) ) ( CEN_RS['R500c'][0] )
        LX_mean_up = interp1d( CEN_RS['x_up'], n.cumsum(Max_values*area_shell) ) ( CEN_RS['R500c'][0] )
        LX_mean_lo = interp1d( CEN_RS['x_up'], n.cumsum(Min_values*area_shell) ) ( CEN_RS['R500c'][0] )
        CEN_RS['total_LX_R500c']    = LX_mean
        CEN_RS['total_LX_R500c_up'] = LX_mean_up
        CEN_RS['total_LX_R500c_lo'] = LX_mean_lo

        # hot gas BC
        Max_values = CEN_BC['dd_profile_up_BGSUB'] - ( CEN_ANY['AGN_profile_normed_BGSUB']*0.75 + CEN_ANY['XRB_profile_normed_BGSUB']*0.75 ) * 4 / 5.
        Min_values = CEN_BC['dd_profile_lo_BGSUB'] - ( CEN_ANY['AGN_profile_normed_BGSUB']*1.25 + CEN_ANY['XRB_profile_normed_BGSUB']*1.25 ) * 4 / 5.
        Mean_value = CEN_BC['dd_profile_BGSUB']    - ( CEN_ANY['AGN_profile_normed_BGSUB']      + CEN_ANY['XRB_profile_normed_BGSUB']      ) * 4 / 5.
        Max_values[Max_values<=1]=1
        Min_values[Min_values<=1]=1
        Mean_value[Mean_value<=1]=1
        #area_shell = ( CEN_ANY['x_up']**2 -  CEN_ANY['x_lo']**2 ) * n.pi
        p.errorbar( CEN_ANY['dd_xb'], Mean_value,
                yerr = [ Mean_value - Min_values, Max_values - Mean_value ],
                xerr = [ CEN_ANY['dd_xb'] - CEN_ANY['x_lo'], CEN_ANY['x_up'] - CEN_ANY['dd_xb']],
                lw=1., ls='', color=cs ['CEN_BC_hotGAS'], label='Hot gas BC')#, L$_X$='+str(n.round(log10_LX_mean,2))+r'$\pm$' +str(n.round(log10_DLX,2)) )

        # hot gas RS
        Max_values = CEN_RS['dd_profile_up_BGSUB'] - ( CEN_ANY['AGN_profile_normed_BGSUB']*0.75 + CEN_ANY['XRB_profile_normed_BGSUB']*0.75 ) * 1 / 5.
        Min_values = CEN_RS['dd_profile_lo_BGSUB'] - ( CEN_ANY['AGN_profile_normed_BGSUB']*1.25 + CEN_ANY['XRB_profile_normed_BGSUB']*1.25 ) * 1 / 5.
        Mean_value = CEN_RS['dd_profile_BGSUB']    - ( CEN_ANY['AGN_profile_normed_BGSUB']      + CEN_ANY['XRB_profile_normed_BGSUB']      ) * 1 / 5.
        Max_values[Max_values<=1]=1
        Min_values[Min_values<=1]=1
        Mean_value[Mean_value<=1]=1
        #area_shell = ( CEN_ANY['x_up']**2 -  CEN_ANY['x_lo']**2 ) * n.pi
        p.errorbar( CEN_ANY['dd_xb'], Mean_value,
                yerr = [ Mean_value - Min_values, Max_values - Mean_value ],
                xerr = [ CEN_ANY['dd_xb'] - CEN_ANY['x_lo'], CEN_ANY['x_up'] - CEN_ANY['dd_xb']],
                lw=1., ls='', color=cs ['CEN_RS_hotGAS'], label='Hot gas RS')#, L$_X$='+str(n.round(log10_LX_mean,2))+r'$\pm$' +str(n.round(log10_DLX,2)) )

        p.xlabel(r'$R_p$ [kpc]')
        p.ylabel(r'$S_X$ $[erg\; s^{-1}\; kpc^{-2}]$')
        p.legend(fontsize=12, loc=3, ncol=2)
        p.xlim((4, 2000))
        p.ylim((CEN_BC['bg'][0]/5000., n.max(CEN_BC['dd_profile_BGSUB'])*1.5))
        p.yscale('log')
        p.xscale('log')
        p.title(str(M_val)+'$<\log_{10}(M_{200m}/M_\odot)<$'+str(M_val+0.5))
        p.tight_layout()
        p.savefig( p2_fig )
        p.clf()
        print(p2_fig, 'written')
        #cumulative LX integral interpolation
        ALLT.append(CEN_RS)
        ALLT.append(CEN_BC)

MERGE = Table(vstack((ALLT)))
MERGE.write(os.path.join(prof_dir, 'FULL_SUMMARY_Mhalobin_centrals_BC_RS_'+Z_SEL+'.fits'), overwrite = True)

##LX_outs = n.transpose(LX_outs)
##t = Table()
##t['M0'] = LX_outs[0]
##t['M1'] = LX_outs[1]
##t['Ng'] = LX_outs[2]
##t['R500c'] = LX_outs[3]
##t['R200b'] = LX_outs[4]
##t['Rvir']  = LX_outs[5]
##t['AGN_LX'] = LX_outs[6]
##t['AGN_frac'] = LX_outs[7]
##t['XRB_LX'] = LX_outs[8]
##t['hotgas_1stBin_LX_mean'] = LX_outs[9]
##t['hotgas_1stBin_LX_mean_up'] = LX_outs[10]
##t['hotgas_1stBin_LX_mean_lo'] = LX_outs[11]
##t['hotgas_LX_mean'] = LX_outs[12]
##t['hotgas_LX_mean_up'] = LX_outs[13]
##t['hotgas_LX_mean_lo'] = LX_outs[14]
##t['M500c_mean'] = LX_outs[15]
##t['M500c_std'] = LX_outs[16]
##t['M500c_Q05'] = LX_outs[17]
##t['M500c_Q95'] = LX_outs[18]
##t['total_LX_mean'] = LX_outs[19]
##t['total_LX_mean_up'] = LX_outs[20]
##t['total_LX_mean_lo'] = LX_outs[21]
##t['z_mean'] = LX_outs[22]
##t['z_std'] = LX_outs[23]
##t.write(os.path.join(prof_dir, 'LX_summary_Mhalobin_centrals_'+Z_SEL+'.fits'), overwrite = True)

##t = Table()
##t['M0'] = LX_outs[0]
##t['M1'] = LX_outs[1]
##t['Ng'] = LX_outs[2].astype('int')
##t['R500c'] = n.round(LX_outs[3],1)
##t['AGN_LX'] = n.round(LX_outs[6],2)
##t['AGN_frac'] = n.round(LX_outs[7],1)
##t['XRB_LX'] = n.round(LX_outs[8],2)
###t['hotgas_1stBin_LX_mean'] = n.round(LX_outs[9]/1e39,2)
###t['hotgas_1stBin_LX_mean_up'] = n.round(LX_outs[10]/1e39,2)
###t['hotgas_1stBin_LX_mean_lo'] = n.round(LX_outs[11]/1e39,2)
##t['hotgas_LX_mean'] = n.round(LX_outs[12]/1e39,2)
##t['hotgas_LX_mean_up'] = n.round(LX_outs[13]/1e39,2)
##t['hotgas_LX_mean_lo'] = n.round(LX_outs[14]/1e39,2)
##t.write(os.path.join(prof_dir, 'LX_summary_Mhalobin_centrals_'+Z_SEL+'.latex'), format='latex', overwrite = True)



