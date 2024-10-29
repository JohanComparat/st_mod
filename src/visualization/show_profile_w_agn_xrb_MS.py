"""
Stacks data around positions and redshifts of GAMA galaxies


python cube_plot_profile_and_spectrum.py  CEN_bins_-14.0_ssfr_-11.0_10.0_mass_11.0 M1
python cube_plot_profile_and_spectrum.py  CEN_bins_-11.0_ssfr_-8.0_10.4_mass_11.0 M1


"""
sys.path.append( os.path.join(os.environ['GIT_STMOD'], 'src') )
from io import *

Z_SEL = sys.argv[1]

fig_dir = os.path.join( os.environ['GIT_STACK'], 'figures', 'Ti20_SDSS_stacked_galaxy_profile' )
os.system( 'mkdir -p ' + fig_dir )
prof_dir = os.path.join( os.environ['GIT_STACK'], 'figures', 'Ti20_SDSS_stacked_galaxy_profile', 'profiles_decomp_AXH_Mstar' )
os.system( 'mkdir -p ' + prof_dir )

p2_summary_file = os.path.join(prof_dir, 'FULL_SUMMARY_Mstarbin_centrals_'+Z_SEL+'.fits')

ALLT = []
M_bins = n.hstack(( n.arange(9, 10, 0.5),  n.arange(10, 12.1, 0.25), 13. ))
t_out = Table()
t_out['Ms0'] = M_bins[:-1]
t_out['Ms1'] = M_bins[1:]
for M_val, M_up in zip(t_out['Ms0'], t_out['Ms1']):
        print('='*100)
        print(M_val, M_up)
        arr = [M_val, M_up]
        try:
                p2_prof = glob.glob(os.path.join(mergedCube_dir,'CEN-ANY-0.01_'+Z_SEL+'_*-'+str(M_val)+'_Mstar_'+str(M_up)+'_STACKEDprofiles.fits' ))[0]
                print(p2_prof)
                CEN_ANY = get_profile_SDSS_Ti21(p2_prof)
        except(FileNotFoundError, IndexError):
                print('missing file')
                continue

        p2_fig = os.path.join(prof_dir, 'CEN_ANY_AGN-XRB-hotGAS_'+Z_SEL+'_Mstar_'+str(M_val)+'.png')
        p.figure(32, (5, 4.5) )
        # ALL
        str_NG = ', N='+str(CEN_ANY['N_gal'][0])
        str_MS = ', M$^*$='+str(n.round(CEN_ANY['MS_mean'][0],1))+'$\pm$'+str(n.round(CEN_ANY['MS_std'][0],1))
        N_all = CEN_ANY['N_gal'][0]
        p.axvline(CEN_ANY['R500c'][0], ls='--', lw=2, c='k')#, label=r'$\bar{R}_{500c}$')#+str(int(CEN_ANY['R500c'][0]))+'kpc')
        p.axhline(CEN_ANY['bg'][0]/300., ls='--', lw=2, c='grey')#, label='BG/300')
        p.errorbar( CEN_ANY['dd_xb'], CEN_ANY['dd_profile_BGSUB'],
                yerr = [ CEN_ANY['dd_profile_BGSUB'] - CEN_ANY['dd_profile_lo_BGSUB'], CEN_ANY['dd_profile_up_BGSUB'] - CEN_ANY['dd_profile_BGSUB'] ],
                xerr = [ CEN_ANY['dd_xb'] - CEN_ANY['x_lo'], CEN_ANY['x_up'] - CEN_ANY['dd_xb']],
                lw=2, ls='', color=cs ['CEN_ANY'], label='CEN ANY')# + str_NG + str_MS )
        # PSF
        p.step(CEN_ANY['dd_xb'], CEN_ANY['ps_profile_normed_BGSUB'], lw=0.8 , color=cs ['CEN_ANY'], where='mid')#, label='PSF' )
        # AGN
        str_AGN = r', L$_X$='+str(n.round(CEN_ANY['AGN_Co19_mean_LX_all'][0],1))+r', $f_{AGN}=$'+str(n.round(CEN_ANY['AGN_Co19_frac'][0]*100,1))+'%'
        p.step(CEN_ANY['dd_xb'], CEN_ANY['AGN_profile_normed_BGSUB'], lw=2, ls='dashed' , color=cs['AGN'], where='mid', label='AGN Co19')# + str_AGN)
        # XRB
        str_XRB = r', L$_X$='+str(n.round(CEN_ANY['XRB_Ai17_med_LX_all'][0],1))
        p.step(CEN_ANY['dd_xb'], CEN_ANY['XRB_profile_normed_BGSUB'], lw=2, ls='dashed' , color=cs['XRB'], where='mid', label='XRB Ai17')# + str_XRB )
        # hot gas
        Max_values = CEN_ANY['dd_profile_up_BGSUB'] - ( CEN_ANY['AGN_profile_normed_BGSUB']*0.75 + CEN_ANY['XRB_profile_normed_BGSUB']*0.75 )
        Min_values = CEN_ANY['dd_profile_lo_BGSUB'] - ( CEN_ANY['AGN_profile_normed_BGSUB']*1.25 + CEN_ANY['XRB_profile_normed_BGSUB']*1.25 )
        Mean_value = CEN_ANY['dd_profile_BGSUB'] - ( CEN_ANY['AGN_profile_normed_BGSUB'] + CEN_ANY['XRB_profile_normed_BGSUB'] )
        Max_values[Max_values<=1]=1
        Min_values[Min_values<=1]=1
        Mean_value[Mean_value<=1]=1
        area_shell = ( CEN_ANY['x_up']**2 -  CEN_ANY['x_lo']**2 ) * n.pi
        # central bin
        LX_mean    = interp1d( CEN_ANY['x_up'], n.cumsum(Mean_value*area_shell) ) ( 20. )
        LX_mean_up = interp1d( CEN_ANY['x_up'], n.cumsum(Max_values*area_shell) ) ( 20. )
        LX_mean_lo = interp1d( CEN_ANY['x_up'], n.cumsum(Min_values*area_shell) ) ( 20. )
        CEN_ANY['hotgas_LX_20kpc']    = LX_mean
        CEN_ANY['hotgas_LX_20kpc_up'] = LX_mean_up
        CEN_ANY['hotgas_LX_20kpc_lo'] = LX_mean_lo
        DLX =  0.5*(LX_mean_up-LX_mean_lo)
        print('central bin', '{:.2e}'.format(LX_mean), '{:.2e}'.format(LX_mean_lo), '{:.2e}'.format(LX_mean_up))#, DLX)
        # Rvir
        LX_mean    = interp1d( CEN_ANY['x_up'], n.cumsum(Mean_value*area_shell) ) ( CEN_ANY['R500c'][0] )
        LX_mean_up = interp1d( CEN_ANY['x_up'], n.cumsum(Max_values*area_shell) ) ( CEN_ANY['R500c'][0] )
        LX_mean_lo = interp1d( CEN_ANY['x_up'], n.cumsum(Min_values*area_shell) ) ( CEN_ANY['R500c'][0] )
        CEN_ANY['hotgas_LX_R500c']    = LX_mean
        CEN_ANY['hotgas_LX_R500c_up'] = LX_mean_up
        CEN_ANY['hotgas_LX_R500c_lo'] = LX_mean_lo
        DLX =  0.5*(LX_mean_up-LX_mean_lo)
        print('within RVIR', '{:.2e}'.format(LX_mean), '{:.2e}'.format(LX_mean_lo), '{:.2e}'.format(LX_mean_up))#, DLX)
        #SN = LX_mean/DLX
        #print('S/N=',SN)
        log10_LX_mean    = n.log10(LX_mean   )
        log10_LX_mean_up = n.log10(LX_mean_up)
        log10_LX_mean_lo = n.log10(LX_mean_lo)
        log10_DLX        =  0.5*(log10_LX_mean_up-log10_LX_mean_lo)
        p.errorbar( CEN_ANY['dd_xb'], Mean_value,
                yerr = [ Mean_value - Min_values, Max_values - Mean_value ],
                xerr = [ CEN_ANY['dd_xb'] - CEN_ANY['x_lo'], CEN_ANY['x_up'] - CEN_ANY['dd_xb']],
                lw=1.5, ls='', color=cs['hotGAS'] , label='Hot gas')#, L$_X$='+str(n.round(log10_LX_mean,2))+r'$\pm$' +str(n.round(log10_DLX,2)) )
        p.xlabel(r'$R_p$ [kpc]')
        p.ylabel(r'$S_X$ $[erg\; s^{-1}\; kpc^{-2}]$')
        p.legend(fontsize=12, loc=3, ncol=2)
        p.xlim((4, 2000))
        p.ylim((CEN_ANY['bg'][0]/5000., n.max(CEN_ANY['dd_profile_BGSUB'])*1.5))
        p.yscale('log')
        p.xscale('log')
        p.title(str(M_val)+'$<\log_{10}(M^*/M_\odot)<$'+str(M_up))
        p.tight_layout()
        p.savefig( p2_fig )
        p.clf()
        print(p2_fig, 'written')
        #cumulative LX integral interpolation
        Max_values = CEN_ANY['dd_profile_up_BGSUB']
        Min_values = CEN_ANY['dd_profile_lo_BGSUB']
        Mean_value = CEN_ANY['dd_profile_BGSUB']
        Max_values[Max_values<=1]=1
        Min_values[Min_values<=1]=1
        Mean_value[Mean_value<=1]=1
        area_shell = ( CEN_ANY['x_up']**2 -  CEN_ANY['x_lo']**2 ) * n.pi
        # central bin
        #r_vir_sel = (CEN_ANY['dd_xb']<CEN_ANY['R500c'][0])
        LX_mean    = interp1d( CEN_ANY['x_up'], n.cumsum(Mean_value*area_shell) ) ( CEN_ANY['R500c'][0] )
        LX_mean_up = interp1d( CEN_ANY['x_up'], n.cumsum(Max_values*area_shell) ) ( CEN_ANY['R500c'][0] )
        LX_mean_lo = interp1d( CEN_ANY['x_up'], n.cumsum(Min_values*area_shell) ) ( CEN_ANY['R500c'][0] )
        CEN_ANY['total_LX_R500c']    = LX_mean
        CEN_ANY['total_LX_R500c_up'] = LX_mean_up
        CEN_ANY['total_LX_R500c_lo'] = LX_mean_lo
        ALLT.append(CEN_ANY)

MERGE = Table(vstack((ALLT)))
MERGE.write(p2_summary_file, overwrite = True)


###t = Table()
###t['M0'] = LX_outs[0]
###t['M1'] = LX_outs[1]
###t['Ng'] = LX_outs[2].astype('int')
###t['R500c'] = n.round(LX_outs[3],1)
###t['AGN_LX'] = n.round(LX_outs[6],2)
###t['AGN_frac'] = n.round(LX_outs[7],1)
###t['XRB_LX'] = n.round(LX_outs[8],2)
####t['hotgas_1stBin_LX_mean'] = n.round(LX_outs[9]/1e39,2)
####t['hotgas_1stBin_LX_mean_up'] = n.round(LX_outs[10]/1e39,2)
####t['hotgas_1stBin_LX_mean_lo'] = n.round(LX_outs[11]/1e39,2)
###t['hotgas_LX_mean'] = n.round(LX_outs[12]/1e39,2)
###t['hotgas_LX_mean_up'] = n.round(LX_outs[13]/1e39,2)
###t['hotgas_LX_mean_lo'] = n.round(LX_outs[14]/1e39,2)
###t.write(os.path.join(prof_dir, 'LX_summary_Mstarbin_centrals_'+Z_SEL+'.latex'), format='latex', overwrite = True)



