sys.path.append( os.path.join(os.environ['GIT_STMOD'], 'src') )
from io import *

fig_dir = os.path.join(os.environ['GIT_STACK'], 'figures', 'Ti20_SDSS_stacked_galaxy_profile', 'scaling_relations')
os.system( 'mkdir -p ' + fig_dir )

Z_SEL = 'Z3P'

prof_dir = os.path.join( os.environ['GIT_STACK'], 'figures', 'Ti20_SDSS_stacked_galaxy_profile', 'profiles_decomp_AXH' )
MERGE_Mh = Table.read(os.path.join(prof_dir, 'FULL_SUMMARY_Mhalobin_centrals_'+Z_SEL+'.fits'))
t_mh = MERGE_Mh[n.unique(MERGE_Mh['file_ID'], return_index=True)[1]]
t_mh = t_mh [(t_mh['N_gal'] > 1000)&(t_mh['M_min'] >= 12.0)&(t_mh['M_max'] <=14.0)]

prof_dir = os.path.join( os.environ['GIT_STACK'], 'figures', 'Ti20_SDSS_stacked_galaxy_profile', 'profiles_decomp_AXH_Mstar' )
MERGE_Ms = Table.read(os.path.join(prof_dir, 'FULL_SUMMARY_Mstarbin_centrals_'+Z_SEL+'.fits'))
t_ms = MERGE_Ms[n.unique(MERGE_Ms['file_ID'], return_index=True)[1]]
t_ms = t_ms [(t_ms['N_gal'] > 1000)&(t_ms['M_min'] >= 10.75)&(t_ms['M_max'] <=11.75)]

MS_A15, LX_max_A15, LX_min_A15 = n.loadtxt(  os.path.join( os.environ['GIT_EVTXCORR'], 'data', 'anderson2015', 'anderson2015_fig5_total.ascii'), unpack = True)
MS_A15 = n.arange(10.1,11.95, 0.1)+0.05


p2_fig = os.path.join(fig_dir, 'LX_Ms_CEN_ANY.png')
p.figure(5, (5,4.5) )
p.errorbar( t_ms['GAL_stellarMass_mean'],
            t_ms['total_LX_R500c'],
            yerr = [ t_ms['total_LX_R500c'] - t_ms['total_LX_R500c_lo'], t_ms['total_LX_R500c_up'] - t_ms['total_LX_R500c'] ],
            xerr = t_ms['GAL_stellarMass_std'],
            lw=2, ls='', color=cs ['CEN_ANY_MS'], label='M*, total' )

p.errorbar( t_ms['GAL_stellarMass_mean'],
            t_ms['hotgas_LX_R500c'],
            yerr = [ t_ms['hotgas_LX_R500c'] - t_ms['hotgas_LX_R500c_lo'], t_ms['hotgas_LX_R500c_up'] - t_ms['hotgas_LX_R500c'] ],
            xerr = t_ms['GAL_stellarMass_std'],
            lw=2, ls='', color=cs['hotGAS_MS'], label='M*, hot gas' )

p.errorbar( t_mh['GAL_stellarMass_mean'],
            t_mh['total_LX_R500c'],
            yerr = [ t_mh['total_LX_R500c'] - t_mh['total_LX_R500c_lo'], t_mh['total_LX_R500c_up'] - t_mh['total_LX_R500c'] ],
            xerr = t_mh['GAL_stellarMass_std'],
            lw=2, ls='', color=cs ['CEN_ANY'], label='Mh, total' )

p.errorbar( t_mh['GAL_stellarMass_mean'],
            t_mh['hotgas_LX_R500c'],
            yerr = [ t_mh['hotgas_LX_R500c'] - t_mh['hotgas_LX_R500c_lo'], t_mh['hotgas_LX_R500c_up'] - t_mh['hotgas_LX_R500c'] ],
            xerr = t_mh['GAL_stellarMass_std'],
            lw=2, ls='', color=cs['hotGAS'], label='Mh, hot gas' )

#y1 = 10** ( (LX_min_A15+LX_max_A15)/2. )
#yerr_up = 10** LX_max_A15 - y1
#yerr_lo = - 10** LX_min_A15 + y1
#p.errorbar( MS_A15-n.log10(0.7), y1 , xerr=0.05, yerr=  [yerr_lo, yerr_up],  label=r'Anderson 2015', ls='', fmt='s',markersize=8,markerfacecolor='none', color='grey')
p.plot([10.75, 10.75], [0.9e39, 3e39], 'k--', alpha=0.8)
p.plot([11., 11]     , [0.9e39, 3e39], 'k--', alpha=0.8)
p.plot([11.25, 11.25], [0.9e39, 3e39], 'k--', alpha=0.8)
p.plot([11.75, 11.75], [0.9e39, 3e39], 'k--', alpha=0.8)
p.text(10.8,  1e39, 'MW' , alpha=0.8)
p.text(11.02, 1e39, 'M31', alpha=0.8)
p.text(11.35, 1e39, 'Groups', alpha=0.8)
p.xlabel(r'$\log_{10}(M_*)$ [M$_\odot$]')
p.ylabel(r'L$^{<R_{500c}}_{X}$ [erg/s]')
p.legend(fontsize=12, loc='upper left', ncol=1)#, title='Centrals')
p.yscale('log')
p.xlim((10.5, 11.8))
p.ylim((0.9e39, 3e42))
p.tight_layout()
p.savefig( p2_fig )
p.clf()
print(p2_fig, 'written')



p2_fig = os.path.join(fig_dir, 'LX_Mhsel_CEN.png')
p.figure(5, (5,4.5) )
p.errorbar( t_ms['GAL_Mhalo_mean'],
            t_ms['total_LX_R500c'],
            yerr = [ t_ms['total_LX_R500c'] - t_ms['total_LX_R500c_lo'], t_ms['total_LX_R500c_up'] - t_ms['total_LX_R500c'] ],
            xerr = t_ms['GAL_Mhalo_std'],
            lw=2, ls='', color=cs ['CEN_ANY_MS'], label='M*, total' )
p.errorbar( t_ms['GAL_Mhalo_mean'],
            t_ms['hotgas_LX_R500c'],
            yerr = [ t_ms['hotgas_LX_R500c'] - t_ms['hotgas_LX_R500c_lo'], t_ms['hotgas_LX_R500c_up'] - t_ms['hotgas_LX_R500c'] ],
            xerr = t_ms['GAL_Mhalo_std'],
            lw=2, ls='', color=cs['hotGAS_MS'], label='M*, hot gas' )

p.errorbar( t_mh['GAL_Mhalo_mean'],
            t_mh['total_LX_R500c'],
            yerr = [ t_mh['total_LX_R500c'] - t_mh['total_LX_R500c_lo'], t_mh['total_LX_R500c_up'] - t_mh['total_LX_R500c'] ],
            xerr = t_mh['GAL_Mhalo_std'],
            lw=2, ls='', color=cs ['CEN_ANY'], label='Mh, total' )
p.errorbar( t_mh['GAL_Mhalo_mean'],
            t_mh['hotgas_LX_R500c'],
            yerr = [ t_mh['hotgas_LX_R500c'] - t_mh['hotgas_LX_R500c_lo'], t_mh['hotgas_LX_R500c_up'] - t_mh['hotgas_LX_R500c'] ],
            xerr = t_mh['GAL_Mhalo_std'],
            lw=2, ls='', color=cs['hotGAS'], label='Mh, hot gas' )

p.plot([12, 12]    , [0.9e39, 3e39], 'k--', alpha=0.8)
p.plot([12.5, 12.5], [0.9e39, 3e39], 'k--', alpha=0.8)
p.plot([13.0, 13.0], [0.9e39, 3e39], 'k--', alpha=0.8)
p.plot([14.0, 14.0], [0.9e39, 3e39], 'k--', alpha=0.8)
p.text(12.1, 1e39, 'MW' , alpha=0.8)
p.text(12.6, 1e39, 'M31', alpha=0.8)
p.text(13.3, 1e39, 'Groups', alpha=0.8)
p.xlabel(r'M$_{200m}$ [M$_\odot$]')
p.ylabel(r'L$^{<R_{500c}}_{X}$ [erg/s]')
p.legend(fontsize=12, loc='upper left', ncol=1)#, title='Centrals')
p.yscale('log')
p.xlim((11.8, 14.1))
p.ylim((0.9e39, 3e42))
p.tight_layout()
p.savefig( p2_fig )
p.clf()
print(p2_fig, 'written')


sys.exit()



prof_dir = os.path.join( os.environ['GIT_STACK'], 'figures', 'Ti20_SDSS_stacked_galaxy_profile', 'profiles_decomp' )
MERGE_Mhbcrs = Table.read(os.path.join(prof_dir, 'FULL_SUMMARY_Mhalobin_centrals_BC_RS_'+Z_SEL+'.fits'))
t_rb = MERGE_Mhbcrs[n.unique(MERGE_Mhbcrs['file_ID'], return_index=True)[1]]

t_rs = t_rb[(t_rb['SFR']=="RS")&(t_rb['N_gal'] > 1000)]
t_bc = t_rb[(t_rb['SFR']=="BC")&(t_rb['N_gal'] > 1000)]


prof_dir = os.path.join( os.environ['GIT_STACK'], 'figures', 'Ti20_SDSS_stacked_galaxy_profile', 'profiles_decomp_Mstar' )
MERGE_Mhbcrs = Table.read(os.path.join(prof_dir, 'FULL_SUMMARY_Mstarbin_centrals_BC_RS_'+Z_SEL+'.fits'))
t_rb = MERGE_Mhbcrs[n.unique(MERGE_Mhbcrs['file_ID'], return_index=True)[1]]

t_Mrs = t_rb[(t_rb['SFR']=="RS")&(t_rb['N_gal'] > 1000)]
t_Mbc = t_rb[(t_rb['SFR']=="BC")&(t_rb['N_gal'] > 1000)]



p2_fig = os.path.join(fig_dir, 'LX_Ms_CEN_ANY+RS+BSsplitMS.png')
p.figure(5, (5,4.5) )
p.errorbar( t_ms['GAL_stellarMass_mean'],
            t_ms['total_LX_R500c'],
            yerr = [ t_ms['total_LX_R500c'] - t_ms['total_LX_R500c_lo'], t_ms['total_LX_R500c_up'] - t_ms['total_LX_R500c'] ],
            xerr = t_ms['GAL_stellarMass_std'],
            lw=2, ls='', color='darkgreen', label='All' )

p.errorbar( t_Mrs['GAL_stellarMass_mean'],
            t_Mrs['total_LX_R500c'],
            yerr = [ t_Mrs['total_LX_R500c'] - t_Mrs['total_LX_R500c_lo'], t_Mrs['total_LX_R500c_up'] - t_Mrs['total_LX_R500c'] ],
            xerr = t_Mrs['GAL_stellarMass_std'],
            lw=2, ls='', color='red', label='RS' )

p.errorbar( t_Mbc['GAL_stellarMass_mean'],
            t_Mbc['total_LX_R500c'],
            yerr = [ t_Mbc['total_LX_R500c'] - t_Mbc['total_LX_R500c_lo'], t_Mbc['total_LX_R500c_up'] - t_Mbc['total_LX_R500c'] ],
            xerr = t_Mbc['GAL_stellarMass_std'],
            lw=2, ls='', color='blue', label='BC' )

y1 = 10** ( (LX_min_A15+LX_max_A15)/2. )
yerr_up = 10** LX_max_A15 - y1
yerr_lo = - 10** LX_min_A15 + y1
p.errorbar( MS_A15, y1 , xerr=0.05, yerr=  [yerr_lo, yerr_up],  label=r'Anderson 2015', ls='', fmt='s',markersize=8,markerfacecolor='none', color='grey')

p.xlabel(r'$\log_{10}(M_*)$ [M$_\odot$]')
p.ylabel(r'L$^{<R_{500c}}_{X}$ [erg/s]')
p.legend(fontsize=12, loc=2, ncol=1, title='Centrals')
p.yscale('log')
p.xlim((9.0, 12.2))
p.ylim((1e39, 1e44))
p.title('Stellar mass selection')
p.tight_layout()
p.savefig( p2_fig )
p.clf()
print(p2_fig, 'written')


p2_fig = os.path.join(fig_dir, 'LX_Ms_CEN_ANY+RS+BSsplit.png')
p.figure(5, (5,4.5) )

p.errorbar( t_mh['GAL_stellarMass_mean'],
            t_mh['total_LX_R500c'],
            yerr = [ t_mh['total_LX_R500c'] - t_mh['total_LX_R500c_lo'], t_mh['total_LX_R500c_up'] - t_mh['total_LX_R500c'] ],
            xerr = t_mh['GAL_stellarMass_std'],
            lw=2, ls='', color='orange', label='All' )

p.errorbar( t_rs['GAL_stellarMass_mean'],
            t_rs['total_LX_R500c'],
            yerr = [ t_rs['total_LX_R500c'] - t_rs['total_LX_R500c_lo'], t_rs['total_LX_R500c_up'] - t_rs['total_LX_R500c'] ],
            xerr = t_rs['GAL_stellarMass_std'],
            lw=2, ls='', color='red', label='RS' )

p.errorbar( t_bc['GAL_stellarMass_mean'],
            t_bc['total_LX_R500c'],
            yerr = [ t_bc['total_LX_R500c'] - t_bc['total_LX_R500c_lo'], t_bc['total_LX_R500c_up'] - t_bc['total_LX_R500c'] ],
            xerr = t_bc['GAL_stellarMass_std'],
            lw=2, ls='', color='blue', label='BC' )


y1 = 10** ( (LX_min_A15+LX_max_A15)/2. )
yerr_up = 10** LX_max_A15 - y1
yerr_lo = - 10** LX_min_A15 + y1
p.errorbar( MS_A15, y1 , xerr=0.05, yerr=  [yerr_lo, yerr_up],  label=r'Anderson 2015', ls='', fmt='s',markersize=8,markerfacecolor='none', color='grey')

p.xlabel(r'$\log_{10}(M_*)$ [M$_\odot$]')
p.ylabel(r'L$^{<R_{500c}}_{X}$ [erg/s]')
p.legend(fontsize=12, loc=2, ncol=1, title='Centrals')
p.yscale('log')
p.xlim((9.0, 12.2))
p.ylim((1e39, 1e44))
p.title('Halo mass selection')
p.tight_layout()
p.savefig( p2_fig )
p.clf()
print(p2_fig, 'written')





p2_fig = os.path.join(fig_dir, 'LX_Mhsel_CEN_ANY+RS+BSsplit.png')
p.figure(5, (5,4.5) )

p.errorbar( t_mh['GAL_Mhalo_mean'],
            t_mh['total_LX_R500c'],
            yerr = [ t_mh['total_LX_R500c'] - t_mh['total_LX_R500c_lo'], t_mh['total_LX_R500c_up'] - t_mh['total_LX_R500c'] ],
            xerr = t_mh['GAL_Mhalo_std'],
            lw=2, ls='', color='orange', label='All' )

p.errorbar( t_rs['GAL_Mhalo_mean'],
            t_rs['total_LX_R500c'],
            yerr = [ t_rs['total_LX_R500c'] - t_rs['total_LX_R500c_lo'], t_rs['total_LX_R500c_up'] - t_rs['total_LX_R500c'] ],
            xerr = t_rs['GAL_Mhalo_std'],
            lw=2, ls='', color='red', label='RS' )

p.errorbar( t_bc['GAL_Mhalo_mean'],
            t_bc['total_LX_R500c'],
            yerr = [ t_bc['total_LX_R500c'] - t_bc['total_LX_R500c_lo'], t_bc['total_LX_R500c_up'] - t_bc['total_LX_R500c'] ],
            xerr = t_bc['GAL_Mhalo_std'],
            lw=2, ls='', color='blue', label='BC' )

p.xlabel(r'M$_{200m}$ [M$_\odot$]')
p.ylabel(r'L$^{<R_{500c}}_{X}$ [erg/s]')
p.legend(fontsize=14, loc=0, ncol=1, title='Centrals')
p.yscale('log')
p.xlim((11.0, 14.5))
p.ylim((1e39, 1e44))
p.title('Halo mass selection')
p.tight_layout()
p.savefig( p2_fig )
p.clf()
print(p2_fig, 'written')



p2_fig = os.path.join(fig_dir, 'LX_Mhsel_CEN_ANY+RS+BSsplitMS.png')
p.figure(5, (5,4.5) )

p.errorbar( t_ms['GAL_Mhalo_mean'],
            t_ms['total_LX_R500c'],
            yerr = [ t_ms['total_LX_R500c'] - t_ms['total_LX_R500c_lo'], t_ms['total_LX_R500c_up'] - t_ms['total_LX_R500c'] ],
            xerr = t_ms['GAL_Mhalo_std'],
            lw=2, ls='', color='orange', label='All' )

p.errorbar( t_Mrs['GAL_Mhalo_mean'],
            t_Mrs['total_LX_R500c'],
            yerr = [ t_Mrs['total_LX_R500c'] - t_Mrs['total_LX_R500c_lo'], t_Mrs['total_LX_R500c_up'] - t_Mrs['total_LX_R500c'] ],
            xerr = t_Mrs['GAL_Mhalo_std'],
            lw=2, ls='', color='red', label='RS' )

p.errorbar( t_Mbc['GAL_Mhalo_mean'],
            t_Mbc['total_LX_R500c'],
            yerr = [ t_Mbc['total_LX_R500c'] - t_Mbc['total_LX_R500c_lo'], t_Mbc['total_LX_R500c_up'] - t_Mbc['total_LX_R500c'] ],
            xerr = t_Mbc['GAL_Mhalo_std'],
            lw=2, ls='', color='blue', label='BC' )

p.xlabel(r'M$_{200m}$ [M$_\odot$]')
p.ylabel(r'L$^{<R_{500c}}_{X}$ [erg/s]')
p.legend(fontsize=14, loc=0, ncol=1, title='Centrals')
p.yscale('log')
p.xlim((11.0, 14.5))
p.ylim((1e39, 1e44))
p.title('Stellar mass selection')
p.tight_layout()
p.savefig( p2_fig )
p.clf()
print(p2_fig, 'written')



