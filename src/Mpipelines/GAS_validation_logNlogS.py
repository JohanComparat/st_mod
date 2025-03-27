"""
Validation of the gas model
Figures :
 * logNlogS (soft band)

python GAS_validation_logNlogS.py LC0002
python GAS_validation_logNlogS.py LC0060
python GAS_validation_logNlogS.py LC1800 12
python GAS_validation_logNlogS.py LC1800 17
python GAS_validation_logNlogS.py LC1800 22
python GAS_validation_logNlogS.py LC1800 27
python GAS_validation_logNlogS.py FullSky 5
python GAS_validation_logNlogS.py FullSky 6
python GAS_validation_logNlogS.py FullSky 7
python GAS_validation_logNlogS.py FullSky 8
"""
import time
t0 = time.time()
import sys, os, glob
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as p
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

from scipy.interpolate import interp1d
from scipy.interpolate import RegularGridInterpolator

from scipy.stats import scoreatpercentile

from astropy.table import Table, vstack
import astropy.io.fits as fits

import numpy as np

cs = {}
cs ['ALL_ANY'] = 'black'
cs ['CEN_ANY'] = 'darkgreen'
cs ['PS_CEN_ANY'] = 'red'
cs ['CEN_RS']  = 'firebrick'
cs ['CEN_BC']  = 'navy'
cs ['CEN_RS_hotGAS']  = 'goldenrod'
cs ['CEN_BC_hotGAS']  = 'aqua'
cs ['SAT_ANY'] = 'aquamarine'
cs ['SAT_BC']  = 'cornflowerblue'
cs ['SAT_RS']  = 'lightcoral'

cs ['ALL_Mstar']  = 'deeppink'
cs ['CEN_Mstar']  = 'fuchsia'
cs ['SAT_Mstar']  = 'plum'
cs ['ALL_Mhalo']  = 'forestgreen'
cs ['CEN_Mhalo']  = 'chartreuse'
cs ['SAT_Mhalo']  = 'turquoise'

cs['AGN'] = 'darkmagenta'
cs['XRB'] = 'hotpink'
cs['hotGAS'] = 'darkorange'
cs ['CEN_ANY_MS'] = 'lime'
cs['hotGAS_MS'] = 'moccasin'


cosmoUNIT = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
h = 0.6774
L_box = 2000.0 / h
cosmo = cosmoUNIT

nl = lambda sel : len(sel.nonzero()[0])
zs = np.arange(0.0000001, 7.1, 0.001)
dm_itp = interp1d(zs, cosmo.distmod(zs).value)

# K-correction to get observed fluxes
#nh_vals = 10**n.arange(-2,4+0.01,0.5)#0.05)
#kT_vals = 10**n.arange(-2.09,1.8,0.01)
#z_vals = n.hstack((n.array([0.]), n.arange(0.001, 0.01, 0.001), 10**n.arange(n.log10(0.01), n.log10(6.1), 0.01)))
#dir_2_result = os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/model_GAS', 'xray_k_correction')
#itp_z, itp_kt, itp_frac_obs = np.loadtxt( os.path.join( dir_2_result, "fraction_observed_clusters_no_nH.txt"), unpack=True )
#YY_z, ZZ_kt = np.meshgrid(z_vals, kT_vals)
#shape_i = YY_z.shape
#matrix_2d = itp_frac_obs.reshape(shape_i)
#from scipy.interpolate import RegularGridInterpolator
#attenuation_2d = RegularGridInterpolator((kT_vals, z_vals), matrix_2d)

LC_dir= sys.argv[1] #'LC0002'
validation_dir       = os.path.join(os.environ['GIT_STMOD_DATA'], 'data', 'validation','validation_GAS')
validation_dir_lNlS = os.path.join(validation_dir, 'logNlogS' )
os.system('mkdir -p ' + validation_dir       )
os.system('mkdir -p ' + validation_dir_lNlS )


z_dir_all = np.array([  'z0p00',
                    'z0p02',
                    'z0p05',
                    'z0p09',
                    'z0p14',
                    'z0p19',
                    'z0p25',
                    'z0p30',
                    'z0p36',
                    'z0p43',
                    'z0p49',
                    'z0p56',
                    'z0p63',
                    'z0p70',
                    'z0p78',
                    'z0p86',
                    'z0p94',
                    'z1p03',
                    'z1p12',
                    'z1p22',
                    'z1p32',
                    'z1p43',
                    'z1p54',
                    'z1p65',
                    'z1p77',
                    'z1p90',
                    'z2p03',
                    'z2p17',
                    'z2p31',
                    'z2p46',
                    'z2p62',
                    'z2p78',
                    'z2p95',
                    'z3p13',
                    'z3p32',
                    'z3p61',
                    'z3p93',
                    'z4p27',
                    'z4p63',
                    'z5p15',
                    'z5p73' ])

if LC_dir == 'LC0002' :
    ra0 = 149.2
    ra1 = 151.1
    de0 = 1.3
    de1 = 3.1
    z_dirs = z_dir_all
    area = (de1-de0)*(ra1-ra0)
    z_max = z_dirs[-1]
elif LC_dir == 'LC0060' :
    ra0 = 129
    ra1 = 141
    de0 = -2
    de1 = 3
    z_dirs = z_dir_all[:37]
    area = (de1-de0)*(ra1-ra0)
    z_max = z_dirs[-1]
elif LC_dir == 'LC1800' :
    ra0 = 0.0
    ra1 = 360.
    de0 = -2.5
    de1 = 2.5
    i_zmax = int(sys.argv[2])
    z_dirs = z_dir_all[:i_zmax]
    #z_dirs = z_dir_all[:27]
    area = (de1-de0)*(ra1-ra0)
    z_max = z_dirs[-1]
elif LC_dir == 'FullSky' :
    ra0 = 0.0
    ra1 = 360.
    de0 = -2.5
    de1 = 2.5
    i_zmax = int(sys.argv[2])
    z_dirs = z_dir_all[:i_zmax]
    area = 129600/np.pi
    z_max = z_dirs[-1]


#p_2_catal_GAS_b06 = np.array( glob.glob(os.path.join(os.environ['UCHUU'], LC_dir, 'z?p??', 'replication_*', 'Xgas_bHS0.6.fits')))
#p_2_catal_GAS_b06.sort()
#within_zmax_b06 = np.isin( np.array([el.split('/')[-3] for el in p_2_catal_GAS_b06 ]), z_dirs )
#p_2_catal_GAS_b06 = p_2_catal_GAS_b06[within_zmax_b06]

p_2_catal_GAS_b08 = np.array( glob.glob(os.path.join(os.environ['UCHUU'], LC_dir, 'z?p??', 'replication_*', 'Xgas_bHS0.8.fits')))
p_2_catal_GAS_b08.sort()
within_zmax_b08 = np.isin( np.array([el.split('/')[-3] for el in p_2_catal_GAS_b08 ]), z_dirs )
p_2_catal_GAS_b08 = p_2_catal_GAS_b08[within_zmax_b08]

#p_2_catal_GAS_b10 = np.array( glob.glob(os.path.join(os.environ['UCHUU'], LC_dir, 'z?p??', 'replication_*', 'Xgas_bHS1.0.fits')))
#p_2_catal_GAS_b10.sort()
#within_zmax_b10 = np.isin( np.array([el.split('/')[-3] for el in p_2_catal_GAS_b10 ]), z_dirs )
#p_2_catal_GAS_b10 = p_2_catal_GAS_b10[within_zmax_b10]

def get_FX(p_2_cata):
    print(p_2_cata)
    XGA = Table.read(p_2_cata)
    print(len(XGA))
    #print(XGA['CLUSTER_kT'].min(), XGA['redshift_S'].min())
    #print(XGA['CLUSTER_kT'].max(), XGA['redshift_S'].max())
    #ZZ = XGA['redshift_S']
    #ZZ[ZZ>=3.98] = 3.98
    # kT_vals = 10**n.arange(-2.09,1.8,0.01)
    #TT = XGA['CLUSTER_kT']
    #TT[TT<-2.09]=-2.09
    #TT[TT>1.8]=1.8
    #k_correction_2d = attenuation_2d( np.transpose([10**TT, ZZ]))
    #DL_z = (cosmo.luminosity_distance(XGA['redshift_S']).to(u.cm)).value
    #dl2_4pi = np.log10(4*np.pi*(DL_z)**2.)
    #LX_obsF_500c = XGA['CLUSTER_LX_soft_RF_R500c'] -np.log10( k_correction_2d ) - dl2_4pi
    #LX_obsF_200c = XGA['CLUSTER_LX_soft_RF_R200c'] -np.log10( k_correction_2d ) - dl2_4pi
    #LX_obsF_vir = XGA['CLUSTER_LX_soft_RF_Rvir']   -np.log10( k_correction_2d ) - dl2_4pi
    t_out = Table()
    t_out['M500c'] = np.log10(XGA['M500c'])
    t_out['z']= XGA['redshift_S']
    #t_out['FX_obsF_500']= XGA['CLUSTER_FX_soft_OBS_R500c_nHattenuated']
    t_out['FX_obsF_500']= XGA['CLUSTER_FX_soft_OBS_R500c']
    #t_out['FX_obsF_200']= LX_obsF_200c
    #t_out['FX_obsF_vir']= LX_obsF_vir
    return t_out

#DATA_b06 = get_FX(p_2_catal_GAS_b06[0])
#for p_2_cata in p_2_catal_GAS_b06[1:]:
    #DATA_b06 = vstack(( DATA_b06, get_FX(p_2_cata)))

DATA_b08 = get_FX(p_2_catal_GAS_b08[0])
for p_2_cata in p_2_catal_GAS_b08[1:]:
    DATA_b08 = vstack(( DATA_b08, get_FX(p_2_cata)))

#DATA_b10 = get_FX(p_2_catal_GAS_b10[0])
#for p_2_cata in p_2_catal_GAS_b10[1:]:
    #DATA_b10 = vstack(( DATA_b10, get_FX(p_2_cata)))

# # #
#
#
# Literature
#
#
# # #

#spiders = fits.open(os.path.join(os.environ['HOME'], 'data/spiders/cluster', 'mastercatalogue_FINAL_CODEXID.fits'))[1].data
#FX_bins = np.arange(8, 18., 0.25)
#out = np.cumsum(np.histogram(-np.log10(spiders['FLUX052']), bins=FX_bins)[0])
#c_out_spiders = out/5128.
#c_err_spiders = c_out_spiders * out**(-0.5)
#x_out_spiders = - 0.5 * (FX_bins[1:] + FX_bins[:-1])
#ok_spiders = (c_err_spiders>0)&(c_out_spiders>0)&(c_out_spiders>2*c_err_spiders)&(x_out_spiders>-14)


#p_2_header = os.path.join(os.environ['GIT_STMOD_DATA'],'logNlogS/spidersLogNLogS','lgnlgs_header.tbl')
#  (1)              (2)        (3)      (4)
#lg(flux) N>flux: 0.1<z<0.3 0.3<z<0.5 0.5<z<0.7
spiders_om205s860 = np.loadtxt( os.path.join(os.environ['GIT_STMOD_DATA'], 'data/validation/validation_GAS/logN_logS/spidersLogNLogS','lgnlgs_om205s860.tbl'), unpack = True )
spiders_om270s792 = np.loadtxt( os.path.join(os.environ['GIT_STMOD_DATA'], 'data/validation/validation_GAS/logN_logS/spidersLogNLogS','lgnlgs_om270s792.tbl'), unpack = True )
spiders_om342s730 = np.loadtxt( os.path.join(os.environ['GIT_STMOD_DATA'], 'data/validation/validation_GAS/logN_logS/spidersLogNLogS','lgnlgs_om342s730.tbl'), unpack = True )
spiders_om307s823_i = np.loadtxt( os.path.join(os.environ['GIT_STMOD_DATA'], 'data/validation/validation_GAS/logN_logS/spidersLogNLogS','lgnlgs_om307s823h67forecast_logz_con.tbl'), unpack = True )
spiders_om307s823 = spiders_om307s823_i.T[4:].T

#
# logN figure logS
#
fig_out = os.path.join(validation_dir_lNlS, LC_dir+'_zlt_'+z_max+'_logNlogS.png' )

p.figure(3, (6., 6.))
e1_F0223, e1_logNlogS = np.loadtxt(os.path.join(    os.environ["GIT_STMOD_DATA"], 'data/validation/validation_GAS/logN_logS', 'erass1.txt'), unpack=True)
efeds_F0223, efeds_logNlogS = np.loadtxt(os.path.join(    os.environ["GIT_STMOD_DATA"], 'data/validation/validation_GAS/logN_logS', 'efeds.txt'), unpack=True)
noras_F0223, noras_logNlogS = np.loadtxt(os.path.join(    os.environ["GIT_STMOD_DATA"], 'data/validation/validation_GAS/logN_logS', 'noras.txt'), unpack=True)
X_shift = np.log10(0.767)
p.fill_between(
    np.log10(e1_F0223)+X_shift,
    y1= e1_logNlogS*0.99,
    y2= e1_logNlogS*1.01,
    rasterized=True,
    alpha=0.5,
    label='eRASS1')
p.fill_between(
    np.log10(efeds_F0223)+X_shift,
    y1= efeds_logNlogS*0.95,
    y2= efeds_logNlogS*1.05,
    rasterized=True,
    alpha=0.5,
    label='eFEDS')
p.fill_between(
    np.log10(noras_F0223)+X_shift,
    y1= noras_logNlogS*0.9,
    y2= noras_logNlogS*1.1,
    rasterized=True,
    alpha=0.5,
    label='NORAS')

path_2_logNlogS_data = os.path.join(    os.environ["GIT_STMOD_DATA"], 'data/validation/validation_GAS/logN_logS', 'logNlogS_Finoguenov_cosmos_2007_clusters.data')
x_data, y_data, y_data_min, y_data_max = np.loadtxt(    path_2_logNlogS_data, unpack=True)
p.fill_between(
    np.log10(x_data),
    y1= y_data_min,
    y2= y_data_max,
    rasterized=True,
    alpha=0.2,
    label='Fi07 COSMOS')
path_2_logNlogS_data = os.path.join(    os.environ["GIT_STMOD_DATA"], 'data/validation/validation_GAS/logN_logS', 'logNlogS_Finoguenov_ecdfs_2015_clusters.data')
x_data, y_data, y_data_min, y_data_max = np.loadtxt(    path_2_logNlogS_data, unpack=True)
p.fill_between(
    np.log10(x_data),
    y1=y_data_min,
    y2=y_data_max,
    rasterized=True,
    alpha=0.2,
    label='Fi15 CDFS')
#path_2_logNlogS_data = os.path.join(os.environ["GIT_STMOD_DATA"],'data/validation/validation_GAS/logN_logS','Boehringer_noras2_2017.data')
#FX_01_24, N_per_str_up, N_per_str_low = np.loadtxt(    path_2_logNlogS_data, unpack=True)
#N_per_str = (N_per_str_up + N_per_str_low)/2.
#x_data = FX_01_24 * 1e-12 * 0.67 / 1
#y_data = N_per_str / (180/np.pi)**2
#y_data_up = N_per_str_up / (180/np.pi)**2
#y_data_lo = N_per_str_low / (180/np.pi)**2
#p.plot(np.log10(x_data)[3:-2], y_data[3:-2], rasterized=True, ls='dashed', label='B17 NORAS 2', lw=2)
#p.fill_between(np.log10(x_data), y1=y_data_lo, y2=y_data_up, rasterized=True, label='Bo17 NORAS 2', color='c')
#
#
# SPIDERS
#
#
#p.fill_between(x_out_spiders[ok_spiders], y1=c_out_spiders[ok_spiders]-c_err_spiders[ok_spiders], y2=c_out_spiders[ok_spiders]+c_err_spiders[ok_spiders], label='JC SPIDERS', alpha=0.3, color='m')
#
x_spiders = spiders_om205s860[0]
y1_spiders = np.min([np.sum(spiders_om205s860[1:,:] , axis=0), np.sum(spiders_om270s792[1:,:] , axis=0), np.sum(spiders_om342s730[1:,:] , axis=0), np.sum(spiders_om307s823[1:,:], axis=0)], axis=0)
y2_spiders = np.max([np.sum(spiders_om205s860[1:,:] , axis=0), np.sum(spiders_om270s792[1:,:] , axis=0), np.sum(spiders_om342s730[1:,:] , axis=0), np.sum(spiders_om307s823[1:,:], axis=0)], axis=0)
y_mean_count = (5100 * (y1_spiders + y2_spiders)/2.).astype('int')
y_err_count = y_mean_count**(-0.5)
#p.fill_between(x_spiders, y1 = y1_spiders*(1-y_err_count), y2 = y2_spiders*(1+y_err_count), label='Fi20 SPIDERS', alpha=0.3, color='c')
#
#p.plot(x_spiders,np.sum(spiders_om270s792[1:,:], axis=0), label='om270s792')
#y_spiders = np.sum(spiders_om307s823[1:,:], axis=0)
#y_mean_count = ( 5100 * y_spiders ).astype('int')
#y_err_count = y_mean_count**(-0.5) * 2
#p.plot(spiders_om307s823[0], y_spiders, label='SPIDERS om307s823')
p.fill_between(x_spiders, y1 = y1_spiders*(1-y_err_count), y2 = y2_spiders*(1+y_err_count), label='Fi20 SPIDERS', alpha=0.2)


#FX_bins = np.arange(8, 18., 0.25)
#out = np.cumsum(np.histogram(- DATA_b06['FX_obsF_500'], bins=FX_bins)[0])
#c_out_mock = out/area
#c_err_mock = c_out_mock * out**(-0.5)
#x_out_mock = - 0.5 * (FX_bins[1:] + FX_bins[:-1])
#ok_mock = (c_err_mock>0)&(c_out_mock>0)&(c_out_mock>2*c_err_mock)
#p.fill_between(x_out_mock[ok_mock], y1=c_out_mock[ok_mock]-c_err_mock[ok_mock], y2=c_out_mock[ok_mock]+c_err_mock[ok_mock], label='Mock b=0.6 R500c', alpha=0.3, color='r')
#p.plot(x_out_mock[ok_mock], c_out_mock[ok_mock], color='r')

FX_bins = np.arange(8, 18., 0.25)
out = np.cumsum(np.histogram(- DATA_b08['FX_obsF_500'], bins=FX_bins)[0])
c_out_mock = out/area
c_err_mock = c_out_mock * out**(-0.5)
x_out_mock = - 0.5 * (FX_bins[1:] + FX_bins[:-1])
ok_mock = (c_err_mock>0)&(c_out_mock>0)&(c_out_mock>2*c_err_mock)
p.fill_between(x_out_mock[ok_mock], y1=c_out_mock[ok_mock]-c_err_mock[ok_mock], y2=c_out_mock[ok_mock]+c_err_mock[ok_mock], label='Mock r<R500c', alpha=0.1, color='k')
p.plot(x_out_mock[ok_mock], c_out_mock[ok_mock], color='k', label='Mock r<R500c')

#FX_bins = np.arange(8, 18., 0.25)
#out = np.cumsum(np.histogram(- DATA_b10['FX_obsF_500'], bins=FX_bins)[0])
#c_out_mock = out/area
#c_err_mock = c_out_mock * out**(-0.5)
#x_out_mock = - 0.5 * (FX_bins[1:] + FX_bins[:-1])
#ok_mock = (c_err_mock>0)&(c_out_mock>0)&(c_out_mock>2*c_err_mock)
#p.fill_between(x_out_mock[ok_mock], y1=c_out_mock[ok_mock]-c_err_mock[ok_mock], y2=c_out_mock[ok_mock]+c_err_mock[ok_mock], label='Mock b=1.0 R500c', alpha=0.3, color='b')
#p.plot(x_out_mock[ok_mock], c_out_mock[ok_mock], color='b')

p.xlabel('log(F[0.5-2.0 keV])')
p.ylabel('log(>F) [/deg2]')
p.legend(frameon=False, loc=3)
p.yscale('log')
p.xlim((-18, -11.))
p.ylim((2e-4, 10000))
p.title(LC_dir+', z<'+z_max)
p.tight_layout()
p.grid()
p.savefig(fig_out)
p.clf()
print(fig_out, 'written')
t = Table()
t['x'] = x_out_mock[ok_mock]
t['y_mean'] = c_out_mock[ok_mock]
t['y_low'] = c_out_mock[ok_mock]-c_err_mock[ok_mock]
t['y_high'] = c_out_mock[ok_mock]+c_err_mock[ok_mock]
t.write(fig_out[:-3]+'fits', overwrite = True)
print(fig_out[:-3]+'fits', 'written')
