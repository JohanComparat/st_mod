import time
#t0 = time.time()
import os, glob, sys
from astropy.table import Table, vstack
import numpy as np
import astropy.io.fits as fits
from scipy.interpolate import interp1d
nl = lambda sel : len(sel.nonzero()[0])

LC_dir = 'LCerass'

sky_map_hdu = Table.read(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS.fits') )

# merge catalog
# for srv_val in sky_map_hdu['SRVMAP'][(sky_map_hdu['OWNER']==1)]:
for srv_val in sky_map_hdu['SRVMAP'][(sky_map_hdu['OWNER']==2)|(sky_map_hdu['OWNER']==0)]:
    t0 = time.time()
    str_field = str(srv_val).zfill(6)
    t_in = Table.read( os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'AGN_list_sigma_0.8_fsat_8.0.fits') )
    print(len(t_in))
    t_in = t_in[(t_in['FX_soft']>-15.7)&(np.log10(t_in['obs_sm'])>10)&(np.log10(t_in['obs_sm'])<11.9)&(t_in['Mvir']<1e14)]
    print(len(t_in), 'after FX, M*, Mhalo cuts')
    #t_in["SPECTRUM"]

    # indexes_all
    # NH
    data_nh = t_in['logNH']
    data_nh = (data_nh * 5).astype('int') / 5.
    data_nh[data_nh < 20.] = 20.0
    data_nh[data_nh > 26.] = 26.0
    #galactic_nh = np.max([hd_all[1].data['nH'], np.ones_like(hd_all[1].data['nH'])*10**19.9], axis=0)
    #galNH = (10*n.log10(galactic_nh[sel_all])).astype('int')/10.
    data_z = np.round(t_in['redshift_S'],1)

    n_e_bins = 2**np.arange(2, 11)

    n_total = int(512e6)  # /(1.1*N_pixels))
    n_allowed = (n_total / n_e_bins).astype('int')[::-1] - 100

    FX_array = t_in['FX_soft']

    fbins = np.arange(-FX_array.max() - 0.1, - FX_array.min() + 0.1, 0.01)
    xf = fbins[:-1] + 0.01 / 2
    hst = np.cumsum(np.histogram(-np.log10(FX_array), bins=fbins)[0])
    itp = interp1d(hst, xf)

    n_allowed_c = np.cumsum(n_allowed)
    n_allowed_t = n_allowed_c[n_allowed_c < itp.x[-1]]
    FX_inner_boundaries = -itp(n_allowed_t)
    FX_boundaries = np.hstack((
            np.log10(FX_array).max(),
            FX_inner_boundaries,
            np.log10(FX_array).min()))

    data_n_e_b = (np.ones(len(data_z)) * n_e_bins[-1]).astype('int')

    #galNH = t_in['logNH']
    #data_z = t_in['logNH']
    # NH24.2_Z3.9_1024.fits
    # 'NH'+str(np.round(nH,1))+'_Z'+str(np.round(z,1))+'_N'+str(int(nb))+'.fits'
    s1 = 'AGNspectra_V2/'
    #str_galnh = np.array(galNH.astype('str'), dtype=np.object)
    s1b = 'NH'
    str_nh = np.array(data_nh.astype('str'))#, dtype=np.object)
    s2 = '_Z'
    str_z = np.array(data_z.astype('str'))#, dtype=np.object)
    s3 = '_N'
    str_n_e_b = np.array(data_n_e_b.astype('str'))#, dtype=np.object)
    s4 = '.fits' + """[SPECTRUM][#row==1]"""
    #for galNH in np.arange(19.0, 22.6, 0.1):
    #for nH in np.arange(20, 26.2, 0.2):
    #for z in np.arange(0.0, 6.2, 0.1):
    #for nb in 2**np.arange(8, 11):
    #filename = 'galNH' + str(np.round(galNH, 1)) +'_NH' + str(np.round(nH, 1)) + '_Z' + str(np.round(z, 1)) + '_N' + str(int(nb)) + '.fits'
    # str_galnh
    tpl = np.char.add(np.char.add(np.char.add(s1 + s1b , str_nh), np.char.add( s2 , str_z)) , np.char.add( np.char.add(s3 , str_n_e_b ), s4))
    print(tpl)
    # 'galNH' + str(np.round(galNH, 1)) +'_NH' + str(np.round(nH, 1)) + '_Z' + str(np.round(z, 1)) + '_N' + str(int(nb)) + '.fits'
    # /afs/mpe/www/people/comparat/eROSITA_AGN_mock/spectra/Xray/spectra/
    # 'agn/spectra/agn_nH_21.6_z_4.1_nEbins_539.fits[SPECTRUM][#row==1]'

    t_in["SPECTRUM"]=tpl
    p2_simput_out = os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'AGN_list_sigma_0.8_fsat_8.0_simput.fits')
    N_agn_all = len(t_in)
    hdu_cols = fits.ColDefs([
        fits.Column(name="SRC_ID",  format='K',    unit='',    array=(np.arange(len(t_in)).astype('int')),
        fits.Column(name="RA",      format='D',    unit='deg', array=t_in["RA"]),
        fits.Column(name="DEC",     format='D',    unit='deg', array=t_in["DEC"]),
        fits.Column(name="E_MIN",   format='D',    unit='keV', array=np.ones(N_agn_all) * 0.5),
        fits.Column(name="E_MAX",   format='D',    unit='keV', array=np.ones(N_agn_all) * 2.0),
        fits.Column(name="FLUX",    format='D',    unit='erg/s/cm**2', array=10**t_in["FX_soft"]),
        fits.Column(name="SPECTRUM",format='256A', unit='', array=t_in["SPECTRUM"])
    ])
    hdu = fits.BinTableHDU.from_columns(hdu_cols)
    hdu.name = 'SRC_CAT'
    hdu.header['HDUCLASS'] = 'HEASARC/SIMPUT'
    hdu.header['HDUCLAS1'] = 'SRC_CAT'
    hdu.header['HDUVERS'] = '1.1.0'
    hdu.header['RADESYS'] = 'FK5'
    hdu.header['EQUINOX'] = 2000.0
    outf = fits.HDUList([fits.PrimaryHDU(), hdu])  # ,  ])
    #if os.path.isfile(p2_simput_out):
        #os.system("rm " + p2_simput_out)
    outf.writeto(p2_simput_out, overwrite=True)
    print(p2_simput_out, 'written', time.time() - t0)

