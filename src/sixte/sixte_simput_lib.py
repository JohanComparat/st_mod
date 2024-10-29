import numpy as np
from scipy.interpolate import interp1d
import astropy.io.fits as fits
import time
import os

def make_agn_simput(t_agn, path_2_SMPT_catalog='test.fits', FX_LIM_value_cen = -16):
    t0 = time.time()
    t_agn = t_agn[ (t_agn['FX_soft'] > FX_LIM_value_cen) ]

    data_nh = t_agn['logNH']
    data_nh = (data_nh * 5).astype('int') / 5.
    data_nh[data_nh < 20.] = 20.0
    data_nh[data_nh > 26.] = 26.0
    galactic_nh = np.max([t_agn['nH'], np.ones_like(t_agn['nH'])*10**19.9], axis=0)
    galNH = (10*np.log10(galactic_nh)).astype('int')/10.
    data_z = np.round(t_agn['redshift_S'],1)
    data_indexes = t_agn['line_ID']
    ra_array = t_agn['RA']
    dec_array = t_agn['DEC']

    n_e_bins = 2**np.arange(2, 11)
    n_total = int(512e6)
    n_allowed = (n_total / n_e_bins).astype('int')[::-1] - 100
    FX_array = 10**t_agn['FX_soft']
    fbins = np.arange(-np.log10(FX_array).max() - 0.1, -  np.log10(FX_array).min() + 0.1, 0.01)
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

    for jj, (f_min, f_max) in enumerate(zip(FX_boundaries[:-1], FX_boundaries[1:])):
        selection = (np.log10(FX_array) <= f_min) & (np.log10(FX_array) >= f_max)
        n_e_val = n_e_bins[::-1][jj]
        data_n_e_b[selection] = n_e_val
        print(f_min,
                f_max,
                n_e_val,
                len(data_n_e_b[selection]) < 512e6 / n_e_val,
                len(data_n_e_b[selection]),
                512e6 / n_e_val,
                data_n_e_b[selection])

    s1 = 'agn_Xspectra/galNH'
    str_galnh = np.array(galNH.astype('str'), dtype=object)
    s1b = '_NH'
    str_nh = np.array(data_nh.astype('str'), dtype=object)
    s2 = '_Z'
    str_z = np.array(data_z.astype('str'), dtype=object)
    s3 = '_N'
    str_n_e_b = np.array(data_n_e_b.astype('str'), dtype=object)
    s4 = '.fits' + """[SPECTRUM][#row==1]"""

    tpl = s1 + str_galnh + s1b + str_nh + s2 + str_z + s3 + str_n_e_b + s4
    print(tpl)
    # 'galNH' + str(np.round(galNH, 1)) +'_NH' + str(np.round(nH, 1)) + '_Z' + str(np.round(z, 1)) + '_N' + str(int(nb)) + '.fits'
    # /afs/mpe/www/people/comparat/eROSITA_AGN_mock/spectra/Xray/spectra/
    # 'agn/spectra/agn_nH_21.6_z_4.1_nEbins_539.fits[SPECTRUM][#row==1]'

    N_agn_out = len(ra_array)

    hdu_cols = fits.ColDefs([
        fits.Column(name="SRC_ID", format='K', unit='', array=data_indexes),
        fits.Column(name="RA", format='D', unit='deg', array=ra_array),
        fits.Column(name="DEC", format='D', unit='deg', array=dec_array),
        fits.Column(name="E_MIN", format='D', unit='keV', array=np.ones(N_agn_out) * 0.5),
        fits.Column(name="E_MAX", format='D', unit='keV', array=np.ones(N_agn_out) * 2.0),
        fits.Column(name="FLUX", format='D', unit='erg/s/cm**2', array=FX_array),
        fits.Column(name="SPECTRUM", format='100A', unit='', array=tpl),
        fits.Column(name="n_energy_bins", format='K', unit='', array=data_n_e_b)
    ])

    hdu = fits.BinTableHDU.from_columns(hdu_cols)

    hdu.name = 'SRC_CAT'
    hdu.header['HDUCLASS'] = 'HEASARC/SIMPUT'
    hdu.header['HDUCLAS1'] = 'SRC_CAT'
    hdu.header['HDUVERS'] = '1.1.0'
    hdu.header['RADESYS'] = 'FK5'
    hdu.header['EQUINOX'] = 2000.0

    outf = fits.HDUList([fits.PrimaryHDU(), hdu])  # ,  ])
    if os.path.isfile(path_2_SMPT_catalog):
        os.system("rm " + path_2_SMPT_catalog)
    outf.writeto(path_2_SMPT_catalog, overwrite=True)
    print(path_2_SMPT_catalog, 'written', time.time() - t0)


