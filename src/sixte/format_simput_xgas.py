import time
#t0 = time.time()
import os, glob, sys
from astropy.table import Table, vstack
import numpy as np
import astropy.io.fits as fits

nl = lambda sel : len(sel.nonzero()[0])

LC_dir = 'LCerass'

sky_map_hdu = Table.read(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS.fits') )

# merge catalog
for srv_val in sky_map_hdu['SRVMAP']:
    t0 = time.time()
    str_field = str(srv_val).zfill(6)
    t_in = Table.read( os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'Xgas_bHS0.8_simput.fits') )
    p2_simput_out = os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'Xgas_bHS0.8_simput.fits')
    N_clu_all = len(t_in)
    hdu_cols = fits.ColDefs([
        fits.Column(name="SRC_ID",  format='K',    unit='',    array=(t_in['SRC_ID']).astype('int')),
        fits.Column(name="RA",      format='D',    unit='deg', array=t_in["RA"]),
        fits.Column(name="DEC",     format='D',    unit='deg', array=t_in["DEC"]),
        fits.Column(name="E_MIN",   format='D',    unit='keV', array=np.ones(N_clu_all) * 0.5),
        fits.Column(name="E_MAX",   format='D',    unit='keV', array=np.ones(N_clu_all) * 2.0),
        fits.Column(name="FLUX",    format='D',    unit='erg/s/cm**2', array=t_in["FLUX"]),
        fits.Column(name="IMAGE",   format='100A', unit='', array=t_in["IMAGE"]),
        fits.Column(name="SPECTRUM",format='100A', unit='', array=t_in["SPECTRUM"]),
        fits.Column(name="IMGROTA", format='D',    unit='deg', array=t_in["IMGROTA"]),
        fits.Column(name="IMGSCAL", format='D',    unit='', array=t_in["IMGSCAL"])
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

