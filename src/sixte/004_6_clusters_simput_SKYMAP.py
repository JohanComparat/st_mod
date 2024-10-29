"""
What it does
------------

Creates a simput catalog fo each healpix pixel (healpix pixel of 13.4 deg2)

References
----------

Command to run
--------------

python3 004_6_clusters_simput.py environmentVAR

arguments
---------

environmentVAR: environment variable linking to the directory where files are e.g. "MD10"
It will then work in the directory : $environmentVAR/hlists/fits/

Dependencies
------------

import time, os, sys, glob, numpy, astropy, scipy, matplotlib

/data40s/erosim/eRASS/eRASS8_cluster/000/erass_ccd1_evt.fits
/data17s/darksim/MD/MD_1.0Gpc/cat_CLU_SIMPUT/
"""
from scipy.special import erf
from astropy.table import Table, Column
import sys, os, time
import astropy.units as u
from astropy_healpix import healpy
import astropy.io.fits as fits
import numpy as n
from scipy.interpolate import interp1d
print('CREATES SIMPUT CLUSTER FILES')
print('------------------------------------------------')
print('------------------------------------------------')
t0 = time.time()

#env = "UNIT_fA1i_DIR" #sys.argv[1] #
#pixel_size_image = 20. #float(sys.argv[2]) # 20 or 2
#b_HS = 8.0 # float(sys.argv[3])
#cov_mat_option = '0' #sys.argv[4]
#logM500c_min = 13.0 # float(sys.argv[5]) # 13.0
#logFX_min = -14.5 # float(sys.argv[6]) # 13.0

#root_dir = os.path.join(os.environ[env])

#path_2_catalog = os.path.join(root_dir, "UNIT_fA1i_DIR_eRO_CLU_b8_CM_0_pixS_20.0_M500c_13.0_FX_-14.5_MGAS_Sept2021.fits")
#dir_2_SMPT = os.path.join(os.environ[env], "SIMPUT_SKYMAP_UNIT_fA1i_DIR_eRO_CLU_b8_CM_0_pixS_20.0_M500c_13.0_FX_-14.5_MGAS_Sept2021" )


env = "MD40" #sys.argv[1] #
pixel_size_image = 10. #float(sys.argv[2]) # 20 or 2
b_HS = 8.0 # float(sys.argv[3])
cov_mat_option = '0' #sys.argv[4]
logM500c_min = 13.0 # float(sys.argv[5]) # 13.0
logFX_min = -14.5 # float(sys.argv[6]) # 13.0

root_dir = os.path.join(os.environ[env])

path_2_catalog = os.path.join(root_dir, "MD40_eRO_CLU_b8_CM_0_pixS_10.0_M500c_13.0_FX_-14.5_25Apr2022_Profiles.fits")
dir_2_SMPT = os.path.join(os.environ[env], "cat_CLU_SIMPUT_b8_CM_0_pixS_10.0_M500c_13.0_FX_-14.5_image_yes" )


dir_2_SMPT_image = os.path.join(dir_2_SMPT, 'cluster_images')
dir_2_SMPT_spectra = os.path.join(dir_2_SMPT, 'cluster_Xspectra')


sky_map_hdu = Table.read(os.path.join(os.environ['GIT_ERASS_SIM'], 'data', 'SKYMAPS.fits'))
import pymangle
eRASS_ply = os.path.join(os.environ['GIT_ERASS_SIM'], 'data', 'eRASS_SKYMAPS.ply' )
mng = pymangle.Mangle( eRASS_ply )

print(path_2_catalog)
print(dir_2_SMPT)
print(dir_2_SMPT_image)
print(dir_2_SMPT_spectra)

if os.path.isdir(dir_2_SMPT) == False:
    os.system('mkdir -p ' + dir_2_SMPT)

# simulation setup
if env[:2] == "MD" : # env == "MD04" or env == "MD40" or env == "MD10" or env == "MD25"
    from astropy.cosmology import FlatLambdaCDM
    import astropy.units as u
    cosmoMD = FlatLambdaCDM(
        H0=67.77 * u.km / u.s / u.Mpc,
        Om0=0.307115)  # , Ob0=0.048206)
    h = 0.6777
    L_box = 1000.0 / h
    cosmo = cosmoMD
if env[:4] == "UNIT" : # == "UNIT_fA1_DIR" or env == "UNIT_fA1i_DIR" or env == "UNIT_fA2_DIR":
    from astropy.cosmology import FlatLambdaCDM
    import astropy.units as u
    cosmoUNIT = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
    h = 0.6774
    L_box = 1000.0 / h
    cosmo = cosmoUNIT

# attenuation of the flux for the groups
#def f_LX_att(log10KT): return (0.51 + 0.5 * erf((log10KT - 0) / 0.2))

t_1 = Table.read(path_2_catalog)
s_1 = ( t_1['HALO_pid']==-1 ) & ( t_1['galaxy_SMHMR_mass']>10. ) & ( n.log10(t_1['HALO_Mvir']) > 13 ) & ( t_1['detectable'] ) & ( t_1['g_lon']>170 )
t = t_1[s_1]
print(len(t_1), len(t))

t['scale_2_2R500c'] = 10**(t['CLUSTER_LX_soft_RF_twiceR500']-t['CLUSTER_LX_soft_RF'])
t['FX_soft_SIMPUT'] = t['CLUSTER_FX_soft'] * t['scale_2_2R500c'] 
N_detec = len(t['CLUSTER_FX_soft'])
print('N total=', N_detec)
print('density=', N_detec*n.pi/129600., '/deg2')

# polygon ID from the mangle mask
polyid = mng.polyid( t['RA'], t['DEC'] )

# link to templates
#def tpl_name(temperature, redshift): return 'cluster_Xspectra/cluster_spectrum_10kT_' + str(int(temperature * 1000)).zfill(4) + '_100z_' + str(int(redshift * 10000)).zfill(4) + '.fits[SPECTRUM][#row==1]'
def tpl_name(temperature, redshift): return 'cluster_Xspectra/cluster_spectrum_10000kT_' + str(int(temperature * 10000)).zfill(7) + '_10000z_' + str(int(redshift * 10000)).zfill(7) + '.fits[SPECTRUM][#row==1]'
#HEALPIX_8_id = 151
def tpl_name(temperature, redshift, nh_val): return  'cluster_Xspectra/galNH_' + str(n.round(nh_val, 3)) +'_10000kT_' + str(int(10000*temperature)) + '_10000z_' + str(int(10000*redshift)) + '.fits'+ """[SPECTRUM][#row==1]"""

Area_per_field =   8.777225797748782

for sky_tile in sky_map_hdu[(sky_map_hdu['OWNER']==2)|(sky_map_hdu['OWNER']==0)] :
    """
    Loops over healpix pixels and writes the files to path_2_eRO_catalog
    """
    sky_tile_id = sky_tile['SRVMAP']
    print(sky_tile_id)
    RA_CEN = sky_tile['RA_CEN']
    DEC_CEN = sky_tile['DE_CEN']
    print('polyids obtained', n.unique(polyid), n.unique(polyid) == sky_tile_id)
    sf = (polyid == sky_tile_id)
    t2 = t[sf]
    N_clu_all = len(t2['RA'])
    print('Number of clusters=',N_clu_all)
    print('density=',N_clu_all/Area_per_field, '/deg2')
    if N_clu_all>0:
        ra_array = t2['RA']
        dec_array = t2['DEC']
        redshift = t2['redshift_R']
        kT = t2['CLUSTER_kT']
        galactic_nh = n.max([t2['nH'], n.ones_like(t2['nH'])*10**19.9], axis=0)
        galNH = (10*n.log10(galactic_nh)).astype('int')/10.
        # size of the pixel in the image written
        # randomize orientations
        rd_all = n.random.rand(N_clu_all)
        orientation = n.random.rand(N_clu_all) * 180.  # IMGROTA
        # scale the image with the size of the cluster
        # all images have 5.5e-04*120*60 = 3.96 arc minute on the side
        # default size 0.033 arcmin/pixel
        pixel_rescaling =  n.ones_like(t2['angularSize_per_pixel'])
        # NOW ASSIGNS TEMPLATES BASED ON THE HALO PROPERTIES
        #template = n.zeros(N_clu_all).astype('U100')
        #template[template == "0.0"] = "cluster_images/elliptical_ba_0p25_cc.fits[SPECTRUM][#row==1]"
        template = n.array([ el.strip()+".fits[IMAGE]" for el in t2['XRAY_image_path'] ])
        template_fullPath = n.array([ os.path.join( dir_2_SMPT, el.strip()+".fits")  for el in t2['XRAY_image_path'] ])
        template_exists = n.array([ os.path.isfile( os.path.join( dir_2_SMPT, el.strip()+".fits") ) for el in t2['XRAY_image_path'] ])
        N_exist = len(template_exists.nonzero()[0])
        N_templates = len(template_exists)
        if N_exist<N_templates:
            print('error, missing images', N_exist, N_templates)
        # NOW links to the grid of SPECTRA
        # kt_arr = n.array([0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])
        # z_arr = n.hstack((n.array([0., 0.05]), n.arange(0.1, 4., 0.1)))
        #kt_arr = 10**n.arange(-1,1.3,0.01)
        #z_arr = n.hstack((n.array([0., 0.01]), 10**n.arange(n.log10(0.02), n.log10(4.), 0.01)))
        #
        kt_arr = 10**n.arange(-1,1.3,0.05)
        z_arr = n.hstack((n.array([0., 0.01]), 10**n.arange(n.log10(0.02), n.log10(4.), 0.05)))
        #galNH = n.arange(19.0, 22.6, 0.1)
        indexes_kt = n.array([(n.abs(kT_val - kt_arr)).argmin() for kT_val in kT])
        kT_values = kt_arr[indexes_kt]
        indexes_z = n.array([(n.abs(z_val - z_arr)).argmin() for z_val in redshift])
        z_values = z_arr[indexes_z]
        spec_names = n.zeros(N_clu_all).astype('U200')
        # "cluster_Xspectra/cluster_spectrum_10kT_0100_100z_0150.fits[SPECTRUM][#row==1]"
        for jj, (kT_values_ii, z_values_ii, galNH_ii) in enumerate(zip(kT_values, z_values, galNH)):
            spec_names[jj] = tpl_name(kT_values_ii, z_values_ii, galNH_ii)

        N_per_simput = 999
        for jj, (id_min, id_max) in enumerate(zip(n.arange(0,N_clu_all,N_per_simput), n.arange(0,N_clu_all,N_per_simput)+N_per_simput)):
            path_2_SMPT_catalog = os.path.join(dir_2_SMPT, 'c_'+str(sky_tile_id).zfill(6) + '_N_'+str(jj)+'.fit')
            hdu_cols = fits.ColDefs([
                fits.Column(name="SRC_ID",  format='K',    unit='',    array=(n.arange(N_clu_all) + 4e8).astype('int')[id_min:id_max]),
                fits.Column(name="RA",      format='D',    unit='deg', array=ra_array[id_min:id_max]),
                fits.Column(name="DEC",     format='D',    unit='deg', array=dec_array[id_min:id_max]),
                fits.Column(name="E_MIN",   format='D',    unit='keV', array=n.ones(N_clu_all)[id_min:id_max] * 0.5),
                fits.Column(name="E_MAX",   format='D',    unit='keV', array=n.ones(N_clu_all)[id_min:id_max] * 2.0),
                fits.Column(name="FLUX",    format='D',    unit='erg/s/cm**2', array=t2['FX_soft_SIMPUT'][id_min:id_max]),
                fits.Column(name="IMAGE",   format='100A', unit='', array=template[id_min:id_max]),
                fits.Column(name="SPECTRUM",format='100A', unit='', array=spec_names[id_min:id_max]),
                fits.Column(name="IMGROTA", format='D',    unit='deg', array=orientation[id_min:id_max]),
                fits.Column(name="IMGSCAL", format='D',    unit='', array=pixel_rescaling[id_min:id_max])
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
