import time
#t0 = time.time()
import os, glob, sys
from astropy.table import Table, vstack
import numpy as np
from scipy.interpolate import interp1d
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as plt

nl = lambda sel : len(sel.nonzero()[0])

LC_dir = 'LCerass'

sky_map_hdu = Table.read(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS.fits') )
validation_dir           = os.path.join(os.environ['GIT_STMOD_DATA'], 'data', 'validation','validation_AGN')
validation_dir_lNlS_soft = os.path.join(validation_dir, 'XraySoftLogNlogS', 'LCerass')
os.system('mkdir -p '+validation_dir_lNlS_soft)
fdex = 0.05
fbins = np.arange(-18, -8, fdex)
# area_per_tile = 3*3

d_dec = (np.sin(sky_map_hdu['DE_MAX']*np.pi/180.)-np.sin(sky_map_hdu['DE_MIN']*np.pi/180.))
sky_map_hdu['AREA'] = (sky_map_hdu['RA_MAX']-sky_map_hdu['RA_MIN']) * d_dec * 180/np.pi

# merge catalog
#for srv_val in sky_map_hdu['SRVMAP'][(sky_map_hdu['OWNER']==1)]:
hists=[]
for srv_val in sky_map_hdu['SRVMAP'][(sky_map_hdu['OWNER']==2)|(sky_map_hdu['OWNER']==0)]:
    t0 = time.time()
    str_field = str(srv_val).zfill(6)
    p_2_AGN = os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'AGN_list_sigma_0.8_fsat_8.0.fits')
    p_2_AGN_lNlS = os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'logNlogS_AGN_list_sigma_0.8_fsat_8.0.fits')
    if os.path.isfile(p_2_AGN):
        AGN = Table.read(p_2_AGN)
        z = AGN['redshift_S']
        #dl2_cm = np.log10(dl_itp(z)) * 2 + np.log10(4*np.pi)
        # fx_hard = AGN['FX_hard'] #- dl2_cm
        fx_soft = AGN['FX_soft'] #- dl2_cm
        t_out = Table()
        t_out['FX_lo'] = fbins[:-1]
        t_out['FX_hi'] = fbins[1:]
        # t_out['N_hard'] = np.histogram(fx_hard, fbins)[0]
        t_out['N_soft'] = np.histogram(fx_soft, fbins)[0]
        #t_out['area'] = area_tile
        # print(t_out)
        t_out.write(p_2_AGN_lNlS, overwrite=True)
        print(p_2_AGN_lNlS, 'written')
        hists.append(t_out['N_soft'])
hists = np.array(hists)


Simput_hists=[]
Simput_areas=[]
for skmid in sky_map_hdu[(sky_map_hdu['OWNER']==2)|(sky_map_hdu['OWNER']==0)]:
    srv_val = skmid['SRVMAP']
    t0 = time.time()
    str_field = str(srv_val).zfill(6)
    p2_simput_out = os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'AGN_list_sigma_0.8_fsat_8.0_simput.fits')
    if os.path.isfile(p2_simput_out):
        AGN = Table.read(p2_simput_out)
        fx_soft = np.log10(AGN['FLUX'])
        Simput_hists.append(np.histogram(fx_soft, fbins)[0])
        Simput_areas.append(skmid['AREA'])

Simput_hists=np.array(Simput_hists)
Simput_areas=np.array(Simput_areas)


hists2=[]
areas = []
for skmid in sky_map_hdu[(sky_map_hdu['OWNER']==2)|(sky_map_hdu['OWNER']==0)]:
    srv_val = skmid['SRVMAP']
    t0 = time.time()
    str_field = str(srv_val).zfill(6)
    p_2_AGN_lNlS = os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'logNlogS_AGN_list_sigma_0.8_fsat_8.0.fits')
    if os.path.isfile(p_2_AGN_lNlS):
        t_out = Table.read(p_2_AGN_lNlS)
        hists2.append(t_out['N_soft'])
        areas.append(skmid['AREA'])
hists2 = np.array(hists2)
areas = np.array(areas)

hists = hists2 # np.array(hists)

x_data1 = (t_out['FX_lo']+t_out['FX_hi'])/2.
y_data1 = np.sum(hists, axis=0)/(np.sum(areas))
Simput_y_data1 = np.sum(Simput_hists, axis=0)/(np.sum(Simput_areas))
ref_line = interp1d( x_data1, y_data1)

plt.figure(1, (6, 6))
for nn, aa in zip(hists, areas):
    plt.plot(x_data1, nn/ aa, lw=0.8, ls='solid', color='k', alpha=0.3, rasterized=True)

# Georgakakis 2008
path_2_logNlogS_data = os.path.join(
    os.environ["GIT_STMOD_DATA"],
    'data/validation/validation_AGN/literature_data',
    'logNlogS_Georgakakis_08_AGN.data')
x_data, y_data = np.loadtxt(path_2_logNlogS_data, unpack=True)
plt.plot(np.log10(x_data),y_data, lw=3, ls='dotted', color='g', label='Ge08')

# Merloni 2012
path_2_logNlogS_data = os.path.join(
    os.environ["GIT_STMOD_DATA"],
    'data/validation/validation_AGN/literature_data',
    'logNlogS_Merloni_12_AGN.data')
x_data, y_data = np.loadtxt(path_2_logNlogS_data, unpack=True)
plt.plot(np.log10(x_data),y_data, lw=3, ls='dotted', color='r', label='Me12')

# Mateos 2008
path_2_logNlogS_data = os.path.join(
    os.environ["GIT_STMOD_DATA"],
    'data/validation/validation_AGN/literature_data',
    'logNlogS_Mateos_08_AGN.data')
x_data, y_data, err = np.loadtxt(path_2_logNlogS_data, unpack=True)
plt.plot(x_data,y_data, lw=3, ls='dotted', color='b', label='Ma08')

plt.xlabel(r'$\log_{10}(F_X$[0.5-2 keV])')
plt.ylabel(r'$N(>F_X)$ [/deg2]')
plt.legend(frameon=False, loc=3)
plt.yscale('log')
plt.xlim((-18, -11))
plt.ylim((1e-3, 1e5))
#plt.grid()
plt.tight_layout()
plt.savefig(os.path.join(validation_dir_lNlS_soft, LC_dir+"_logN_logS_soft_AGN_pertile.png"))
plt.clf()
print(os.path.join(validation_dir_lNlS_soft, LC_dir+"_logN_logS_soft_AGN_pertile.png"), 'written')

plt.figure(1, (6, 6))
plt.plot(x_data1, y_data1, lw=2, ls='solid', label='AGN, this work')
plt.plot(x_data1, y_data1*8, lw=1, ls='solid', label='AGN, this work *8')
plt.plot(x_data1, Simput_y_data1, lw=2, ls='solid', label='AGN, Simput')

# Georgakakis 2008
path_2_logNlogS_data = os.path.join(
    os.environ["GIT_STMOD_DATA"],
    'data/validation/validation_AGN/literature_data',
    'logNlogS_Georgakakis_08_AGN.data')
x_data, y_data = np.loadtxt(path_2_logNlogS_data, unpack=True)
plt.plot(np.log10(x_data),y_data, lw=3, ls='dotted', color='g', label='Ge08')

# Merloni 2012
path_2_logNlogS_data = os.path.join(
    os.environ["GIT_STMOD_DATA"],
    'data/validation/validation_AGN/literature_data',
    'logNlogS_Merloni_12_AGN.data')
x_data, y_data = np.loadtxt(path_2_logNlogS_data, unpack=True)
plt.plot(np.log10(x_data),y_data, lw=3, ls='dotted', color='r', label='Me12')

# Mateos 2008
path_2_logNlogS_data = os.path.join(
    os.environ["GIT_STMOD_DATA"],
    'data/validation/validation_AGN/literature_data',
    'logNlogS_Mateos_08_AGN.data')
x_data, y_data, err = np.loadtxt(path_2_logNlogS_data, unpack=True)
plt.plot(x_data,y_data, lw=3, ls='dotted', color='b', label='Ma08')

plt.xlabel(r'$\log_{10}(F_X$[0.5-2 keV])')
plt.ylabel(r'$N(>F_X)$ [/deg2]')
plt.legend(frameon=False, loc=3)
plt.yscale('log')
plt.xlim((-18, -11))
plt.ylim((1e-3, 1e5))
#plt.grid()
plt.tight_layout()
plt.savefig(os.path.join(validation_dir_lNlS_soft, LC_dir+"_logN_logS_soft_AGN.png"))
plt.clf()
print(os.path.join(validation_dir_lNlS_soft, LC_dir+"_logN_logS_soft_AGN.png"), 'written')


plt.figure(1, (6, 6))


# Georgakakis 2008
path_2_logNlogS_data = os.path.join(
    os.environ["GIT_STMOD_DATA"],
    'data/validation/validation_AGN/literature_data',
    'logNlogS_Georgakakis_08_AGN.data')
x_data, y_data = np.loadtxt(path_2_logNlogS_data, unpack=True)
plt.plot(np.log10(x_data),y_data/ref_line(np.log10(x_data)), lw=3, ls='dotted', color='g', label='G08')

# Merloni 2012
path_2_logNlogS_data = os.path.join(
    os.environ["GIT_STMOD_DATA"],
    'data/validation/validation_AGN/literature_data',
    'logNlogS_Merloni_12_AGN.data')
x_data, y_data = np.loadtxt(path_2_logNlogS_data, unpack=True)
plt.plot(np.log10(x_data),y_data/ref_line(np.log10(x_data)), lw=3, ls='dotted', color='r', label='M12')

# Mateos 2008
path_2_logNlogS_data = os.path.join(
    os.environ["GIT_STMOD_DATA"],
    'data/validation/validation_AGN/literature_data/logNlogS_Mateos_08_AGN.data')
x_data, y_data, err = np.loadtxt(path_2_logNlogS_data, unpack=True)
plt.plot(x_data,y_data/ref_line(x_data), lw=3, ls='dotted', color='b', label='M08')

plt.xlabel('log10(F_X[0.5-2 keV])')
plt.ylabel('N(>F_X)/(N AGN,>F_X) ')
plt.legend(frameon=False, loc=0)
plt.xlim((-19, -11.5))
plt.ylim((0.3, 1.7))
plt.tight_layout()
#plt.grid()
plt.savefig(os.path.join(validation_dir_lNlS_soft, LC_dir+"_logN_logS_soft_AGN_ratio.png"))
plt.clf()


z_bins = np.arange(0, 4.1, 0.01)
# merge catalog
#for srv_val in sky_map_hdu['SRVMAP'][(sky_map_hdu['OWNER']==1)]:
zhists=[]
for srv_val in sky_map_hdu['SRVMAP'][(sky_map_hdu['OWNER']==2)|(sky_map_hdu['OWNER']==0)]:
    t0 = time.time()
    str_field = str(srv_val).zfill(6)
    p_2_AGN = os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'AGN_list_sigma_0.8_fsat_8.0.fits')
    p_2_AGN_lNlS = os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'zhist_AGN_list_sigma_0.8_fsat_8.0.fits')
    if os.path.isfile(p_2_AGN):
        AGN = Table.read(p_2_AGN)
        t_out = Table()
        t_out['N_z'] = np.histogram(AGN['redshift_S'], z_bins)[0]
        #t_out['area'] = area_tile
        # print(t_out)
        t_out.write(p_2_AGN_lNlS, overwrite=True)
        print(p_2_AGN_lNlS, 'written')
        zhists.append(t_out['N_z'])
zhists = np.array(zhists)


areas = []
for skmid in sky_map_hdu[(sky_map_hdu['OWNER']==2)|(sky_map_hdu['OWNER']==0)]:
    srv_val = skmid['SRVMAP']
    t0 = time.time()
    str_field = str(srv_val).zfill(6)
    p_2_AGN = os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'AGN_list_sigma_0.8_fsat_8.0.fits')
    if os.path.isfile(p_2_AGN):
        areas.append(skmid['AREA'])
areas = np.array(areas)



x_data1 = (z_bins[:-1]+z_bins[1:])/2.
y_data1 = np.sum(zhists, axis=0)/(np.sum(areas))

plt.figure(1, (6, 6))
for nn, aa in zip(zhists, areas):
    plt.plot(x_data1, nn/ aa, lw=0.1, ls='solid', color='k', alpha=0.3, rasterized=True)


plt.plot(x_data1, y_data1, lw=2, ls='solid', label='AGN, this work')#suffix)

plt.xlabel(r'$z$')
plt.ylabel(r'$N$ [/deg2 dz=0.01]')
plt.legend(frameon=False, loc=3)
plt.yscale('log')
plt.xlim((0,4))
#plt.ylim((1e-3, 1e5))
#plt.grid()
plt.tight_layout()
plt.savefig(os.path.join(validation_dir_lNlS_soft, LC_dir+"_zhist_AGN_pertile.png"))
plt.clf()
print(os.path.join(validation_dir_lNlS_soft, LC_dir+"_zhist_AGN_pertile.png"), 'written')

