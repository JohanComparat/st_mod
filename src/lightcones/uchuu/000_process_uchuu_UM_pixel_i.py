"""

Replicates the snapshot over the sky and extracts halos in a given shell of comoving distance
Does full sky replication

Parameters
==========

 * l_box: length of the box in Mpc/h, example 400
 * h: unitless hubble parameter h=H0/100mk/s/Mpc, example 0.6777
 * N_skip: how many lines of header should be skipped while reading the halo file
 * env: environment variable linking to the simulation directory, example MD04 MD10, UNIT_fA1_DIR ...
 * path_2_header: path to the header of the output file
 * path_2_in_0: path to the rockstar hlist file, example /path/2/simulation/hlist_1.00000.h5

/home/users/dae/comparat/Uchuu/hlists/halodir_050/halolist_z0p00_0.h5
/home/users/dae/comparat/Uchuu/LC/halodir_050/halolist_z0p00_0.h5

Processing
==========

Generates and executes commands to replicate and filter the halo files. Lines 120 - 129 are the most important lines of the script.

outputs
=======

 * directory(path_2_in_0)/replication_$expansionParameter/all the files: contains a single file per replication. Deleted after the computation finishes.
 * directory(path_2_in_0)/replicated_$expansionParameter/one summary file: contains one file with all replication concatenated. Could be deleted after the computation is done.
 * $env/fits/all_'+$expansionParameter+'.fits') the former file converted to fits staged into the final directory. Keep this one ! It is the shell containing all the halos


import numpy as n
import os, sys, glob
import healpy


all_hlist_snap = n.array([ glob.glob( os.path.join( os.environ['UCHUU'], 'hlists', 'halodir_'+str(snap_id).zfill(3), 'halolist_*.h5' ) ) for snap_id in n.arange(10,51,1) ])
n.ravel(all_hlist_snap)
print(all_hlist_snap.shape)
NSIDE=16
all_npix = n.arange( healpy.nside2npix( NSIDE ) )
for pix in all_npix[:10]:
	f=open("pixel_batches_uchuu/pixel_batch_"+str(pix).zfill(5)+".sh", 'w')
	f.write("#!/bin/bash")
	f.write("\n")
	f.write("cd /home/users/dae/comparat/software/lss_mock_dev/python/DM_LC/")
        f.write("\n")
	for file_name in n.ravel(all_hlist_snap):	
		f.write("python 000_process_uchuu_pixel_i.py "+file_name+" "+str(pix)+" \n")
	f.close()

cd pixel_batches_uchuu
ipython

import numpy as n
import os, sys, glob
import healpy

batch_list = n.array(glob.glob("pixel_batch_?????.sh"))
batch_list.sort()

for ii in n.arange(31):
	os.system('rm '+"batch_"+str(ii).zfill(ii)+".sh")
        f=open("batch_"+str(ii).zfill(2)+".sh", 'w')
        f.write("#!/bin/bash")
        f.write("\n")
	print(len(batch_list[100*ii:100*(ii+1)]))
	for todo in batch_list[100*ii:100*(ii+1)]:
		f.write("sh "+todo+" \n")
	f.close()

nohup sh batch_00.sh > log_batch_00.log &
nohup sh batch_01.sh > log_batch_01.log &
nohup sh batch_02.sh > log_batch_02.log &
nohup sh batch_03.sh > log_batch_03.log &
nohup sh batch_04.sh > log_batch_04.log &
nohup sh batch_05.sh > log_batch_05.log &
nohup sh batch_06.sh > log_batch_06.log &
nohup sh batch_07.sh > log_batch_07.log &
nohup sh batch_08.sh > log_batch_08.log &
nohup sh batch_09.sh > log_batch_09.log &
nohup sh batch_10.sh > log_batch_10.log &
nohup sh batch_11.sh > log_batch_11.log &
nohup sh batch_12.sh > log_batch_12.log &
nohup sh batch_13.sh > log_batch_13.log &
nohup sh batch_14.sh > log_batch_14.log &
nohup sh batch_15.sh > log_batch_15.log &
nohup sh batch_16.sh > log_batch_16.log &
nohup sh batch_17.sh > log_batch_17.log &
nohup sh batch_18.sh > log_batch_18.log &
nohup sh batch_19.sh > log_batch_19.log &
nohup sh batch_20.sh > log_batch_20.log &
nohup sh batch_21.sh > log_batch_21.log &
nohup sh batch_22.sh > log_batch_22.log &
nohup sh batch_23.sh > log_batch_23.log &
nohup sh batch_24.sh > log_batch_24.log &
nohup sh batch_25.sh > log_batch_25.log &
nohup sh batch_26.sh > log_batch_26.log &
nohup sh batch_27.sh > log_batch_27.log &
nohup sh batch_28.sh > log_batch_28.log &
nohup sh batch_29.sh > log_batch_29.log &
nohup sh batch_30.sh > log_batch_30.log &


z_names = n.array([  "z0p00"
					,"z0p02"
					,"z0p05"
					,"z0p09"
					,"z0p14"
					,"z0p19"
					,"z0p25"
					,"z0p30"
					,"z0p36"
					,"z0p43"
					,"z0p49"
					,"z0p56"
					,"z0p63"
					,"z0p70"
					,"z0p78"
					,"z0p86"
					,"z0p94"
					,"z1p03"
					,"z1p12"
					,"z1p22"
					,"z1p32"
					,"z1p43"
					,"z1p54"
					,"z1p65"
					,"z1p77"
					,"z1p90"
					,"z2p03"
					,"z2p17"
					,"z2p31"
					,"z2p46"
					,"z2p62"
					,"z2p78"
					,"z2p95"
					,"z3p13"
					,"z3p32"
					,"z3p61"
					,"z3p93"
					,"z4p27"
					,"z4p63"
					,"z5p15"
					,"z5p73"
					,"z6p35"
					,"z7p03"
					,"z7p76"
					,"z8p59"
					,"z9p47"
					,"z10p44"
					,"z11p51"
					,"z12p69"
					,"z13p96" ])

"""
print('runs 001_process_hlists.py with arguments' )
import time
t0=time.time()
import numpy as n
import os, sys, subprocess
print(sys.argv)
import h5py
from astropy.table import Table, Column
from astropy.coordinates import SkyCoord
import healpy
from scipy.interpolate import interp1d
import astropy.io.fits as fits
from scipy.special import erf
from scipy.stats import norm
import scipy
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

from dustmaps.planck import PlanckQuery
planck = PlanckQuery()

path_2_NH_map = '/ptmp/joco/h1_maps/H14PI/asu.fit'
NH_DATA2 = fits.open(path_2_NH_map)[1].data
nh_val = NH_DATA2['NHI']

cosmoUNIT = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
h = 0.6774
cosmo = cosmoUNIT

z_array = n.arange(0, 7.5, 0.001)
dc_to_z = interp1d(cosmo.comoving_distance(z_array), z_array)

d_L = cosmo.luminosity_distance(z_array)
dl_cm = (d_L.to(u.cm)).value
dL_interpolation = interp1d(z_array, dl_cm)

l_box = 2000.0
L_box = l_box/h

env = "UCHUU"
Mmin = float(7.0e9) # 20 times the resolution
ID_in_0 = int(sys.argv[1]) #"/home/users/dae/comparat/Uchuu/hlists/halodir_050/halolist_z0p00_0.h5"

z_names = n.array([  "z0p00" #  0
					,"z0p02" #  1
					,"z0p05" #  2
					,"z0p09" #  3
					,"z0p14" #  4
					,"z0p19" #  5
					,"z0p25" #  6
					,"z0p30" #  7
					,"z0p36" #  8
					,"z0p43" #  9
					,"z0p49" # 10
					,"z0p56" # 11
					,"z0p63" # 12
					,"z0p70" # 13
					,"z0p78" # 14
					,"z0p86" # 15
					,"z0p94" # 16
					,"z1p03" # 17
					,"z1p12" # 18
					,"z1p22" # 19
					,"z1p32" # 20
					,"z1p43" # 21
					,"z1p54" # 22
					,"z1p65" # 23
					,"z1p77" # 24
					,"z1p90" # 25
					,"z2p03" # 26
					,"z2p17" # 27
					,"z2p31" # 28
					,"z2p46" # 29
					,"z2p62" # 30
					,"z2p78" # 31
					,"z2p95" # 32
					,"z3p13" # 33
					,"z3p32" # 34
					,"z3p61" # 35
					,"z3p93" # 36
					,"z4p27" # 37
					,"z4p63" # 38
					,"z5p15" # 39
					,"z5p73" # 40
					,"z6p35" # 41
					,"z7p03" #
					,"z7p76" #
					,"z8p59" #
					,"z9p47" #
					,"z10p44"
					,"z11p51"
					,"z12p69"
					,"z13p96" ])

z_name = z_names[ID_in_0]
id_snap = n.arange(len(z_names))[::-1][ID_in_0]
print(z_name, id_snap)

#N_snap_i = float(path_2_in_0.split('/')[-2][-3:])

# reads the list of snapshot available and the boundaries to be applied for each
# previously computed with 000_geometry.py
N_snap, Z_snap, A_snap, DC_max, DC_min, NN_replicas = n.loadtxt(os.path.join(os.environ[env], 'snap_list_with_border.txt'), unpack=True).T[id_snap].T
print(N_snap, Z_snap, A_snap, DC_max, DC_min, NN_replicas)

# Retrieve the snapshot and gets the boundaries
#sel = (N_snap == N_snap_i)
# Get the NN parameter: very important
# NN: number of times the box is replicated by 1 Lbox positively and negatively. Example: NN=3 , the box gets replicated by -3xL_box, -2xL_box, -1xL_box, 0xL_box, 1xL_box, 2xL_box in all X,Y,Z directions, (2*NN)^3 = 6^3 = 216 replications.
NN = NN_replicas
D2_min = DC_min**2
D2_max = DC_max**2

# specifyng the output directory 
LC_dir = os.path.join( os.environ[env], "GPX8" )
out_dir = os.path.join(LC_dir, z_name )
if os.path.isdir(out_dir) == False:
	os.system("mkdir -p " + out_dir)

path_2_in_1 = os.path.join(os.environ[env], 'Uchuu_UM', 'Uchuu_UM_'+z_name+'_data1.h5' )
path_2_in_2 = os.path.join(os.environ[env], 'Uchuu_UM', 'Uchuu_UM_'+z_name+'_data2.h5' )
path_2_in_3 = os.path.join(os.environ[env], 'Uchuu_UM', 'Uchuu_UM_'+z_name+'_data3.h5' )

# STELLAR MASS
# Equations 1 of Comparat et al. 2019
def meanSM(Mh, z): return n.log10(Mh * 2. * (0.0351 - 0.0247 * z / (1. + z)) / ((Mh / (10**(11.79 + 1.5 * z / (1. + z))))** (- 0.9 + 0.5 * z / (1. + z)) + (Mh / (10**(11.79 + 1.5 * z / (1. + z))))**(0.67 + 0.2 * z / (1. + z))))

# QUENCHED fraction
def scatter_Qf(z): return - 0.45 * (1 + z) + 1.54

def log10_M0_Qf(z): return 9.71 + 0.78 * (1 + z)

def fraction_Qf(mass, z): return 0.5 + 0.5 * erf((mass - log10_M0_Qf(z)) / scatter_Qf(z))

def beta_z(z): return -0.57 * z + 1.43

def alpha_z(z): return 6.32 * z - 16.26

# SFR
def mean_SFR_Q(mass, z): return mass * beta_z(z) + alpha_z(z)

def scale_z(z): return -0.34 * z + 0.99

# Hard X-ray emission, after Aird et al. 2018
def galaxy_lx(redshift, mass, sfr):
	return 10**(28.81) * (1 + redshift)**(3.9) * mass + \
        10**(39.5) * (1 + redshift)**(0.67) * sfr**(0.86)


# reading input file
hf1 = h5py.File( path_2_in_1, 'r')
hf2 = h5py.File( path_2_in_2, 'r')
hf3 = h5py.File( path_2_in_3, 'r')

# hf.keys()
#In [24]: hf1.keys()
#Out[24]: <KeysViewHDF5 ['Mvir', 'icl', 'id', 'obs_sfr', 'obs_sm', 'obs_uv', 'sfr', 'sm', 'upid']>
#In [25]: hf2.keys()
#Out[25]: <KeysViewHDF5 ['x', 'y', 'z']>
#In [26]: hf3.keys()
#Out[26]: <KeysViewHDF5 ['A_UV', 'Mpeak', 'Vmax_Mpeak', 'desc_id', 'lvmp', 'vmax', 'vx', 'vy', 'vz']>

# Creates the regular pavement for the replication
pts_i = n.arange(-1 * NN, NN, 1)
ix_i, iy_i, iz_i = n.meshgrid(pts_i, pts_i, pts_i)
ix, iy, iz = n.ravel(ix_i), n.ravel(iy_i), n.ravel(iz_i)
# print(ix,iy,iz)
pts = pts_i * L_box
x0_i, y0_i, z0_i = n.meshgrid(pts, pts, pts)
x0, y0, z0 = n.ravel(x0_i), n.ravel(y0_i), n.ravel(z0_i)
XT = n.transpose([x0, y0, z0])
iXT = n.transpose([ix, iy, iz])

######################
# REPLICATE over the sky
######################
mvir = n.array( hf1['Mvir'] )
x = n.array( hf2['x'] )/h
y = n.array( hf2['y'] )/h
z = n.array( hf2['z'] )/h
m_selection = ( mvir > Mmin )

# Loop over all (2*NN)^3 configuration to do all the replications
for jj in range(len(iXT)):
	jx,jy,jz=iXT[jj]
	out_dir2 = os.path.join(  out_dir, 'replication_'+str(jx)+'_'+str(jy)+'_'+str(jz))
	if os.path.isdir(out_dir2) == False:
		os.system("mkdir -p " + out_dir2)
	base_2_out = os.path.join(  out_dir2, 'glist')
	print(base_2_out)
	path_2_out = base_2_out + '.fits'
	if os.path.isfile(path_2_out)==False :
		xt = XT[jj]
		x0, y0, z0 = xt
		d2 = ( x + x0 )**2 + ( y + y0 )**2 + ( z + z0 )**2
		d_selection = ( d2 > D2_min ) & ( d2 < D2_max ) & ( m_selection )
		N_selected = len(x[d_selection])
		if N_selected>0:
			# now computing coordinates
			xn = x[d_selection] + x0
			yn = y[d_selection] + y0
			zn = z[d_selection] + z0
			rr = d2[d_selection]**0.5
			vxn = hf3['vx'][d_selection]
			vyn = hf3['vy'][d_selection]
			vzn = hf3['vz'][d_selection]
			# angular coordinates
			theta = n.arccos(zn / rr) * 180 / n.pi
			phi = n.arctan2(yn, xn) * 180 / n.pi
			ra = phi + 180.
			dec = theta - 90.
			# galactic and ecliptic coordinates
			coords = SkyCoord(ra, dec, unit='deg', frame='icrs')
			bb_gal = coords.galactic.b.value
			ll_gal = coords.galactic.l.value
			bb_ecl = coords.barycentrictrueecliptic.lat
			ll_ecl = coords.barycentrictrueecliptic.lon
			# extinction
			ebv = planck(coords)
			# line of sight, redshift
			redshift_R = dc_to_z(rr)
			vPara = (vxn * xn + vyn * yn + vzn * zn) / rr
			rr_s = rr + vPara / cosmo.H(redshift_R).value
			rr_s[rr_s <= 0] = rr[rr_s <= 0]
			redshift_S = dc_to_z(rr_s)
			# NH
			HEALPIX = healpy.ang2pix(1024, n.pi / 2. - bb_gal * n.pi / 180., 2 * n.pi - ll_gal * n.pi / 180.)
			NH = nh_val[HEALPIX]
			# final table
			t = []
			for kk in hf2.keys():
				if kk=='x':
					t.append(fits.Column(name=kk, array = xn, unit='Mpc', format=hf2[kk].dtype ))
				elif kk=='y':
					t.append(fits.Column(name=kk, array = yn, unit='Mpc', format=hf2[kk].dtype))
				elif kk=='z':
					t.append(fits.Column(name=kk, array = zn, unit='Mpc', format=hf2[kk].dtype))
			for kk in hf1.keys():
				t.append(fits.Column(name=kk, array=hf1[kk][d_selection], unit='', format=hf1[kk].dtype))
			for kk in hf3.keys():
				print(kk, hf3[kk].dtype, hf3[kk][d_selection])
				t.append(fits.Column(name=kk, array=hf3[kk][d_selection], unit='', format=hf3[kk].dtype))

			# luminosity distance
			dL_cm = dL_interpolation(redshift_R)
			# coordinate columns
			t.append(fits.Column(name='RA', array=ra, unit='deg', format = 'D'))
			t.append(fits.Column(name='DEC', array=dec, unit='deg', format = 'D'))
			#
			t.append(fits.Column(name='g_lat', array=bb_gal, unit='deg', format = 'D'))
			t.append(fits.Column(name='g_lon', array=ll_gal, unit='deg', format = 'D'))
			#
			t.append(fits.Column(name='ecl_lat', array=bb_ecl, unit='deg', format = 'D'))
			t.append(fits.Column(name='ecl_lon', array=ll_ecl, unit='deg', format = 'D'))
			#
			t.append(fits.Column(name='redshift_R', array=redshift_R, unit='', format = 'D'))
			t.append(fits.Column(name='redshift_S', array=redshift_S, unit='', format = 'D'))
			#
			t.append(fits.Column(name='dL', array=dL_cm, unit='cm', format = 'D'))
			#
			t.append(fits.Column(name='nH', array=NH, unit='cm**(-2)', format = 'D'))
			t.append(fits.Column(name='ebv', array=ebv, unit='mag', format = 'D'))
			# writing
			print('writing', path_2_out)
			tbhdu = fits.BinTableHDU.from_columns(fits.ColDefs(t))
			tbhdu.header['author'] = 'JC'
			if os.path.isfile(path_2_out):
				os.remove(path_2_out)
			tbhdu.writeto(path_2_out)

