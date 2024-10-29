"""

Replicates the snapshot over the sky and extracts halos in a given shell of comoving distance
DONE
#nohup python 000_process_uchuu_DMH_pixel_i.py 00 > log_batch_00.log &
#nohup python 000_process_uchuu_DMH_pixel_i.py 01 > log_batch_01.log &
#nohup python 000_process_uchuu_DMH_pixel_i.py 02 > log_batch_02.log &
#nohup python 000_process_uchuu_DMH_pixel_i.py 03 > log_batch_03.log &
#nohup python 000_process_uchuu_DMH_pixel_i.py 04 > log_batch_04.log &
#nohup python 000_process_uchuu_DMH_pixel_i.py 05 > log_batch_05.log &
DONE
ds43
#nohup python 000_process_uchuu_DMH_pixel_i.py 06 > log_batch_06.log &
#nohup python 000_process_uchuu_DMH_pixel_i.py 07 > log_batch_07.log &
#nohup python 000_process_uchuu_DMH_pixel_i.py 08 > log_batch_08.log &
#nohup python 000_process_uchuu_DMH_pixel_i.py 09 > log_batch_09.log &
#nohup python 000_process_uchuu_DMH_pixel_i.py 10 > log_batch_10.log &
#nohup python 000_process_uchuu_DMH_pixel_i.py 11 > log_batch_11.log &
ds54
#nohup python 000_process_uchuu_DMH_pixel_i.py 12 > log_batch_12.log &
#nohup python 000_process_uchuu_DMH_pixel_i.py 13 > log_batch_13.log &
#nohup python 000_process_uchuu_DMH_pixel_i.py 14 > log_batch_14.log &
#nohup python 000_process_uchuu_DMH_pixel_i.py 15 > log_batch_15.log &
ds52
#nohup python 000_process_uchuu_DMH_pixel_i.py 16 > log_batch_16.log &
#nohup python 000_process_uchuu_DMH_pixel_i.py 17 > log_batch_17.log &
#nohup python 000_process_uchuu_DMH_pixel_i.py 18 > log_batch_18.log &
#nohup python 000_process_uchuu_DMH_pixel_i.py 19 > log_batch_19.log &
#nohup python 000_process_uchuu_DMH_pixel_i.py 20 > log_batch_20.log &
#nohup python 000_process_uchuu_DMH_pixel_i.py 21 > log_batch_21.log &
#nohup python 000_process_uchuu_DMH_pixel_i.py 22 > log_batch_22.log &
#nohup python 000_process_uchuu_DMH_pixel_i.py 23 > log_batch_23.log &
#nohup python 000_process_uchuu_DMH_pixel_i.py 24 > log_batch_24.log &
#nohup python 000_process_uchuu_DMH_pixel_i.py 25 > log_batch_25.log &
#nohup python 000_process_uchuu_DMH_pixel_i.py 26 > log_batch_26.log &

"""
import time
t0=time.time()
import numpy as np
import os, sys, glob
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

cosmoUNIT = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
h = 0.6774
cosmo = cosmoUNIT

z_array = np.arange(0, 7.5, 0.001)
dc_to_z = interp1d(cosmo.comoving_distance(z_array), z_array)

l_box = 2000.0
L_box = l_box/h

env = "UCHUU"
Mmin = 9.0 
ID_in_0 = int(sys.argv[1]) #"/home/users/dae/comparat/Uchuu/hlists/halodir_050/halolist_z0p00_0.h5"

z_names = np.array([  "z0p00" #  0     "halodir_050",
                    ,"z0p02" #  1    "halodir_049",
                    ,"z0p05" #  2    "halodir_048",
                    ,"z0p09" #  3    "halodir_047",
                    ,"z0p14" #  4    "halodir_046",
                    ,"z0p19" #  5    "halodir_045",
                    ,"z0p25" #  6    "halodir_044",
                    ,"z0p30" #  7    "halodir_043",
                    ,"z0p36" #  8    "halodir_042",
                    ,"z0p43" #  9    "halodir_041",
                    ,"z0p49" # 10    "halodir_040",
                    ,"z0p56" # 11    "halodir_039",
                    ,"z0p63" # 12    "halodir_038",
                    ,"z0p70" # 13    "halodir_037",
                    ,"z0p78" # 14    "halodir_036",
                    ,"z0p86" # 15    "halodir_035",
                    ,"z0p94" # 16    "halodir_034",
                    ,"z1p03" # 17    "halodir_033",
                    ,"z1p12" # 18    "halodir_032",
                    ,"z1p22" # 19    "halodir_031",
                    ,"z1p32" # 20    "halodir_030",
                    ,"z1p43" # 21    "halodir_029",
                    ,"z1p54" # 22    "halodir_028",
                    ,"z1p65" # 23    "halodir_027",
                    ,"z1p77" # 24    "halodir_026",
                    ,"z1p90" # 25    "halodir_025",
                    ,"z2p03" # 26    "halodir_024"
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

H_names = np.array([  "halodir_050", #   "z0p00" #  0   
                    "halodir_049",  #  ,"z0p02"     1   
                    "halodir_048",  #  ,"z0p05"     2   
                    "halodir_047",  #  ,"z0p09"     3   
                    "halodir_046",  #  ,"z0p14"     4   
                    "halodir_045",  #  ,"z0p19"     5   
                    "halodir_044",  #  ,"z0p25"     6   
                    "halodir_043",  #  ,"z0p30"     7   
                    "halodir_042",  #  ,"z0p36"     8   
                    "halodir_041",  #  ,"z0p43"     9   
                    "halodir_040",  #  ,"z0p49"    10   
                    "halodir_039",  #  ,"z0p56"    11   
                    "halodir_038",  #  ,"z0p63"    12   
                    "halodir_037",  #  ,"z0p70"    13   
                    "halodir_036",  #  ,"z0p78"    14   
                    "halodir_035",  #  ,"z0p86"    15   
                    "halodir_034",  #  ,"z0p94"    16   
                    "halodir_033",  #  ,"z1p03"    17   
                    "halodir_032",  #  ,"z1p12"    18   
                    "halodir_031",  #  ,"z1p22"    19   
                    "halodir_030",  #  ,"z1p32"    20   
                    "halodir_029",  #  ,"z1p43"    21   
                    "halodir_028",  #  ,"z1p54"    22   
                    "halodir_027",  #  ,"z1p65"    23   
                    "halodir_026",  #  ,"z1p77"    24   
                    "halodir_025",  #  ,"z1p90"    25   
                    "halodir_024"   #  ,"z2p03"    26   
                    ])


z_name = z_names[ID_in_0]
halo_dir = H_names[ID_in_0]
id_snap = np.arange(len(z_names))[::-1][ID_in_0]
print(z_name, id_snap)

# reads the list of snapshot available and the boundaries to be applied for each
# previously computed with 000_geometry.py
N_snap, Z_snap, A_snap, DC_max, DC_min, NN_replicas = np.loadtxt('snap_list_with_border.txt', unpack=True).T[id_snap].T
print(N_snap, Z_snap, A_snap, DC_max, DC_min, NN_replicas)

# Retrieve the snapshot and gets the boundaries
#sel = (N_snap == N_snap_i)
# Get the NN parameter: very important
# NN: number of times the box is replicated by 1 Lbox positively and negatively. Example: NN=3 , the box gets replicated by -3xL_box, -2xL_box, -1xL_box, 0xL_box, 1xL_box, 2xL_box in all X,Y,Z directions, (2*NN)^3 = 6^3 = 216 replications.
NN = NN_replicas
D2_min = DC_min**2
D2_max = DC_max**2

# specifyng the output directory 
LC_dir = os.path.join( os.environ[env], "HaloSky" )
out_dir = os.path.join(LC_dir, z_name )
if os.path.isdir(out_dir) == False:
	os.system("mkdir -p " + out_dir)

path_2_in_files = np.array( glob.glob( os.path.join(os.environ[env], 'RockstarSmall', halo_dir, 'halolist_'+z_name+'_*.fits' ) ) )
path_2_in_files.sort()
for path_2_in in path_2_in_files:
    file_Bname = os.path.basename(path_2_in)[:-5]
    file_N = file_Bname.split('_')[-1]
    # reading input file
    hf1 = Table.read( path_2_in )
    # Creates the regular pavement for the replication
    pts_i = np.arange(-1 * NN, NN, 1)
    ix_i, iy_i, iz_i = np.meshgrid(pts_i, pts_i, pts_i)
    ix, iy, iz = np.ravel(ix_i), np.ravel(iy_i), np.ravel(iz_i)
    # print(ix,iy,iz)
    pts = pts_i * L_box
    x0_i, y0_i, z0_i = np.meshgrid(pts, pts, pts)
    x0, y0, z0 = np.ravel(x0_i), np.ravel(y0_i), np.ravel(z0_i)
    XT = np.transpose([x0, y0, z0])
    iXT = np.transpose([ix, iy, iz])

    ######################
    # REPLICATE over the sky
    ######################
    x = np.array( hf1['x'] )/h
    y = np.array( hf1['y'] )/h
    z = np.array( hf1['z'] )/h

    # Loop over all (2*NN)^3 configuration to do all the replications
    for jj in range(len(iXT)):
        jx,jy,jz=iXT[jj]
        out_dir2 = os.path.join(  out_dir, 'replication_'+str(jx)+'_'+str(jy)+'_'+str(jz))
        if os.path.isdir(out_dir2) == False:
                os.system("mkdir -p " + out_dir2)
        base_2_out = os.path.join(  out_dir2, 'hlist_'+file_N)
        print(base_2_out)
        path_2_out = base_2_out + '.fits'
        if os.path.isfile(path_2_out)==False :
            xt = XT[jj]
            x0, y0, z0 = xt
            d2 = ( x + x0 )**2 + ( y + y0 )**2 + ( z + z0 )**2
            d_selection = ( d2 >= D2_min ) & ( d2 < D2_max ) 
            N_selected = len(x[d_selection])
            if N_selected>0:
                print(N_selected)
                # now computing coordinates
                xn = x[d_selection] + x0
                yn = y[d_selection] + y0
                zn = z[d_selection] + z0
                rr = d2[d_selection]**0.5
                vxn = hf1['vx'][d_selection]
                vyn = hf1['vy'][d_selection]
                vzn = hf1['vz'][d_selection]
                # angular coordinates
                theta = np.arccos(zn / rr) * 180 / np.pi
                phi = np.arctan2(yn, xn) * 180 / np.pi
                ra = phi + 180.
                dec = theta - 90.
                # galactic and ecliptic coordinates
                coords = SkyCoord(ra, dec, unit='deg', frame='icrs')
                bb_gal = coords.galactic.b.value
                ll_gal = coords.galactic.l.value
                # extinction
                #ebv = planck(coords)
                # line of sight, redshift
                redshift_R = dc_to_z(rr)
                vPara = (vxn * xn + vyn * yn + vzn * zn) / rr
                rr_s = rr + vPara / cosmo.H(redshift_R).value
                rr_s[rr_s <= 0] = rr[rr_s <= 0]
                redshift_S = dc_to_z(rr_s)
                # final table
                t = hf1[d_selection]
                t['x'] = xn
                t['x'].unit='Mpc'
                t['x'].description='Halo position (Mpc comoving)'   
                t['y'] = yn
                t['y'].unit='Mpc'
                t['y'].description='Halo position (Mpc comoving)'   
                t['z'] = zn
                t['z'].unit='Mpc'
                t['z'].description='Halo position (Mpc comoving)'   

                # coordinate columns
                t['RA']  =ra
                t['DEC'] =dec
                t['RA']  .unit='deg'
                t['DEC'] .unit='deg'

                t['g_lat']=bb_gal
                t['g_lon']=ll_gal
                t['g_lat'].unit='deg'
                t['g_lon'].unit='deg'

                t['redshift_R']=redshift_R
                t['redshift_S']=redshift_S
    
                # writing
                print('writing', path_2_out)
                t.write(path_2_out, overwrite = True)
            
        

