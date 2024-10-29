import time
t0=time.time()
import numpy as n
import os, sys, glob
print(sys.argv)
import h5py
from astropy.table import Table, Column, hstack, vstack

#h_dir = 'halodir_050'
h_dir = sys.argv[1]
Mvir_min = 11.0

out_dir = os.path.join( os.environ['HOME'], 'data7/UCHUU/RockstarSmall', h_dir)
os.system('mkdir -p '+out_dir)
ext_dir = os.path.join( os.environ['HOME'], 'data7/UCHUU/RockstarExtended', h_dir)
all_in = n.array(glob.glob(os.path.join(ext_dir, "halolist_*.h5" )))
col_list = ['id', 'pid', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'M200c', 'M500c', 'Mvir', 'Mvir_all', 'Rvir','rs', 'Spin','Spin_Bullock','T_U', 'Xoff', 'b_to_a', 'c_to_a','scale_of_last_MM']

all_in.sort()
for path_2_in_0 in all_in :
    print(path_2_in_0)
    #p_2_out = os.path.join( out_dir, os.path.basename(path_2_in_0) )
    p_2_out = os.path.join( out_dir, os.path.basename(path_2_in_0)[:-2]+'fits' )
    if os.path.isfile(p_2_out)==False:
        #print(p_2_out)
        hf = h5py.File(path_2_in_0, "r")
        selection = ( n.log10(n.array(hf['Mvir'    ])) > Mvir_min )
        t = Table()
        t.add_column(Column(name='id' , data=n.array(hf['id']) [selection], unit='', dtype=n.int64, description='ID of halo (unique across entire simulation)') )
        t.add_column(Column(name='pid', data=n.array(hf['pid'])[selection], unit='', dtype=n.int64, description='ID of least massive host halo (-1 if distinct halo)') )

        t.add_column(Column(name='x' , data=n.array(hf['x']) [selection], unit='Mpc/h', dtype=n.float32, description='Halo position (Mpc/h comoving)' ) )
        t.add_column(Column(name='y' , data=n.array(hf['y']) [selection], unit='Mpc/h', dtype=n.float32, description='Halo position (Mpc/h comoving)' ) )
        t.add_column(Column(name='z' , data=n.array(hf['z']) [selection], unit='Mpc/h', dtype=n.float32, description='Halo position (Mpc/h comoving)' ) )
        t.add_column(Column(name='vx', data=n.array(hf['vx'])[selection], unit='km/s', dtype=n.float32, description='Halo velocity (km/s physical)' ) )
        t.add_column(Column(name='vy', data=n.array(hf['vy'])[selection], unit='km/s', dtype=n.float32, description='Halo velocity (km/s physical)' ) )
        t.add_column(Column(name='vz', data=n.array(hf['vz'])[selection], unit='km/s', dtype=n.float32, description='Halo velocity (km/s physical)' ) )

        t.add_column(Column(name='M200c'   , data=n.log10(n.array(hf['M200c'   ]))[selection], unit='log10(M[Msun/h])', dtype=n.float32, description='Mass enclosed within 200c overdensity') )
        t.add_column(Column(name='M500c'   , data=n.log10(n.array(hf['M500c'   ]))[selection], unit='log10(M[Msun/h])', dtype=n.float32, description='Mass enclosed within 500c overdensity') )
        t.add_column(Column(name='Mvir'    , data=n.log10(n.array(hf['Mvir'    ]))[selection], unit='log10(M[Msun/h])', dtype=n.float32, description='Mass enclosed within vir overdensity') )
        t.add_column(Column(name='Mvir_all', data=n.log10(n.array(hf['Mvir_all']))[selection], unit='log10(M[Msun/h])', dtype=n.float32, description='Mass enclosed within vir overdensity, including unbound particles') )

        t.add_column(Column(name='Rvir', data=n.array(hf['Rvir'])[selection], unit='kpc/h', dtype=n.float32, description='Halo radius (kpc/h comoving)') )
        t.add_column(Column(name='rs'  , data=n.array(hf['rs'])  [selection], unit='kpc/h', dtype=n.float32, description='Scale radius (kpc/h comoving)') )
        t.add_column(Column(name='Spin', data=n.array(hf['Spin'])[selection], unit='', dtype=n.float32, description='Halo spin parameter') )
        t.add_column(Column(name='Spin_Bullock', data=n.array(hf['Spin_Bullock'])[selection], unit='', dtype=n.float32, description='Bullock spin parameter (J/(sqrt(2)*GMVR))') )
        t.add_column(Column(name='Xoff'  , data=n.array(hf['Xoff'])  [selection], unit='kpc/h', dtype=n.float32, description='Offset of density peak from average particle position (kpc/h comoving)') )
        t.add_column(Column(name='b_to_a', data=n.array(hf['b_to_a'])[selection], unit='', dtype=n.float32, description='Ratio of second largest shape ellipsoid axes (B) to largest shape ellipsoid axis (A) (dimensionless)') )
        t.add_column(Column(name='c_to_a', data=n.array(hf['c_to_a'])[selection], unit='', dtype=n.float32, description='Ratio of third largest shape ellipsoid axes (C) to largest shape ellipsoid axis (A) (dimensionless)') )
        t.add_column(Column(name='scale_of_last_MM', data=n.array(hf['scale_of_last_MM'])[selection], unit='', dtype=n.float32, description='scale factor of the last major merger (Mass ratio > 0.3)') )
        t.add_column(Column(name='T_U'   , data=n.array(hf['T_U'])[selection], unit='', dtype=n.float32, description='ratio of kinetic to potential energies') )
        #t.write(p_2_out, overwrite=True)
        #print(p_2_out,'written')
        t.write(p_2_out, overwrite=True)
        print(p_2_out,'written',time.time()-t0,'s')
    else:
        print('done')

#ID: ID of halo (unique across entire simulation).
#Desc_Scale: Scale of descendant halo, if applicable.
#Descid: ID of descendant halo, if applicable.
#Num_prog: Number of progenitors.
#Pid: ID of least massive host halo (-1 if distinct halo).
#Upid: ID of most massive host halo (different from Pid when the halo is within two or more larger halos).
#Desc_pid: Pid of descendant halo (if applicable).
#Phantom: Nonzero for halos interpolated across timesteps.
#SAM_Mvir: Halo mass, smoothed across accretion history; always greater than sum of halo masses of contributing progenitors (Msun/h).  Only for use with select semi-analytic models.
#Mvir: Halo mass (Msun/h).
#Rvir: Halo radius (kpc/h comoving).
#Rs: Scale radius (kpc/h comoving).
#Vrms: Velocity dispersion (km/s physical).
#mmp?: whether the halo is the most massive progenitor or not.
#scale_of_last_MM: scale factor of the last major merger (Mass ratio > 0.3).
#Vmax: Maxmimum circular velocity (km/s physical).
#X/Y/Z: Halo position (Mpc/h comoving).
#VX/VY/VZ: Halo velocity (km/s physical).
#JX/JY/JZ: Halo angular momenta ((Msun/h) * (Mpc/h) * km/s (physical)).
#Spin: Halo spin parameter.
#Breadth_first_ID: breadth-first ordering of halos within a tree.
#Depth_first_ID: depth-first ordering of halos within a tree.
#Tree_root_ID: ID of the halo at the last timestep in the tree.
#Orig_halo_ID: Original halo ID from halo finder.
#Snap_num: Snapshot number from which halo originated.
#Next_coprogenitor_depthfirst_ID: Depthfirst ID of next coprogenitor.
#Last_progenitor_depthfirst_ID: Depthfirst ID of last progenitor.
#Last_mainleaf_depthfirst_ID: Depthfirst ID of last progenitor on main progenitor branch.
#Tidal_Force: Strongest tidal force from any nearby halo, in dimensionless units (Rhalo / Rhill).
#Tidal_ID: ID of halo exerting strongest tidal force.
#Rs_Klypin: Scale radius determined using Vmax and Mvir (see Rockstar paper)
#Mvir_all: Mass enclosed within the specified overdensity, including unbound particles (Msun/h)
#M200b--M2500c: Mass enclosed within specified overdensities (Msun/h)
#Xoff: Offset of density peak from average particle position (kpc/h comoving)
#Voff: Offset of density peak from average particle velocity (km/s physical)
#Spin_Bullock: Bullock spin parameter (J/(sqrt(2)*GMVR))
#b_to_a, c_to_a: Ratio of second and third largest shape ellipsoid axes (B and C) to largest shape ellipsoid axis (A) (dimensionless).
#  Shapes are determined by the method in Allgood et al. (2006).
#  (500c) indicates that only particles within R500c are considered.
#A[x],A[y],A[z]: Largest shape ellipsoid axis (kpc/h comoving)
#T/|U|: ratio of kinetic to potential energies
#M_pe_*: Pseudo-evolution corrected masses (very experimental)
#Consistent Trees Version 1.01
#Includes fix for Rockstar spins & T/|U| (assuming T/|U| = column 53)
#Macc,Vacc: Mass and Vmax at accretion.
#Mpeak,Vpeak: Peak mass and Vmax over mass accretion history.
#Halfmass_Scale: Scale factor at which the MMP reaches 0.5*Mpeak.
#Acc_Rate_*: Halo mass accretion rates in Msun/h/yr.
#            Inst: instantaneous; 100Myr: averaged over past 100Myr,
#            X*Tdyn: averaged over past X*virial dynamical time.
#            Mpeak: Growth Rate of Mpeak, averaged from current z to z+0.5
#Mpeak_Scale: Scale at which Mpeak was reached.
#Acc_Scale: Scale at which satellites were (last) accreted.
#First_Acc_Scale: Scale at which current and former satellites first passed through a larger halo.
#First_Acc_(Mvir|Vmax): Mvir and Vmax at First_Acc_Scale.
#Vmax@Mpeak: Halo Vmax at the scale at which Mpeak was reached.
#Tidal_Force_Tdyn: Dimensionless tidal force averaged over past dynamical time.
