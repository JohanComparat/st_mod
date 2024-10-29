"""
What it does
------------

Computes the Host Mvir for the galaxies in the light cone.
 * for a central galaxy: CENTRAL_Mvir = Mvir
 * for a satellite galaxy, CENTRAL_Mvir = Mvir of the host halo found by matching upid with id

conda activate stmod
cd $GIT_STMOD/src/lightcones/uchuu
nohup python compute_HOD.py z0p00 > logs/log_computeHOD_uchuu_z0p00.log & # DONE
nohup python compute_HOD.py z0p02 > logs/log_computeHOD_uchuu_z0p02.log & # DONE
nohup python compute_HOD.py z0p05 > logs/log_computeHOD_uchuu_z0p05.log & # DONE
nohup python compute_HOD.py z0p09 > logs/log_computeHOD_uchuu_z0p09.log & # DONE
nohup python compute_HOD.py z0p14 > logs/log_computeHOD_uchuu_z0p14.log & # DONE
nohup python compute_HOD.py z0p19 > logs/log_computeHOD_uchuu_z0p19.log & # DONE
nohup python compute_HOD.py z0p25 > logs/log_computeHOD_uchuu_z0p25.log & # DONE
nohup python compute_HOD.py z0p30 > logs/log_computeHOD_uchuu_z0p30.log & # DONE
nohup python compute_HOD.py z0p36 > logs/log_computeHOD_uchuu_z0p36.log & # TODO
nohup python compute_HOD.py z0p43 > logs/log_computeHOD_uchuu_z0p43.log & # TODO
nohup python compute_HOD.py z0p49 > logs/log_computeHOD_uchuu_z0p49.log & # TODO
nohup python compute_HOD.py z0p56 > logs/log_computeHOD_uchuu_z0p56.log & # TODO

"""
import sys, glob, os, time
from astropy.table import Table, Column
import numpy as n
import astropy.io.fits as fits

print('Add central mass ')
print('------------------------------------------------')
print('------------------------------------------------')
z_dir = sys.argv[1] # "z0p05"
env = 'UCHUU'
p_2_catalogues = n.array( glob.glob( os.path.join(os.environ[env], 'GPX8', z_dir, 'replication_*_*_*', 'glist.fits') ) )

def write_mCEN(p_2_catalogue, p_2_catalogue_out, path_2_HistOut_file,path_2_HistSatOut_file):
	# input catalogue
	t_in = fits.open(p_2_catalogue)
	upid = t_in[1].data['upid']
	zz_1 = t_in[1].data['redshift_R']
	hid = t_in[1].data['id']
	M_in = n.log10(t_in[1].data['Mvir'])
	N_gals = len(M_in)
	#
	cen = (upid==-1)
	sat = (cen==False)
	all_ids = n.arange(len(zz_1))
	# sat for which the host halo is in the light cone (>99.5%)
	is_avail = n.isin(upid, hid)
	sat2 = (sat) & (is_avail)
	# table to be output
	t_out = Table.read(p_2_catalogue_out)

	# histogram2d of redshift and mass
	H_2d = n.histogram2d( t_out['CENTRAL_Mvir'][cen], zz_1[cen], bins=[n.arange(9.0, 16.0, 0.05), n.arange(0,6.2, 0.01)])[0]
	H_2dS = n.histogram2d( t_out['CENTRAL_Mvir'][sat2], zz_1[sat2], bins=[n.arange(9.0, 16.0, 0.05), n.arange(0,6.2, 0.01)])[0]
	n.savetxt(path_2_HistOut_file, H_2d)
	n.savetxt(path_2_HistSatOut_file, H_2dS)


for p_2_catalogue in p_2_catalogues:
	p_2_catalogue_out = os.path.join( os.path.dirname(p_2_catalogue), 'mcen_list.fits')
	path_2_HistOut_file = os.path.join( os.path.dirname(p_2_catalogue), 'HistCen_Mvir_zz.ascii')
	path_2_HistSatOut_file = os.path.join( os.path.dirname(p_2_catalogue), 'HistSat_HostMvir_zz.ascii' )
	print(p_2_catalogue, p_2_catalogue_out, path_2_HistOut_file, path_2_HistSatOut_file)
	write_mCEN(p_2_catalogue, p_2_catalogue_out, path_2_HistOut_file, path_2_HistSatOut_file)


