"""
What it does
------------

Computes the Host Mvir for the galaxies in the light cone.
 * for a central galaxy: CENTRAL_Mvir = Mvir
 * for a satellite galaxy, CENTRAL_Mvir = Mvir of the host halo found by matching upid with id

conda activate stmod
cd $GIT_STMOD/src/lightcones/uchuu
nohup python add_central_mvir.py z0p00 > logs/log_add_CENTRAL_mvir_uchuu_z0p00.log & # DONE
nohup python add_central_mvir.py z0p02 > logs/log_add_CENTRAL_mvir_uchuu_z0p02.log & # DONE
nohup python add_central_mvir.py z0p05 > logs/log_add_CENTRAL_mvir_uchuu_z0p05.log & # DONE
nohup python add_central_mvir.py z0p09 > logs/log_add_CENTRAL_mvir_uchuu_z0p09.log & # DONE
nohup python add_central_mvir.py z0p14 > logs/log_add_CENTRAL_mvir_uchuu_z0p14.log & # DONE
nohup python add_central_mvir.py z0p19 > logs/log_add_CENTRAL_mvir_uchuu_z0p19.log & # DONE
nohup python add_central_mvir.py z0p25 > logs/log_add_CENTRAL_mvir_uchuu_z0p25.log & # DONE
nohup python add_central_mvir.py z0p30 > logs/log_add_CENTRAL_mvir_uchuu_z0p30.log & # DONE
nohup python add_central_mvir.py z0p36 > logs/log_add_CENTRAL_mvir_uchuu_z0p36.log & # raven06
nohup python add_central_mvir.py z0p43 > logs/log_add_CENTRAL_mvir_uchuu_z0p43.log & # raven06
nohup python add_central_mvir.py z0p49 > logs/log_add_CENTRAL_mvir_uchuu_z0p49.log & # raven06
nohup python add_central_mvir.py z0p56 > logs/log_add_CENTRAL_mvir_uchuu_z0p56.log & # raven06

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
	set_hosts = upid[sat2]
	selection_hosts = n.isin(hid, n.unique(set_hosts))
	host_ids = all_ids[selection_hosts]
	all_mvir = M_in[host_ids]
	all_halo_id = hid[host_ids]

	host_mvir_i = n.array([all_mvir[all_halo_id==element][0] for element in upid[sat2] ])
	host_mvir = M_in
	host_mvir[sat2] = host_mvir_i
	# table to be output
	t_out = Table()
	t_out['CENTRAL_Mvir'] = host_mvir
	t_out.write(p_2_catalogue_out, overwrite=True)
	print(p_2_catalogue_out, 'written')

	# histogram2d of redshift and mass
	H_2d = n.histogram2d( host_mvir[cen], zz_1[cen], bins=[n.arange(9.0, 15.5, 0.1), n.arange(0,6.2, 0.1)])[0]
	H_2dS = n.histogram2d( host_mvir[sat2], zz_1[sat2], bins=[n.arange(9.0, 15.5, 0.1), n.arange(0,6.2, 0.1)])[0]
	n.savetxt(path_2_HistOut_file, H_2d)
	n.savetxt(path_2_HistSatOut_file, H_2dS)


for p_2_catalogue in p_2_catalogues:
	p_2_catalogue_out = os.path.join( os.path.dirname(p_2_catalogue), 'mcen_list.fits')
	path_2_HistOut_file = os.path.join( os.path.dirname(p_2_catalogue), 'HistCen_Mvir_zz.ascii')
	path_2_HistSatOut_file = os.path.join( os.path.dirname(p_2_catalogue), 'HistSat_HostMvir_zz.ascii' )
	print(p_2_catalogue, p_2_catalogue_out, path_2_HistOut_file, path_2_HistSatOut_file)
	#if os.path.isfile(p_2_catalogue_out)==False:
	write_mCEN(p_2_catalogue, p_2_catalogue_out, path_2_HistOut_file, path_2_HistSatOut_file)

