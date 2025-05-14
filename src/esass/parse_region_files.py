import time
t0=time.time()
import os
import numpy as n

from astropy.table import Table, Column

sky_map_hdu = Table.read(os.path.join(os.environ['GIT_ERASS_SIM'], 'data', 'SKYMAPS.fits'))

ttt = n.arange(len(sky_map_hdu[(sky_map_hdu['OWNER'] == 2) | (sky_map_hdu['OWNER'] == 0)]))
sky_tile_id = sky_map_hdu['SRVMAP'][(sky_map_hdu['OWNER'] == 2) | (sky_map_hdu['OWNER'] == 0)]

# all_div = divisors(len(ttt))
# tRS = ttt.reshape((int(len(ttt)/271),271))
for ii in n.arange(len(ttt)):
	tile_id = str(sky_tile_id[ii]).zfill(6)
	#inputdir = "/data26s/comparat/simulations/erosim/eRASS1_UNIT_fA1i_2021_10_12_SKYMAP/"+tile_id
	inputdir = os.path.join("/home/idies/workspace/erosim/Uchuu/LCerass/", tile_id, 'sim_evt_e4_merge', 'eSASS')
	SrcReg = os.path.join(inputdir, tile_id + "__02_23_srcAUTO.reg")
	SrcRad = os.path.join(inputdir, tile_id + "_srcRAD.fits")
	print(tile_id, os.path.isfile(SrcReg), os.path.isfile(SrcRad))
	if os.path.isfile(SrcReg) and os.path.isfile(SrcRad) == False:
		f=open(SrcReg, 'r')
		out = f.readlines()
		N_entries = len(out)
		R0, Rmin, Rmax = n.zeros(N_entries), n.zeros(N_entries), n.zeros(N_entries)
		for jj in n.arange(len(out)):
			circles = out[jj].split(';')[1::2]
			radii = n.array([float(el.split(' ')[-1]) for el in circles ])
			R0   [jj] = radii[0]
			Rmin [jj] = radii.min()
			Rmax [jj] = radii.max()

		t = Table()
		t.add_column(Column(name='R0', data=R0, unit='deg', dtype=n.float32 ) )
		t.add_column(Column(name='Rmin', data=Rmin, unit='deg', dtype=n.float32 ) )
		t.add_column(Column(name='Rmax', data=Rmax, unit='deg', dtype=n.float32 ) )
		t.write(SrcRad, overwrite=True)
		print(SrcRad, N_entries)#, (time.time()-t0)/N_entries)
