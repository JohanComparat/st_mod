"""
HAM procedure for Uchuu between tabulated AGN files and the light cone


# testing commands
conda activate stmod
cd $GIT_STMOD/src/Mpipelines

python AGN_pipeline.py z0p30 LC0002  0.8 8
python AGN_pipeline.py z0p30 LC0060  0.8 8

python AGN_validation.py z0p30 LC0002
python AGN_validation.py z0p30 LC0060

python AGN_pipeline.py z1p03 LC0002  0.8 8
python AGN_pipeline.py z1p03 LC0060  0.8 8

python AGN_validation.py z1p03 LC0002
python AGN_validation.py z1p03 LC0060

python AGN_validation.py z0p94 LC0002
python AGN_validation.py z0p94 LC0060

cd $UCHUU/LC0060/

"""

import time
t0 = time.time()
import sys, os, glob
from astropy.table import Table, vstack
import numpy as np

sys.path.append( os.path.join(os.environ['GIT_STMOD'], 'src') )
from models import AGN as GG

z_dir = sys.argv[1]
LC_dir = sys.argv[2] # 'FullSky'
scatter_0 = float(sys.argv[3])
f_sat = float(sys.argv[4])

C_AGN = GG.AGN(z_dir, LC_dir=LC_dir, scatter_0 = scatter_0, f_sat = f_sat )

DZ = 0.01
z_bins = np.arange(0.0, 6., DZ)

print(len(C_AGN.p_2_catalogues), 'catalogues to loop over', ', Dt=', time.time()-t0, 'seconds')
is_cat = np.array([os.path.isfile(el) for el in C_AGN.p_2_catalogues])
for p_2_catalogue in C_AGN.p_2_catalogues[is_cat]:
    print('='*100)
    ##
    ## retrieve the area subtended by the catalog, skip step if area=0
    ##
    print(p_2_catalogue)
    replication_dir = p_2_catalogue.split('/')[-2]
    print(replication_dir)
    path_2_output = os.path.join( os.path.dirname(p_2_catalogue),'AGN_list_sigma_'+C_AGN.str_scatter_0+'_fsat_'+C_AGN.str_fsat+'.fits' )
    print('already done ? ', os.path.isfile( path_2_output ))
    #if os.path.isfile( path_2_output )==False:
    ix = float(replication_dir.split('_')[1])
    iy = float(replication_dir.split('_')[2])
    iz = float(replication_dir.split('_')[3])
    #print(ix, iy, iz)
    selection_area = ( C_AGN.LC_MetaData['jx'] == ix ) & ( C_AGN.LC_MetaData['jy'] == iy ) & ( C_AGN.LC_MetaData['jz'] == iz )
    #print(C_AGN.LC_MetaData[selection_area])
    sky_frac_Dmin = C_AGN.LC_MetaData['area_DC_min'][selection_area].sum() / C_AGN.LC_MetaData['area_DC_min'].sum()
    sky_frac_Dmax = C_AGN.LC_MetaData['area_DC_max'][selection_area].sum() / C_AGN.LC_MetaData['area_DC_max'].sum()
    C_AGN.sky_frac = sky_frac_Dmax * C_AGN.LC_MetaData['N_frac_'+LC_dir][selection_area]
    print('sky fraction', C_AGN.sky_frac)
    if C_AGN.sky_frac==0.:
        print('area = 0 deg2')
        continue
    print('area = ', C_AGN.sky_frac, ' deg2')

    ##
    ## OPENS Galaxy light cone
    ##
    print('opens', p_2_catalogue)
    t_sim_gal_full = Table.read(p_2_catalogue)
    t_sim_gal_full['ID_glist'] = np.arange(len(t_sim_gal_full))
    ZZ = t_sim_gal_full['redshift_R']
    C_AGN.z_min = np.min(ZZ)
    C_AGN.z_max = np.max(ZZ)
    C_AGN.z_mean = np.mean(ZZ)
    print(C_AGN.z_min,'<z<',C_AGN.z_max)
    j0 = np.searchsorted(z_bins, C_AGN.z_min, side='right')-1
    j1 = np.searchsorted(z_bins, C_AGN.z_max, side='right')-1
    print(C_AGN.z_min, C_AGN.z_max)
    all_cat_outputs = []
    for jj in np.arange(j0, j1+1, 1):
        print('=='*50)
        z_bins_val = z_bins[jj]
        print(z_bins_val)
        p_2_catalogue_out = os.path.join( os.path.dirname(p_2_catalogue),
                        'AGN_list_z_'+str(np.round(z_bins_val,2))+'_sigma_'+C_AGN.str_scatter_0+'_fsat_'+C_AGN.str_fsat+'.fits' )
        t_sim_gal = t_sim_gal_full[(t_sim_gal_full['redshift_R']>=z_bins_val)&(t_sim_gal_full['redshift_R']<z_bins_val+DZ)]
        if len(t_sim_gal)==0:
            continue
        C_AGN.z_min = np.min(t_sim_gal['redshift_R'])
        C_AGN.z_max = np.max(t_sim_gal['redshift_R'])
        C_AGN.z_mean = np.mean(t_sim_gal['redshift_R'])
        dCx_min, dCx_max = np.min(t_sim_gal['x']), np.max(t_sim_gal['x'])
        dCy_min, dCy_max = np.min(t_sim_gal['y']), np.max(t_sim_gal['y'])
        dCz_min, dCz_max = np.min(t_sim_gal['z']), np.max(t_sim_gal['z'])
        #print(t_sim_gal.info())
        t_sim_gal.remove_columns([
                    'x',
                    'y',
                    'z',
                    #'Mvir',
                    'icl',
                    'id',
                    #'obs_sfr',
                    #'obs_sm',
                    'obs_uv',
                    'sfr',
                    'sm',
                    #'upid',
                    'A_UV',
                    'Mpeak',
                    'Vmax_Mpeak',
                    'desc_id',
                    'vmax',
                    'vx',
                    'vy',
                    'vz',
                    #'RA',
                    #'DEC',
                    'g_lat',
                    'g_lon',
                    'ecl_lat',
                    'ecl_lon',
                    'redshift_R',
                    #'redshift_S',
                    #'dL',
                    #'nH',
                    #'ebv'
                    ])
        #print(t_sim_gal.info())
        ##
        ## retrieves tabulated AGN FILE
        ##
        try:
            C_AGN.get_tabulated_AGN(z_bins_val, dCx_min, dCx_max, dCy_min, dCy_max, dCz_min, dCz_max)
            ##
            ## applies f_sat and duty cycle
            ##
            C_AGN.get_z_mass_id(t_sim_gal)
            ## scatter parameter is used here :
            C_AGN.abundance_matching()
            #
            C_AGN.get_obscured_fractions()
            C_AGN.compute_fluxes()
            C_AGN.AGN['agn_type'] = C_AGN.compute_agn_type(C_AGN.AGN['redshift_S'], C_AGN.AGN['LX_hard'], C_AGN.AGN['logNH'])
            C_AGN.compute_r_mag()

            t_out = C_AGN.AGN
            #print(t_out.info())
            t_out.remove_columns([
                        'dL',
                        'nH',
                        'ebv'
                        ])

            #print(t_out.info())

            t_out.write(p_2_catalogue_out, overwrite = True)
            print(p_2_catalogue_out, 'written', time.time()-t0)
            all_cat_outputs.append(p_2_catalogue_out)
        except(ValueError, ZeroDivisionError):
            print('Value Error, ZeroDivisionError')

    file_list = np.array(all_cat_outputs)
    t = Table.read(file_list[0])

    for p_2_file in file_list[1:]:
        print(p_2_file)
        t1 = Table.read(p_2_file)
        t = vstack(( t, t1))

    t_out = Table(t)
    t_out.write(path_2_output, overwrite=True)
    print(path_2_output, 'written')
    for p_2_file in file_list:
        os.remove(p_2_file)
