import os, sys
import numpy as np
import astropy.table as tbl
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as plt

#Set environmental variables
os.environ['UCHUU']='/home/idies/workspace/erosim/Uchuu'
os.environ['GIT_STMOD']='/home/idies/workspace/erosim/software/st_mod'
os.environ['GIT_STMOD_DATA']='/home/idies/workspace/erosim/software/st_mod_data'

#Read skymap data model
sky_map_hdu = tbl.Table.read(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS.fits') )

#Read list of tiles in eROSITA DE good footprint
good_tiles_list = np.loadtxt('/home/idies/workspace/erosim/center_in_foot.txt', dtype = str, unpack = True)

#Define seed list
clu_seed = list(range(19,118))
agn_seed = list(range(1,10))*11
t_tot_seed = list(zip(agn_seed, clu_seed))
tot_seed = [(str(seed[0]).zfill(3), str(seed[1]).zfill(3)) for seed in t_tot_seed]

#List of experiment names
GE_names = ['GE_e5_merge_AGNseed{0}_SimBKG_CLUseed{1}'.format(seed[0], seed[1]) for seed in tot_seed]+['GE_e4_merge_AGNseed{0}_SimBKG_CLUseed{1}'.format(seed[0], seed[1]) for seed in tot_seed]

SKYMAP = {}
for GE_name in GE_names:
    SKYMAP[GE_name] = tbl.Table.read(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS_'+GE_name+'.fits'))

#Set number of tiles to be done per batch
N_per_batch = int(sys.argv[1])

#Populate array
statusrows = []
to_do_global = []
already_done_global = []
for GE_name in GE_names:
    sky_map_hdu = SKYMAP[GE_name]

    #Add column to sky_map_hdu with tile id
    str_field_list = [str(sky_tile['SRVMAP']).zfill(6) for sky_tile in sky_map_hdu]
    sky_map_hdu.add_column(tbl.Column(str_field_list, name = 'tile_id'))

    #Select tiles in DE sky for simplicity
    skm_hdu_in_de_sky = sky_map_hdu[np.where((sky_map_hdu['OWNER']==2) | (sky_map_hdu['OWNER']==0))[0]]

    #Indexes of sky tiles in priority list
    skm_hdu_priority_idx = []
    skm_hdu_not_priority_idx = []
    for tid in skm_hdu_in_de_sky['tile_id']:
        if tid in list(good_tiles_list):    
            skm_hdu_priority_idx.append(list(skm_hdu_in_de_sky['tile_id']).index(tid))
        else:
            skm_hdu_not_priority_idx.append(list(skm_hdu_in_de_sky['tile_id']).index(tid))
    skm_hdu_priority_idx = np.array(skm_hdu_priority_idx)
    skm_hdu_not_priority_idx = np.array(skm_hdu_not_priority_idx)

    #Separate priority from not in priority
    skm_hdu_priority = skm_hdu_in_de_sky[skm_hdu_priority_idx]
    skm_hdu_not_priority = skm_hdu_in_de_sky[skm_hdu_not_priority_idx]

    #Identify those already done and those to be done
    to_process_priority = skm_hdu_priority[np.where((skm_hdu_priority['has_merged_events'])&(skm_hdu_priority['has_Sc1Cat']==False))[0]]
    already_done_priority = skm_hdu_priority[np.where((skm_hdu_priority['has_merged_events'])&(skm_hdu_priority['has_Sc1Cat']))[0]]
    to_process_not_priority = skm_hdu_not_priority[np.where((skm_hdu_not_priority['has_merged_events'])&(skm_hdu_not_priority['has_Sc1Cat']==False))[0]]
    already_done_not_priority = skm_hdu_not_priority[np.where((skm_hdu_not_priority['has_merged_events'])&(skm_hdu_not_priority['has_Sc1Cat']))[0]]

    #If there are no more priority tiles to be done, then set variable
    priority_done = len(to_process_priority) == 0

    #Stack back tables
    already_done_all = tbl.vstack([already_done_priority, already_done_not_priority], join_type = 'exact')
    to_process_all = tbl.vstack([to_process_priority, to_process_not_priority], join_type = 'exact')

    #Print information
    if len(already_done_all)+len(to_process_all) > 0:
        statusrows.append([len(to_process_all), len(already_done_all), len(already_done_all)/(len(already_done_all)+len(to_process_all)), (len(already_done_all)/(len(already_done_all)+len(to_process_all)))*100, priority_done, GE_name])
    else:
        statusrows.append([len(to_process_all), len(already_done_all), 0, 0, priority_done, GE_name])

    #Global numbers
    to_do_global.append(len(to_process_all))
    already_done_global.append(len(already_done_all))

    #If N_per_batch is not zero, then actually write the files
    if (N_per_batch > 0) & ((GE_name == 'GE_e4_merge_AGNseed005_SimBKG_CLUseed095') 
                          | (GE_name == 'GE_e4_merge_AGNseed006_SimBKG_CLUseed096') 
                          | (GE_name == 'GE_e4_merge_AGNseed007_SimBKG_CLUseed097')
                          | (GE_name == 'GE_e4_merge_AGNseed008_SimBKG_CLUseed098')
                          | (GE_name == 'GE_e4_merge_AGNseed009_SimBKG_CLUseed099')):
        
        #Do figure
        p2fig = os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'ra-dec-SKYMAPS_' + GE_name + '.png')
        plt.plot(skm_hdu_in_de_sky['RA_CEN'], skm_hdu_in_de_sky['DE_CEN'], 'k,', label='eRO DE')
        plt.plot(already_done_all['RA_CEN'], already_done_all['DE_CEN'], 'k+', label=str(len(already_done_all))+' done')
        plt.plot(to_process_all['RA_CEN'], to_process_all['DE_CEN'], 'rx', label=str(len(to_process_all))+' todo')
        plt.legend(loc=0)
        plt.title(GE_name)
        plt.savefig(p2fig)
        plt.clf()

        #If there are tiles left to process
        if len(to_process_all)>0:
            for kk in np.arange(0, len(to_process_all), N_per_batch):
                out_im1 = os.path.join('/home/idies/workspace/erosim', 'runs', GE_name + '_processing_'+str(kk).zfill(4)+'.sh')
                f_out = open(out_im1, 'w')
                f_out.write("""#!/bin/bash/ \n""")
                f_out.write("source activate heasoft \n")
                f_out.write("[ -r /home/idies/.healpix/3_50_Linux/config ] && . /home/idies/.healpix/3_50_Linux/config \n")
                f_out.write("source /opt/esass/bin/esass-init.sh \n")

                #Cycle tiles to process
                for sky_tile in to_process_all[kk: kk+N_per_batch]:
                    indir = os.path.join("/home/idies/workspace/erosim/Uchuu/LCerass/", sky_tile['tile_id'], GE_name)
                    esass_dir = os.path.join(indir, 'eSASS')
                    git_dir = os.path.join(os.environ['GIT_STMOD'], 'src/esass' )
                    path_2_event_file = os.path.join(indir, 'evt_'+sky_tile['tile_id']+'.fits')
                    if os.path.isfile(path_2_event_file) and os.path.isfile(os.path.join(esass_dir,sky_tile['tile_id']+"_pipeline_img1.sh")):
                        f_out.write ("cd "+esass_dir +" \n")
                        f_out.write ("sh "+sky_tile['tile_id']+"_pipeline_img1.sh"+" \n")
                        f_out.write ("sh "+sky_tile['tile_id']+"_pipeline_det1.sh"+" \n")
                        f_out.write ("sh "+sky_tile['tile_id']+"_pipeline_Src1.sh"+" \n")
                        f_out.write ("cd "+git_dir+" \n")
                        if not GE_name=='GE_e4_merge_SimBKG':
                            f_out.write ("python photon_matching_RS.py "+GE_name+" "+sky_tile['tile_id']+" \n")
                        f_out.write('# ====='+" \n")
                f_out.close()
                print(out_im1, 'written')

#Create table
statustab = tbl.Table(rows = statusrows, names = ['To_do', 'Done', 'Fraction', 'Percentage', 'Priority_done', 'Experiment_name'])
statustab['Fraction'].info.format = '.2f'
statustab['Percentage'].info.format = '.2f'
statustab.pprint_all()
print('='*100)
print('fields')
print(np.sum(to_do_global), 'todo')
print(np.sum(already_done_global), 'done')