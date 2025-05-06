# installation of sixte done on sciserver

# makes the image library for the galaxy cluster Uchuu run
export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/sixte_GAS/

export OMP_NUM_THREADS=8

# python create_image_library_hotgas.py # DONE no need to REDO

# re-cast cluster simput files into erosita tiles
# implement flux cut
# import numpy as np
# for kk in np.arange(24):
#     print("nohup python extract_erosita_tile.py "+str(kk)+" > logs/extract_erosita_tile_"+str(kk)+".log &")

# takes <1h
#
# # nohup python extract_erosita_tile.py 0 > logs/extract_erosita_tile_0.log & # DO NOT extract z=0. Skip it !
# nohup python extract_erosita_tile.py 1  > logs/extract_erosita_tile_1.log  & # DONE
# nohup python extract_erosita_tile.py 2  > logs/extract_erosita_tile_2.log  & # DONE
# nohup python extract_erosita_tile.py 3  > logs/extract_erosita_tile_3.log  & # DONE
# nohup python extract_erosita_tile.py 4  > logs/extract_erosita_tile_4.log  & # DONE
# nohup python extract_erosita_tile.py 5  > logs/extract_erosita_tile_5.log  & # DONE
# nohup python extract_erosita_tile.py 6  > logs/extract_erosita_tile_6.log  & # DONE
# nohup python extract_erosita_tile.py 7  > logs/extract_erosita_tile_7.log  & # DONE
# nohup python extract_erosita_tile.py 8  > logs/extract_erosita_tile_8.log  & # DONE
# nohup python extract_erosita_tile.py 9  > logs/extract_erosita_tile_9.log  & # DONE
# nohup python extract_erosita_tile.py 10 > logs/extract_erosita_tile_10.log & # DONE
# nohup python extract_erosita_tile.py 11 > logs/extract_erosita_tile_11.log & # DONE
# nohup python extract_erosita_tile.py 12 > logs/extract_erosita_tile_12.log & # DONE
# nohup python extract_erosita_tile.py 13 > logs/extract_erosita_tile_13.log & # DONE
# nohup python extract_erosita_tile.py 14 > logs/extract_erosita_tile_14.log & # DONE
# nohup python extract_erosita_tile.py 15 > logs/extract_erosita_tile_15.log & # DONE
# nohup python extract_erosita_tile.py 16 > logs/extract_erosita_tile_16.log & # DONE
# nohup python extract_erosita_tile.py 17 > logs/extract_erosita_tile_17.log & # DONE
# nohup python extract_erosita_tile.py 18 > logs/extract_erosita_tile_18.log & # DONE
# nohup python extract_erosita_tile.py 19 > logs/extract_erosita_tile_19.log & # DONE
# nohup python extract_erosita_tile.py 20 > logs/extract_erosita_tile_20.log & # DONE
# nohup python extract_erosita_tile.py 21 > logs/extract_erosita_tile_21.log & # DONE
# nohup python extract_erosita_tile.py 22 > logs/extract_erosita_tile_22.log & # DONE
# installation of sixte done on sciserver

# makes the image library for the galaxy cluster Uchuu run
export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/sixte_GAS/

export OMP_NUM_THREADS=8

# nohup python extract_erosita_tile_glist.py 1  > logs/extract_erosita_tile_1.log  & # DONE  'z0p00',
# nohup python extract_erosita_tile_glist.py 2  > logs/extract_erosita_tile_2.log  & # DONE  'z0p02',
# nohup python extract_erosita_tile_glist.py 3  > logs/extract_erosita_tile_3.log  & # DONE  'z0p05',
# nohup python extract_erosita_tile_glist.py 4  > logs/extract_erosita_tile_4.log  & # DONE  'z0p09',
# nohup python extract_erosita_tile_glist.py 5  > logs/extract_erosita_tile_5.log  & # DONE  'z0p14',
# nohup python extract_erosita_tile_glist.py 6  > logs/extract_erosita_tile_6.log  & # DONE  'z0p19',
# nohup python extract_erosita_tile_glist.py 7  > logs/extract_erosita_tile_7.log  & # DONE  'z0p25',
# nohup python extract_erosita_tile_glist.py 8  > logs/extract_erosita_tile_8.log  & # DONE  'z0p30',
# nohup python extract_erosita_tile_glist.py 9  > logs/extract_erosita_tile_9.log  & # DONE  'z0p36',
# nohup python extract_erosita_tile_glist.py 10 > logs/extract_erosita_tile_10.log & # DONE  'z0p43',
# nohup python extract_erosita_tile_glist.py 11 > logs/extract_erosita_tile_11.log & # DONE  'z0p49',
# nohup python extract_erosita_tile_glist.py 12 > logs/extract_erosita_tile_12.log & # DONE  'z0p56',
# nohup python extract_erosita_tile_glist.py 13 > logs/extract_erosita_tile_13.log & # DONE 'z0p63',
# nohup python extract_erosita_tile_glist.py 14 > logs/extract_erosita_tile_14.log & # DONE 'z0p70',
# nohup python extract_erosita_tile_glist.py 15 > logs/extract_erosita_tile_15.log & # DONE 'z0p78',
# nohup python extract_erosita_tile_glist.py 16 > logs/extract_erosita_tile_16.log & # DONE 'z0p86',
# nohup python extract_erosita_tile_glist.py 17 > logs/extract_erosita_tile_17.log & # DONE 'z0p94',
# nohup python extract_erosita_tile_glist.py 18 > logs/extract_erosita_tile_18.log & # DONE 'z1p03',
# nohup python extract_erosita_tile_glist.py 19 > logs/extract_erosita_tile_19.log & # DONE 'z1p12',
# nohup python extract_erosita_tile_glist.py 20 > logs/extract_erosita_tile_20.log & # DONE 'z1p22',
# nohup python extract_erosita_tile_glist.py 21 > logs/extract_erosita_tile_21.log & # DONE 'z1p32',
nohup python extract_erosita_tile_glist.py 22 > logs/extract_erosita_tile_22.log & # ONGOING 'z1p43',
nohup python extract_erosita_tile_glist.py 23 > logs/extract_erosita_tile_23.log & # ONGOING 'z1p54',
nohup python extract_erosita_tile_glist.py 24 > logs/extract_erosita_tile_24.log & # ONGOING 'z1p65',
nohup python extract_erosita_tile_glist.py 25 > logs/extract_erosita_tile_25.log & # TODO 'z1p77',
nohup python extract_erosita_tile_glist.py 26 > logs/extract_erosita_tile_26.log & # TODO 'z1p90',
nohup python extract_erosita_tile_glist.py 27 > logs/extract_erosita_tile_27.log & # TODO 'z2p03',
nohup python extract_erosita_tile_glist.py 28 > logs/extract_erosita_tile_28.log & # TODO 'z2p17',
nohup python extract_erosita_tile_glist.py 29 > logs/extract_erosita_tile_29.log & # TODO 'z2p31',
nohup python extract_erosita_tile_glist.py 30 > logs/extract_erosita_tile_30.log & # TODO 'z2p46',
nohup python extract_erosita_tile_glist.py 31 > logs/extract_erosita_tile_31.log & # TODO 'z2p62',
nohup python extract_erosita_tile_glist.py 32 > logs/extract_erosita_tile_32.log & # TODO 'z2p78',
nohup python extract_erosita_tile_glist.py 33 > logs/extract_erosita_tile_33.log & # TODO 'z2p95',
nohup python extract_erosita_tile_glist.py 34 > logs/extract_erosita_tile_34.log & # TODO 'z3p13',
nohup python extract_erosita_tile_glist.py 35 > logs/extract_erosita_tile_35.log & # TODO 'z3p32',
nohup python extract_erosita_tile_glist.py 36 > logs/extract_erosita_tile_36.log & # TODO 'z3p61',
nohup python extract_erosita_tile_glist.py 37 > logs/extract_erosita_tile_37.log & # TODO 'z3p93',

TODO : extract files with the cluster model computed

# after the scripts above are finished
python merge_erosita_tile.py # DONE ERO_DE DONE ERO_RU
python merge_erosita_tile_glist.py # TODO
TODO : merge files with the cluster model computed

# after the merging is finished (it deletes temporary files !)
# python clean_erosita_tile.py # TODO AT THE END, or not

# create links for the images and spectra in each folder
# python create_links_per_tile.py # DONE
# /home/idies/workspace/erosim/Uchuu/cluster_images
# /home/idies/workspace/erosim/Uchuu/cluster_Xspectra

python format_simput_xgas.py # DONE ERO_DE DONE  ERO_RU
# python format_simput_xgas_rewrite_images_radec.py # NOT NEEDED

# retrieve real events from complete repo
# ~/workspace/Storage/comparat/persistent/data/data_s4_c030 # rsync finished
# ~/workspace/Storage/comparat/persistent/data/data_s5_c030 # rsync finished
# rsync DONE
# copy c030 s4 and s5 into erosim
# python copy_events.py # DONE

# simulates events with seed fixed to 001
# only does the first tile : 121048
python simulate_cluster_only_SEED_SKYMAP.py # # DONE ERO_DE DONE ERO_RU v2.7 sixte runs
# SEED 1 done
# SEED 2 ONGOING
# SEED 3 ONGOING


ls /home/idies/workspace/erosim/Uchuu/LCerass/??????/Att_eRASS8_sixte_v27_SEED_001_events_cluster/t0erass_ccd1_evt.fits > list_clu_ccd1.list
ls /home/idies/workspace/erosim/Uchuu/LCerass/??????/Att_eRASS8_sixte_v27_SEED_001_events_cluster/t0erass_ccd1_raw.fits > list_clu_ccd1_raw.list
wc -l list*

python simulate_eRASS45_cluster_only_SEED_SKYMAP.py # v2.7 sixte runs TODO

python plot_events_with_catalog.py 164087 "eRASS8_SEED_001_events_cluster_2025_04"

/opt/sixte/share/sixte/instruments/srg/erosita/erosita_1.xml
/opt/sixte/share/sixte/instruments/srg/erosita/sixte_erormf_normalized_singles_20170725.rmf
/opt/sixte/share/sixte/instruments/srg/erosita/erosita_pha2pi_v20210719.fits

# simulate data with sixte and events using attitude fiels in the events !!!

for eRASS:5,
let's use this field as for testing:
81141
data_s5_c030/081141/c030/sm05_081141_020_FlareGTI_c030.fits

for eRASS:4,
let's use this field as for testing:
164087
data_s4_c030/164087/c030/sm04_164087_020_FlareGTI_c030.fits


# what to simulate
# take attitude file from the true events
take simput and quantities from $UCHUU DIR

# simulation setup
Make a long eRASS:8 sixte simulation

Make a long eRASS:4 or 5 sixte simulation

export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/sixte_GAS/


for jj in n.arange(100):
    print("nohup python simulate_cluster_only_SEED_SKYMAP.py "+str(jj).zfill(3)+" > sim_cluster_seed"+str(jj).zfill(3)+".log &")


cd $GIT_ERASS_SIM/sixte
nohup python simulate_cluster_only_SEED_SKYMAP.py 001 > sim_cluster_seed001.log &
nohup python simulate_cluster_only_SEED_SKYMAP.py 002 > sim_cluster_seed002.log &
nohup python simulate_cluster_only_SEED_SKYMAP.py 003 > sim_cluster_seed003.log &
