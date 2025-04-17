# installation of sixte done on sciserver

# makes the image library for the galaxy cluster Uchuu run
export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/sixte_AGN/

# re-cast cluster simput files into erosita tiles
# implement flux cut
# import numpy as np
# for kk in np.arange(237):
#     print("nohup python extract_erosita_tile.py "+str(kk)+" > logs/extract_erosita_tile_"+str(kk)+".log &")

# takes <1h
#
# # nohup python extract_erosita_tile.py 0 > logs/extract_erosita_tile_0.log & # DO NOT extract z=0. Skip it !
nohup python extract_erosita_tile.py 1  > logs/extract_erosita_tile_1.log  & # TODO   'z0p00',
nohup python extract_erosita_tile.py 2  > logs/extract_erosita_tile_2.log  & # TODO   'z0p02',
nohup python extract_erosita_tile.py 3  > logs/extract_erosita_tile_3.log  & # TODO   'z0p05',
nohup python extract_erosita_tile.py 4  > logs/extract_erosita_tile_4.log  & # TODO   'z0p09',
nohup python extract_erosita_tile.py 5  > logs/extract_erosita_tile_5.log  & # TODO   'z0p14',
nohup python extract_erosita_tile.py 6  > logs/extract_erosita_tile_6.log  & # TODO   'z0p19',
nohup python extract_erosita_tile.py 7  > logs/extract_erosita_tile_7.log  & # TODO   'z0p25',
nohup python extract_erosita_tile.py 8  > logs/extract_erosita_tile_8.log  & # TODO   'z0p30',
nohup python extract_erosita_tile.py 9  > logs/extract_erosita_tile_9.log  & # TODO   'z0p36',
nohup python extract_erosita_tile.py 10 > logs/extract_erosita_tile_10.log & # TODO   'z0p43',
nohup python extract_erosita_tile.py 11 > logs/extract_erosita_tile_11.log & # TODO   'z0p49',
nohup python extract_erosita_tile.py 12 > logs/extract_erosita_tile_12.log & # TODO   'z0p56',
nohup python extract_erosita_tile.py 13 > logs/extract_erosita_tile_13.log & # TODO   'z0p63',
nohup python extract_erosita_tile.py 14 > logs/extract_erosita_tile_14.log & # TODO   'z0p70',
nohup python extract_erosita_tile.py 15 > logs/extract_erosita_tile_15.log & # TODO   'z0p78',
nohup python extract_erosita_tile.py 16 > logs/extract_erosita_tile_16.log & # TODO   'z0p86',
nohup python extract_erosita_tile.py 17 > logs/extract_erosita_tile_17.log & # TODO   'z0p94',
nohup python extract_erosita_tile.py 18 > logs/extract_erosita_tile_18.log & # TODO   'z1p03',
nohup python extract_erosita_tile.py 19 > logs/extract_erosita_tile_19.log & # TODO   'z1p12',
nohup python extract_erosita_tile.py 20 > logs/extract_erosita_tile_20.log & # TODO   'z1p22',
nohup python extract_erosita_tile.py 21 > logs/extract_erosita_tile_21.log & # TODO   'z1p32',
nohup python extract_erosita_tile.py 22 > logs/extract_erosita_tile_22.log & # TODO   'z1p43',
nohup python extract_erosita_tile.py 23 > logs/extract_erosita_tile_23.log & # TODO   'z1p54',
nohup python extract_erosita_tile.py 24 > logs/extract_erosita_tile_24.log & # TODO   'z1p65',
nohup python extract_erosita_tile.py 25 > logs/extract_erosita_tile_25.log & # TODO   'z1p77',
nohup python extract_erosita_tile.py 26 > logs/extract_erosita_tile_26.log & # TODO   'z1p90',
nohup python extract_erosita_tile.py 27 > logs/extract_erosita_tile_27.log & # TODO   'z2p03',
nohup python extract_erosita_tile.py 28 > logs/extract_erosita_tile_28.log & # TODO   'z2p17',
nohup python extract_erosita_tile.py 29 > logs/extract_erosita_tile_29.log & # TODO   'z2p31',
nohup python extract_erosita_tile.py 30 > logs/extract_erosita_tile_30.log & # TODO   'z2p46',
nohup python extract_erosita_tile.py 31 > logs/extract_erosita_tile_31.log & # TODO   'z2p62',
nohup python extract_erosita_tile.py 32 > logs/extract_erosita_tile_32.log & # TODO   'z2p78',
nohup python extract_erosita_tile.py 33 > logs/extract_erosita_tile_33.log & # TODO   'z2p95',
nohup python extract_erosita_tile.py 34 > logs/extract_erosita_tile_34.log & # TODO   'z3p13',
nohup python extract_erosita_tile.py 35 > logs/extract_erosita_tile_35.log & # TODO   'z3p32',
# nohup python extract_erosita_tile.py 36 > logs/extract_erosita_tile_36.log & # NO DATA 'z3p61',
# nohup python extract_erosita_tile.py 37 > logs/extract_erosita_tile_37.log & # NO DATA 'z3p93',
# nohup python extract_erosita_tile.py 38 > logs/extract_erosita_tile_38.log & # NO DATA 'z4p27',
# nohup python extract_erosita_tile.py 39 > logs/extract_erosita_tile_39.log & # NO DATA 'z4p63',
# nohup python extract_erosita_tile.py 40 > logs/extract_erosita_tile_10.log & # NO DATA 'z5p15',
# nohup python extract_erosita_tile.py 41 > logs/extract_erosita_tile_11.log & # NO DATA 'z5p73'
# nohup python extract_erosita_tile.py 42 > logs/extract_erosita_tile_12.log & # NO DATA
# nohup python extract_erosita_tile.py 43 > logs/extract_erosita_tile_13.log & # NO DATA
# nohup python extract_erosita_tile.py 44 > logs/extract_erosita_tile_14.log & # NO DATA
# nohup python extract_erosita_tile.py 45 > logs/extract_erosita_tile_15.log & # NO DATA
# nohup python extract_erosita_tile.py 46 > logs/extract_erosita_tile_16.log & # NO DATA
# nohup python extract_erosita_tile.py 47 > logs/extract_erosita_tile_17.log & # NO DATA
# nohup python extract_erosita_tile.py 48 > logs/extract_erosita_tile_18.log & # NO DATA
# nohup python extract_erosita_tile.py 49 > logs/extract_erosita_tile_19.log & # NO DATA

# TODO
# TODO

# after the scripts above are finished
python merge_erosita_tile.py # DONE
# after the merging is finished (it deletes temporary files !)
python clean_erosita_tile.py # TODO AFTER

# create links for the images and spectra in each folder
# python create_links_per_tile.py # DONE
# /home/idies/workspace/erosim/Uchuu/cluster_images
# /home/idies/workspace/erosim/Uchuu/cluster_Xspectra

python format_simput_xgas.py # DONE


# retrieve real events from complete repo
# ~/workspace/Storage/comparat/persistent/data/data_s4_c030 # rsync finished
# ~/workspace/Storage/comparat/persistent/data/data_s5_c030 # rsync finished
# rsync DONE
# copy c030 s4 and s5 into erosim
# python copy_events.py # DONE

# simulates events with seed fixed to 001
# only does the first tile : 121048
python simulate_cluster_only_SEED_SKYMAP.py # ONGOING v2.7 sixte runs
python simulate_eRASS45_cluster_only_SEED_SKYMAP.py # ONGOING v2.7 sixte runs





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
cd $GIT_STMOD/src/sixte_AGN/


for jj in n.arange(100):
    print("nohup python simulate_cluster_only_SEED_SKYMAP.py "+str(jj).zfill(3)+" > sim_cluster_seed"+str(jj).zfill(3)+".log &")


cd $GIT_ERASS_SIM/sixte
nohup python simulate_cluster_only_SEED_SKYMAP.py 001 > sim_cluster_seed001.log &
nohup python simulate_cluster_only_SEED_SKYMAP.py 002 > sim_cluster_seed002.log &
nohup python simulate_cluster_only_SEED_SKYMAP.py 003 > sim_cluster_seed003.log &
