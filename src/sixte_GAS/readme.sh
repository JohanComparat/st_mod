# installation of sixte done on sciserver

# makes the image library for the galaxy cluster Uchuu run
export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/sixte_GAS/

# python create_image_library_hotgas.py # DONE no need to REDO

# re-cast cluster simput files into erosita tiles
# implement flux cut
# import numpy as np
# for kk in np.arange(24):
#     print("nohup python extract_erosita_tile.py "+str(kk)+" > logs/extract_erosita_tile_"+str(kk)+".log &")

# takes <1h
#
# # nohup python extract_erosita_tile.py 0 > logs/extract_erosita_tile_0.log & # DO NOT extract z=0. Skip it !
nohup python extract_erosita_tile.py 1  > logs/extract_erosita_tile_1.log  & # DONE
nohup python extract_erosita_tile.py 2  > logs/extract_erosita_tile_2.log  & # DONE
nohup python extract_erosita_tile.py 3  > logs/extract_erosita_tile_3.log  & # DONE
nohup python extract_erosita_tile.py 4  > logs/extract_erosita_tile_4.log  & # DONE
nohup python extract_erosita_tile.py 5  > logs/extract_erosita_tile_5.log  & # DONE
nohup python extract_erosita_tile.py 6  > logs/extract_erosita_tile_6.log  & # DONE
nohup python extract_erosita_tile.py 7  > logs/extract_erosita_tile_7.log  & # DONE
nohup python extract_erosita_tile.py 8  > logs/extract_erosita_tile_8.log  & # DONE
nohup python extract_erosita_tile.py 9  > logs/extract_erosita_tile_9.log  & # DONE
nohup python extract_erosita_tile.py 10 > logs/extract_erosita_tile_10.log & # DONE
nohup python extract_erosita_tile.py 11 > logs/extract_erosita_tile_11.log & # DONE
nohup python extract_erosita_tile.py 12 > logs/extract_erosita_tile_12.log & # DONE
nohup python extract_erosita_tile.py 13 > logs/extract_erosita_tile_13.log & # DONE
nohup python extract_erosita_tile.py 14 > logs/extract_erosita_tile_14.log & # DONE
nohup python extract_erosita_tile.py 15 > logs/extract_erosita_tile_15.log & # DONE
nohup python extract_erosita_tile.py 16 > logs/extract_erosita_tile_16.log & # DONE
nohup python extract_erosita_tile.py 17 > logs/extract_erosita_tile_17.log & # DONE
nohup python extract_erosita_tile.py 18 > logs/extract_erosita_tile_18.log & # DONE
nohup python extract_erosita_tile.py 19 > logs/extract_erosita_tile_19.log & # DONE
nohup python extract_erosita_tile.py 20 > logs/extract_erosita_tile_20.log & # DONE
nohup python extract_erosita_tile.py 21 > logs/extract_erosita_tile_21.log & # DONE
nohup python extract_erosita_tile.py 22 > logs/extract_erosita_tile_22.log & # DONE

# after the scripts above are finished
python merge_erosita_tile.py # DONE
# after the merging is finished (it deletes temporary files !)
python clean_erosita_tile.py # DONE

# create links for the images and spectra in each folder
python create_links_per_tile.py # DONE
/home/idies/workspace/erosim/Uchuu/cluster_images
/home/idies/workspace/erosim/Uchuu/cluster_Xspectra

python format_simput_xgas.py # DONE


# retrieve real events from complete repo
# ~/workspace/Storage/comparat/persistent/data/data_s4_c030 # rsync finished
# ~/workspace/Storage/comparat/persistent/data/data_s5_c030 # rsync finished
# rsync DONE

# copy c030 s4 and s5 into erosim
python copy_events.py # ONGOING




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
