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
#       6ujYKT&D1Kw7b&$1V$gvv!cPM
# takes <1h
#
# # nohup python extract_erosita_tile.py 0 > logs/extract_erosita_tile_0.log & # DO NOT extract z=0. Skip it !
# nohup python extract_erosita_tile.py 1  > logs/extract_erosita_tile_1.log  & # DONE   'z0p00',
# nohup python extract_erosita_tile.py 2  > logs/extract_erosita_tile_2.log  & # DONE   'z0p02',
# nohup python extract_erosita_tile.py 3  > logs/extract_erosita_tile_3.log  & # DONE   'z0p05',
# nohup python extract_erosita_tile.py 4  > logs/extract_erosita_tile_4.log  & # DONE   'z0p09',
# nohup python extract_erosita_tile.py 5  > logs/extract_erosita_tile_5.log  & # DONE   'z0p14',
# nohup python extract_erosita_tile.py 6  > logs/extract_erosita_tile_6.log  & # DONE   'z0p19',
# nohup python extract_erosita_tile.py 7  > logs/extract_erosita_tile_7.log  & # DONE   'z0p25',
# nohup python extract_erosita_tile.py 8  > logs/extract_erosita_tile_8.log  & # DONE   'z0p30',
# nohup python extract_erosita_tile.py 9  > logs/extract_erosita_tile_9.log  & # DONE   'z0p36',
# nohup python extract_erosita_tile.py 10 > logs/extract_erosita_tile_10.log & # DONE   'z0p43',
# nohup python extract_erosita_tile.py 11 > logs/extract_erosita_tile_11.log & # DONE   'z0p49',
# nohup python extract_erosita_tile.py 12 > logs/extract_erosita_tile_12.log & # DONE   'z0p56',
# nohup python extract_erosita_tile.py 13 > logs/extract_erosita_tile_13.log & # DONE   'z0p63',
# nohup python extract_erosita_tile.py 14 > logs/extract_erosita_tile_14.log & # DONE   'z0p70',
# nohup python extract_erosita_tile.py 15 > logs/extract_erosita_tile_15.log & # DONE   'z0p78',
# nohup python extract_erosita_tile.py 16 > logs/extract_erosita_tile_16.log & # DONE   'z0p86',
# nohup python extract_erosita_tile.py 17 > logs/extract_erosita_tile_17.log & # DONE   'z0p94',
# nohup python extract_erosita_tile.py 18 > logs/extract_erosita_tile_18.log & # DONE   'z1p03',
# nohup python extract_erosita_tile.py 19 > logs/extract_erosita_tile_19.log & # DONE   'z1p12',
# nohup python extract_erosita_tile.py 20 > logs/extract_erosita_tile_20.log & # DONE   'z1p22',
# nohup python extract_erosita_tile.py 21 > logs/extract_erosita_tile_21.log & # DONE   'z1p32',
# nohup python extract_erosita_tile.py 22 > logs/extract_erosita_tile_22.log & # DONE   'z1p43',
# nohup python extract_erosita_tile.py 23 > logs/extract_erosita_tile_23.log & # DONE   'z1p54',
# nohup python extract_erosita_tile.py 24 > logs/extract_erosita_tile_24.log & # DONE   'z1p65',
# nohup python extract_erosita_tile.py 25 > logs/extract_erosita_tile_25.log & # DONE   'z1p77',
# nohup python extract_erosita_tile.py 26 > logs/extract_erosita_tile_26.log & # DONE   'z1p90',
# nohup python extract_erosita_tile.py 27 > logs/extract_erosita_tile_27.log & # DONE   'z2p03',
# nohup python extract_erosita_tile.py 28 > logs/extract_erosita_tile_28.log & # DONE   'z2p17',
# nohup python extract_erosita_tile.py 29 > logs/extract_erosita_tile_29.log & # DONE   'z2p31',
# nohup python extract_erosita_tile.py 30 > logs/extract_erosita_tile_30.log & # DONE   'z2p46',
# nohup python extract_erosita_tile.py 31 > logs/extract_erosita_tile_31.log & # DONE   'z2p62',
# nohup python extract_erosita_tile.py 32 > logs/extract_erosita_tile_32.log & # DONE   'z2p78',
# nohup python extract_erosita_tile.py 33 > logs/extract_erosita_tile_33.log & # DONE   'z2p95',
# nohup python extract_erosita_tile.py 34 > logs/extract_erosita_tile_34.log & # DONE   'z3p13',
# nohup python extract_erosita_tile.py 35 > logs/extract_erosita_tile_35.log & # DONE   'z3p32',
# nohup python extract_erosita_tile.py 36 > logs/extract_erosita_tile_36.log & # DONE   'z3p61',
# nohup python extract_erosita_tile.py 37 > logs/extract_erosita_tile_37.log & # DONE   'z3p93',


# TODO
# TODO

# after the scripts above are finished
nohup python merge_erosita_tile.py > logs/merge_erosita_tile.log & # DONE
# after the merging is finished (it deletes temporary files !)
# python clean_erosita_tile.py # not necessary, only at the end

# create links for the images and spectra in each folder
# python create_links_per_tile.py # DONE
# /home/idies/workspace/erosim/simput/AGNspectra_V2

python format_simput_AGN.py # DONE



# simulates events with seed fixed to 001
# only does the first tile : 121048
python simulate_AGN_only_SEED_SKYMAP.py # ONGOING v2.7 sixte runs
# python simulate_eRASS45_cluster_only_SEED_SKYMAP.py # ONGOING v2.7 sixte runs

ls /home/idies/workspace/erosim/Uchuu/LCerass/??????/eRASS8_SEED_001_events_AGN_2025_04/t0erass_ccd1_evt.fits > list_agn_ccd1.list
ls /home/idies/workspace/erosim/Uchuu/LCerass/??????/eRASS8_SEED_001_events_AGN_2025_04/t0erass_ccd1_raw.fits > list_agn_ccd1_raw.list
wc -l list*


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
