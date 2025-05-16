# installation of sixte done on sciserver

# makes the image library for the galaxy cluster Uchuu run
export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/sixte_GAS/

export OMP_NUM_THREADS=8

# DONE no need to REDO
# python create_image_library_hotgas.py # DONE no need to REDO
# DONE no need to REDO
result stored here :
/home/idies/workspace/erosim/software/st_mod_data/data/models/model_GAS/profiles_010z015_1e14M2e14.fits


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
python merge_erosita_tile.py # DONE

# create links for the images and spectra in each folder
# python create_links_per_tile.py # DONE
# /home/idies/workspace/erosim/Uchuu/cluster_images
# /home/idies/workspace/erosim/Uchuu/cluster_Xspectra

python format_simput_xgas.py # DONE
# python format_simput_xgas_rewrite_images_radec.py # NOT NEEDED

export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/sixte_GAS/

export OMP_NUM_THREADS=8

nohup python extract_erosita_tile_XGAS.py 1  > logs/extract_erosita_tile_XGAS_1.log  & # DONE
nohup python extract_erosita_tile_XGAS.py 2  > logs/extract_erosita_tile_XGAS_2.log  & # DONE
nohup python extract_erosita_tile_XGAS.py 3  > logs/extract_erosita_tile_XGAS_3.log  & # DONE
nohup python extract_erosita_tile_XGAS.py 4  > logs/extract_erosita_tile_XGAS_4.log  & # DONE
nohup python extract_erosita_tile_XGAS.py 5  > logs/extract_erosita_tile_XGAS_5.log  & # DONE
nohup python extract_erosita_tile_XGAS.py 6  > logs/extract_erosita_tile_XGAS_6.log  & # DONE
nohup python extract_erosita_tile_XGAS.py 7  > logs/extract_erosita_tile_XGAS_7.log  & # ONGOING
nohup python extract_erosita_tile_XGAS.py 8  > logs/extract_erosita_tile_XGAS_8.log  & # ONGOING
nohup python extract_erosita_tile_XGAS.py 9  > logs/extract_erosita_tile_XGAS_9.log  & # ONGOING
nohup python extract_erosita_tile_XGAS.py 10 > logs/extract_erosita_tile_XGAS_10.log & # ONGOING
nohup python extract_erosita_tile_XGAS.py 11 > logs/extract_erosita_tile_XGAS_11.log & # ONGOING
nohup python extract_erosita_tile_XGAS.py 12 > logs/extract_erosita_tile_XGAS_12.log & # ONGOING
nohup python extract_erosita_tile_XGAS.py 13 > logs/extract_erosita_tile_XGAS_13.log & # ONGOING
nohup python extract_erosita_tile_XGAS.py 14 > logs/extract_erosita_tile_XGAS_14.log & # ONGOING
nohup python extract_erosita_tile_XGAS.py 15 > logs/extract_erosita_tile_XGAS_15.log & # ONGOING
nohup python extract_erosita_tile_XGAS.py 16 > logs/extract_erosita_tile_XGAS_16.log & # TODO
nohup python extract_erosita_tile_XGAS.py 17 > logs/extract_erosita_tile_XGAS_17.log & # TODO
nohup python extract_erosita_tile_XGAS.py 18 > logs/extract_erosita_tile_XGAS_18.log & # TODO
nohup python extract_erosita_tile_XGAS.py 19 > logs/extract_erosita_tile_XGAS_19.log & # TODO
nohup python extract_erosita_tile_XGAS.py 20 > logs/extract_erosita_tile_XGAS_20.log & # TODO
nohup python extract_erosita_tile_XGAS.py 21 > logs/extract_erosita_tile_XGAS_21.log & # TODO
nohup python extract_erosita_tile_XGAS.py 22 > logs/extract_erosita_tile_XGAS_22.log & # TODO
python merge_erosita_tile_XGAS.py # TODO

# makes the image library for the galaxy cluster Uchuu run
export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/sixte_GAS/

export OMP_NUM_THREADS=8

nohup python extract_erosita_tile_glist.py 1  > logs/extract_erosita_tile_1.log  & # DONE  'z0p00',
nohup python extract_erosita_tile_glist.py 2  > logs/extract_erosita_tile_2.log  & # DONE  'z0p02',
nohup python extract_erosita_tile_glist.py 3  > logs/extract_erosita_tile_3.log  & # DONE  'z0p05',
nohup python extract_erosita_tile_glist.py 4  > logs/extract_erosita_tile_4.log  & # DONE  'z0p09',
nohup python extract_erosita_tile_glist.py 5  > logs/extract_erosita_tile_5.log  & # DONE  'z0p14',
nohup python extract_erosita_tile_glist.py 6  > logs/extract_erosita_tile_6.log  & # DONE  'z0p19',
nohup python extract_erosita_tile_glist.py 7  > logs/extract_erosita_tile_7.log  & # DONE  'z0p25',
nohup python extract_erosita_tile_glist.py 8  > logs/extract_erosita_tile_8.log  & # DONE  'z0p30',
nohup python extract_erosita_tile_glist.py 9  > logs/extract_erosita_tile_9.log  & # DONE  'z0p36',
nohup python extract_erosita_tile_glist.py 10 > logs/extract_erosita_tile_10.log & # DONE  'z0p43',
nohup python extract_erosita_tile_glist.py 11 > logs/extract_erosita_tile_11.log & # DONE  'z0p49',
nohup python extract_erosita_tile_glist.py 12 > logs/extract_erosita_tile_12.log & # DONE  'z0p56',
nohup python extract_erosita_tile_glist.py 13 > logs/extract_erosita_tile_13.log & # DONE 'z0p63',
nohup python extract_erosita_tile_glist.py 14 > logs/extract_erosita_tile_14.log & # DONE 'z0p70',
nohup python extract_erosita_tile_glist.py 15 > logs/extract_erosita_tile_15.log & # DONE 'z0p78',
nohup python extract_erosita_tile_glist.py 16 > logs/extract_erosita_tile_16.log & # DONE 'z0p86',
nohup python extract_erosita_tile_glist.py 17 > logs/extract_erosita_tile_17.log & # DONE 'z0p94',
nohup python extract_erosita_tile_glist.py 18 > logs/extract_erosita_tile_18.log & # DONE 'z1p03',
nohup python extract_erosita_tile_glist.py 19 > logs/extract_erosita_tile_19.log & # DONE 'z1p12',
nohup python extract_erosita_tile_glist.py 20 > logs/extract_erosita_tile_20.log & # DONE 'z1p22',
nohup python extract_erosita_tile_glist.py 21 > logs/extract_erosita_tile_21.log & # DONE 'z1p32',
nohup python extract_erosita_tile_glist.py 22 > logs/extract_erosita_tile_22.log & # DONE 'z1p43',
nohup python extract_erosita_tile_glist.py 23 > logs/extract_erosita_tile_23.log & # DONE 'z1p54',
nohup python extract_erosita_tile_glist.py 24 > logs/extract_erosita_tile_24.log & # DONE 'z1p65',
nohup python extract_erosita_tile_glist.py 25 > logs/extract_erosita_tile_25.log & # DONE 'z1p77',
nohup python extract_erosita_tile_glist.py 26 > logs/extract_erosita_tile_26.log & # DONE 'z1p90',
nohup python extract_erosita_tile_glist.py 27 > logs/extract_erosita_tile_27.log & # DONE 'z2p03',
nohup python extract_erosita_tile_glist.py 28 > logs/extract_erosita_tile_28.log & # DONE 'z2p17',
nohup python extract_erosita_tile_glist.py 29 > logs/extract_erosita_tile_29.log & # DONE 'z2p31',
nohup python extract_erosita_tile_glist.py 30 > logs/extract_erosita_tile_30.log & # ONGOING 'z2p46',
nohup python extract_erosita_tile_glist.py 31 > logs/extract_erosita_tile_31.log & # ONGOING 'z2p62',
nohup python extract_erosita_tile_glist.py 32 > logs/extract_erosita_tile_32.log & # ONGOING 'z2p78',
nohup python extract_erosita_tile_glist.py 33 > logs/extract_erosita_tile_33.log & # ONGOING 'z2p95',
nohup python extract_erosita_tile_glist.py 34 > logs/extract_erosita_tile_34.log & # ONGOING 'z3p13',
nohup python extract_erosita_tile_glist.py 35 > logs/extract_erosita_tile_35.log & # ONGOING 'z3p32',
nohup python extract_erosita_tile_glist.py 36 > logs/extract_erosita_tile_36.log & # ONGOING 'z3p61',
# nohup python extract_erosita_tile_glist.py 37 > logs/extract_erosita_tile_37.log & # DONE 'z3p93',

python merge_erosita_tile_glist.py # TODO

export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/sixte_GAS/
export OMP_NUM_THREADS=8

# after the merging is finished (it deletes temporary files, careful !)
# python clean_erosita_tile.py # TODO AT THE END, or not at all.


# simulates events with seed fixed
python simulate_cluster_only_SEED_SKYMAP.py # DONE
# DONE for seeds 1,2,3,4,5,6,7,8,9


Now go to these folders :
sixte_AGN
sixte_fgbg
esass
