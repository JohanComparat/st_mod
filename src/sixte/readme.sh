# installation of sixte done on sciserver

# makes the image library for the galaxy cluster Uchuu run
export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/sixte
python create_image_library_hotgas.py # DONE

# re-cast cluster simput files into erosita tiles
# implement flux cut
# import numpy as np
# for kk in np.arange(24):
#     print("nohup python extract_erosita_tile.py "+str(kk)+" > logs/extract_erosita_tile_"+str(kk)+".log &")

# takes <1h
#
nohup python extract_erosita_tile.py 0 > logs/extract_erosita_tile_0.log & # TODO
nohup python extract_erosita_tile.py 1 > logs/extract_erosita_tile_1.log & # TODO
nohup python extract_erosita_tile.py 2 > logs/extract_erosita_tile_2.log & # TODO
nohup python extract_erosita_tile.py 3 > logs/extract_erosita_tile_3.log &   # TODO
nohup python extract_erosita_tile.py 4 > logs/extract_erosita_tile_4.log &   # TODO
nohup python extract_erosita_tile.py 5 > logs/extract_erosita_tile_5.log &   # TODO
nohup python extract_erosita_tile.py 6 > logs/extract_erosita_tile_6.log &   # TODO
nohup python extract_erosita_tile.py 7 > logs/extract_erosita_tile_7.log &   # TODO
nohup python extract_erosita_tile.py 8 > logs/extract_erosita_tile_8.log &   # TODO
nohup python extract_erosita_tile.py 9 > logs/extract_erosita_tile_9.log &   # TODO
nohup python extract_erosita_tile.py 10 > logs/extract_erosita_tile_10.log & # TODO
nohup python extract_erosita_tile.py 11 > logs/extract_erosita_tile_11.log & # TODO
nohup python extract_erosita_tile.py 12 > logs/extract_erosita_tile_12.log & # TODO
nohup python extract_erosita_tile.py 13 > logs/extract_erosita_tile_13.log & # TODO
nohup python extract_erosita_tile.py 14 > logs/extract_erosita_tile_14.log & # TODO
nohup python extract_erosita_tile.py 15 > logs/extract_erosita_tile_15.log & # TODO
nohup python extract_erosita_tile.py 16 > logs/extract_erosita_tile_16.log & # TODO
nohup python extract_erosita_tile.py 17 > logs/extract_erosita_tile_17.log & # TODO
nohup python extract_erosita_tile.py 18 > logs/extract_erosita_tile_18.log & # TODO
nohup python extract_erosita_tile.py 19 > logs/extract_erosita_tile_19.log & # TODO
nohup python extract_erosita_tile.py 20 > logs/extract_erosita_tile_20.log & # TODO
nohup python extract_erosita_tile.py 21 > logs/extract_erosita_tile_21.log & # TODO
nohup python extract_erosita_tile.py 22 > logs/extract_erosita_tile_22.log & # TODO

# after the scripts above are finished
python merge_erosita_tile.py # TODO

# after the merging is finished (it deletes temporary files !)
python clean_erosita_tile.py # TODO

# create links for the images and spectra in each folder
python create_links_per_tile.py # DONE

/home/idies/workspace/erosim/Uchuu/cluster_images
/home/idies/workspace/erosim/Uchuu/cluster_Xspectra


# retrieve real events
# rsync ONGOING
~/workspace/Storage/comparat/persistent/data/data_s4_c030
~/workspace/Storage/comparat/persistent/data/data_s5_c030

szr16jdjsgpd:q5=v0jJ

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



choose a field ID on the erosita tiles and its metadata

# hot gas
haloes that comply to
z_max(M500c) = z at which M500c = 20'
will have an image for the simulation

haloes that do not comply to this threshold are simulated as PSF


In a given redshift bin of width 0.01
create a set of images (on sky pixels) that represent the CC/NCC and the extent and possible R500c.
Say take 20 profile images in each R500c bin (log10 bins of 0.2 or 0.1)
Depending on the total number on has to write.

Use the mean ellipticity in each R500c bin ? (is iti available for Uchuu ?)

Single temperature apec model


# AGN

C019 model extended to lower stellar masses
Add dependence on SFR ?


# XRB

Lehmer model ?


# write all the sixte simputs directly

/data56s/comparat/erosim/data_s5/??????/Uchuu/simput_gas
/data56s/comparat/erosim/data_s5/??????/Uchuu/simput_gas/images
/data56s/comparat/erosim/data_s5/??????/Uchuu/simput_agn
/data56s/comparat/erosim/data_s5/??????/Uchuu/simput_xrb
/data56s/comparat/erosim/data_s5/??????/Uchuu/simput_bkg/images


# simulation setup

Make a long eRASS:8 sixte simulation



