# installation of sixte done on sciserver

# makes the image library for the galaxy cluster Uchuu run
export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/sixte
python create_image_library_hotgas.py # ONGOING sciserver

# images and spectra are here :
/home/idies/workspace/erosim/Uchuu/cluster_images
/home/idies/workspace/erosim/Uchuu/cluster_Xspectra

# re-cast cluster simput files into erosita tiles
# implement flux cut
import numpy as np
for kk in np.arange(24):
    print("nohup python extract_erosita_tile.py "+str(kk)+" > logs/extract_erosita_tile_"+str(kk)+".log &")

nohup python extract_erosita_tile.py 0 > logs/extract_erosita_tile_0.log &
nohup python extract_erosita_tile.py 1 > logs/extract_erosita_tile_1.log &
nohup python extract_erosita_tile.py 2 > logs/extract_erosita_tile_2.log &
nohup python extract_erosita_tile.py 3 > logs/extract_erosita_tile_3.log &
nohup python extract_erosita_tile.py 4 > logs/extract_erosita_tile_4.log &
nohup python extract_erosita_tile.py 5 > logs/extract_erosita_tile_5.log &
nohup python extract_erosita_tile.py 6 > logs/extract_erosita_tile_6.log &
nohup python extract_erosita_tile.py 7 > logs/extract_erosita_tile_7.log &
nohup python extract_erosita_tile.py 8 > logs/extract_erosita_tile_8.log &
nohup python extract_erosita_tile.py 9 > logs/extract_erosita_tile_9.log &
nohup python extract_erosita_tile.py 10 > logs/extract_erosita_tile_10.log &
nohup python extract_erosita_tile.py 11 > logs/extract_erosita_tile_11.log &
nohup python extract_erosita_tile.py 12 > logs/extract_erosita_tile_12.log &
nohup python extract_erosita_tile.py 13 > logs/extract_erosita_tile_13.log &
nohup python extract_erosita_tile.py 14 > logs/extract_erosita_tile_14.log &
nohup python extract_erosita_tile.py 15 > logs/extract_erosita_tile_15.log &
nohup python extract_erosita_tile.py 16 > logs/extract_erosita_tile_16.log &
nohup python extract_erosita_tile.py 17 > logs/extract_erosita_tile_17.log &
nohup python extract_erosita_tile.py 18 > logs/extract_erosita_tile_18.log &
nohup python extract_erosita_tile.py 19 > logs/extract_erosita_tile_19.log &
nohup python extract_erosita_tile.py 20 > logs/extract_erosita_tile_20.log &
nohup python extract_erosita_tile.py 21 > logs/extract_erosita_tile_21.log &
nohup python extract_erosita_tile.py 22 > logs/extract_erosita_tile_22.log &
nohup python extract_erosita_tile.py 23 > logs/extract_erosita_tile_23.log &

# after the scripts above are finished
python merge_erosita_tile.py

# after the merging is finished (it deletes temporary files !)
python clean_erosita_tile.py
















# what to simulate

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



