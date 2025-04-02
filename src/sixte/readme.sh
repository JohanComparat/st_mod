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

#
# TODO
# re-cast cluster simput files into erosita tiles
# implement flux cut
import numpy as np
for jj in np.arange(0,47,1):
    f=open('extract_erosita_tile_'+str(jj).zfill(3)+'.sh', 'w')
    f.write('#!/bin/bash \n')
    for kk in np.arange(100*jj, 100*(jj+1), 1):
        f.write("python extract_erosita_tile.py "+str(kk)+'\n')
    f.close()
    print('nohup sh extract_erosita_tile_'+str(jj).zfill(3)+'.sh > logs/extract_erosita_tile_'+str(jj).zfill(3)+'.log &')

nohup sh extract_erosita_tile_000.sh > logs/extract_erosita_tile_000.log &
nohup sh extract_erosita_tile_001.sh > logs/extract_erosita_tile_001.log &
nohup sh extract_erosita_tile_002.sh > logs/extract_erosita_tile_002.log &
nohup sh extract_erosita_tile_003.sh > logs/extract_erosita_tile_003.log &
nohup sh extract_erosita_tile_004.sh > logs/extract_erosita_tile_004.log &
nohup sh extract_erosita_tile_005.sh > logs/extract_erosita_tile_005.log &
nohup sh extract_erosita_tile_006.sh > logs/extract_erosita_tile_006.log &
nohup sh extract_erosita_tile_007.sh > logs/extract_erosita_tile_007.log &
nohup sh extract_erosita_tile_008.sh > logs/extract_erosita_tile_008.log &
nohup sh extract_erosita_tile_009.sh > logs/extract_erosita_tile_009.log &
nohup sh extract_erosita_tile_010.sh > logs/extract_erosita_tile_010.log &
nohup sh extract_erosita_tile_011.sh > logs/extract_erosita_tile_011.log &
nohup sh extract_erosita_tile_012.sh > logs/extract_erosita_tile_012.log &
nohup sh extract_erosita_tile_013.sh > logs/extract_erosita_tile_013.log &
nohup sh extract_erosita_tile_014.sh > logs/extract_erosita_tile_014.log &
nohup sh extract_erosita_tile_015.sh > logs/extract_erosita_tile_015.log &
nohup sh extract_erosita_tile_016.sh > logs/extract_erosita_tile_016.log &
nohup sh extract_erosita_tile_017.sh > logs/extract_erosita_tile_017.log &
nohup sh extract_erosita_tile_018.sh > logs/extract_erosita_tile_018.log &
nohup sh extract_erosita_tile_019.sh > logs/extract_erosita_tile_019.log &
nohup sh extract_erosita_tile_020.sh > logs/extract_erosita_tile_020.log &
nohup sh extract_erosita_tile_021.sh > logs/extract_erosita_tile_021.log &
nohup sh extract_erosita_tile_022.sh > logs/extract_erosita_tile_022.log &
nohup sh extract_erosita_tile_023.sh > logs/extract_erosita_tile_023.log &
nohup sh extract_erosita_tile_024.sh > logs/extract_erosita_tile_024.log &
nohup sh extract_erosita_tile_025.sh > logs/extract_erosita_tile_025.log &
nohup sh extract_erosita_tile_026.sh > logs/extract_erosita_tile_026.log &
nohup sh extract_erosita_tile_027.sh > logs/extract_erosita_tile_027.log &
nohup sh extract_erosita_tile_028.sh > logs/extract_erosita_tile_028.log &
nohup sh extract_erosita_tile_029.sh > logs/extract_erosita_tile_029.log &
nohup sh extract_erosita_tile_030.sh > logs/extract_erosita_tile_030.log &
nohup sh extract_erosita_tile_031.sh > logs/extract_erosita_tile_031.log &
nohup sh extract_erosita_tile_032.sh > logs/extract_erosita_tile_032.log &
nohup sh extract_erosita_tile_033.sh > logs/extract_erosita_tile_033.log &
nohup sh extract_erosita_tile_034.sh > logs/extract_erosita_tile_034.log &
nohup sh extract_erosita_tile_035.sh > logs/extract_erosita_tile_035.log &
nohup sh extract_erosita_tile_036.sh > logs/extract_erosita_tile_036.log &
nohup sh extract_erosita_tile_037.sh > logs/extract_erosita_tile_037.log &
nohup sh extract_erosita_tile_038.sh > logs/extract_erosita_tile_038.log &
nohup sh extract_erosita_tile_039.sh > logs/extract_erosita_tile_039.log &
nohup sh extract_erosita_tile_040.sh > logs/extract_erosita_tile_040.log &
nohup sh extract_erosita_tile_041.sh > logs/extract_erosita_tile_041.log &
nohup sh extract_erosita_tile_042.sh > logs/extract_erosita_tile_042.log &
nohup sh extract_erosita_tile_043.sh > logs/extract_erosita_tile_043.log &
nohup sh extract_erosita_tile_044.sh > logs/extract_erosita_tile_044.log &
nohup sh extract_erosita_tile_045.sh > logs/extract_erosita_tile_045.log &
nohup sh extract_erosita_tile_046.sh > logs/extract_erosita_tile_046.log &





















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



