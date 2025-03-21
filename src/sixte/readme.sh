# installation of sixte done on sciserver

# makes the image library for the galaxy cluster Uchuu run
export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/sixte
python create_image_library_hotgas.py # DONE

# images and spectra are here :
/home/idies/workspace/erosim/Uchuu/cluster_images
/home/idies/workspace/erosim/Uchuu/cluster_Xspectra

#
# TODO
# re-cast cluster simput files into erosita tiles
# implement flux cut
#

#
# TODO
# re-cast cluster simput files into erosita tiles
# implement flux cut
#






















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



