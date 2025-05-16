# installation of sixte done on sciserver

# makes the image library for the galaxy cluster Uchuu run
export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/sixte_fgbg/

# ONGOING
# recast events into the erosita fields
~/workspace/erosim/sixte_output_events/eRASS8_stars
python extract_erosita_tile_stars.py # DONE
# # re-casts these files :
# simulated_photons_ccd1.fits
# simulated_photons_ccd2.fits
# simulated_photons_ccd3.fits
# simulated_photons_ccd4.fits
# simulated_photons_ccd5.fits
# simulated_photons_ccd6.fits
# simulated_photons_ccd7.fits

# ONGOING
# recast events into the erosita fields
# need to work per field, summary file would be too big
cd ~/workspace/erosim/sixte_output_events/eRASS8_ParticleBackground
rsync -avz joco@raven.mpcdf.mpg.de:/ptmp/joco/erosim/sixte_output_events/eRASS8_ParticleBackground/evt_particle_???.fits . # ONGOING
python extract_erosita_tile_BG.py  # Deprecated
python extract_erosita_tile_BG_v2.py  # DONE


#
# ONGOING, generation of FG/BG events by re-sampling
#
~/workspace/erosim/simput/bkg_erosita_simput_full_sky$
export UCHUU='/home/idies/workspace/erosim/Uchuu'
export GIT_STMOD='/home/idies/workspace/erosim/software/st_mod'
export GIT_STMOD_DATA='/home/idies/workspace/erosim/software/st_mod_data'
cd $GIT_STMOD/src/sixte_fgbg/
# event generation ongoing
# see readme_sixte_bg.sh
