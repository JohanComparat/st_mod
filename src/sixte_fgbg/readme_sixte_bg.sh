
# on ds machines, done for the eRASS:1 twin :
cd $GIT_STMOD/src/sixte_fgbg/
python make_BG_simput_spectrum.py # DONE (previously)
# creates a spectrum here : ~/workspace/erosim/simput/bkg_erosita_simput_full_sky/spectra/mean_bg_spectrum.fits
python create_images.py # DONE (previously)
# creates images here :
# ls ~/workspace/erosim/simput/bkg_erosita_simput_full_sky/images/ | wc -l
# 196608 images
python map_to_simput_catalogue.py # DONE (previously)
# simput catalogues are here :
~/workspace/erosim/simput/bkg_erosita_simput_full_sky
catalog.fits has the full map


# on sciserver
# execute commands in 
cd $GIT_STMOD/src/sixte_fgbg/
nohup sh run_bg_sixte_1.sh > run_bg_sixte_1.log & # DONE
nohup sh run_bg_sixte_2.sh > run_bg_sixte_2.log & # DONE
nohup sh run_bg_sixte_3.sh > run_bg_sixte_3.log & # DONE
nohup sh run_bg_sixte_4.sh > run_bg_sixte_4.log & # DONE
nohup sh run_bg_sixte_5.sh > run_bg_sixte_5.log & # DONE
nohup sh run_bg_sixte_6.sh > run_bg_sixte_6.log & # DONE
nohup sh run_bg_sixte_7.sh > run_bg_sixte_7.log & # DONE
# simulated events are here :
ls /home/idies/workspace/erosim/eRASS8_events_BG_FG/???/t0erass_ccd1_evt.fits | wc -l

# re cast into the skymap fields for the eRASS light cone. DONE
python extract_erosita_tile_BG_v2.py  # DONE
# in the pBG2 directories

# DONE
