
# on laptop :
cd $GIT_SXRBG/sixte
python make_BG_simput_spectrum.py
# creates a spectrum here : /data40s/erosim/eRASS/bkg_erosita_simput_full_sky/spectra/mean_bg_spectrum.fits
python create_images.py
# creates images here :
# ls /data40s/erosim/eRASS/bkg_erosita_simput_full_sky/images/ | wc -l
# 196608
python map_to_simput_catalogue.py
# simput catalogues are here :
# /data40s/erosim/eRASS/bkg_erosita_simput_full_sky

# on DS43
# execute commands in 
cd $GIT_SXRBG/sixte
nohup sh run_bg_sixte_1.sh > run_bg_sixte_1.log &
nohup sh run_bg_sixte_2.sh > run_bg_sixte_2.log &
nohup sh run_bg_sixte_3.sh > run_bg_sixte_3.log &
nohup sh run_bg_sixte_4.sh > run_bg_sixte_4.log &
nohup sh run_bg_sixte_5.sh > run_bg_sixte_5.log &
nohup sh run_bg_sixte_6.sh > run_bg_sixte_6.log &
nohup sh run_bg_sixte_7.sh > run_bg_sixte_7.log &

# event files are written here 
NSIDE = 64
/data26s/comparat/simulations/erosim/eRASS1_events_BG_FG_Jun2021/
NSIDE = 128
/data26s/comparat/simulations/erosim/eRASS1_events_BG_FG_Sept2021/

