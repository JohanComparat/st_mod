"""
This script resambles https://gitlab.mpcdf.mpg.de/joco/erosita_sxrbg/-/blob/main/esass/erassX_write_scripts.py
It creates 3 scripts:
1) create images
2) start detection: 3 loops of erbox and erbackmap
3) finish detection: ermldet, apetool, srctool
/data26s/mpecl/eRASS1/??????/c946
/data26s/mpecl/eRASS1/358144/c946
/data26s/mpecl/eRASS1/358144/c946/*events*

nohup bash /data26s/comparat/simulations/erosim/eRASS1_UNIT_fA1i_2021_10_12_SKYMAP/086066/086066_pipeline_img.sh > /data26s/comparat/simulations/erosim/eRASS1_UNIT_fA1i_2021_10_12_SKYMAP/086066/086066_pipeline_img_RS.log &
nohup bash /data26s/comparat/simulations/erosim/eRASS1_UNIT_fA1i_2021_10_12_SKYMAP/086066/086066_pipeline_det1.sh > /data26s/comparat/simulations/erosim/eRASS1_UNIT_fA1i_2021_10_12_SKYMAP/086066/086066_pipeline_det1_RS.log &
nohup bash /data26s/comparat/simulations/erosim/eRASS1_UNIT_fA1i_2021_10_12_SKYMAP/086066/086066_pipeline_det2.sh > /data26s/comparat/simulations/erosim/eRASS1_UNIT_fA1i_2021_10_12_SKYMAP/086066/086066_pipeline_det2_RS.log &

nohup bash /data26s/comparat/simulations/erosim/eRASS1_UNIT_fA1i_2021_10_12_SKYMAP/092084/092084_pipeline_img.sh > /data26s/comparat/simulations/erosim/eRASS1_UNIT_fA1i_2021_10_12_SKYMAP/092084/092084_pipeline_img_RS.log &
nohup bash /data26s/comparat/simulations/erosim/eRASS1_UNIT_fA1i_2021_10_12_SKYMAP/092084/092084_pipeline_det1.sh > /data26s/comparat/simulations/erosim/eRASS1_UNIT_fA1i_2021_10_12_SKYMAP/092084/092084_pipeline_det1_RS.log &
nohup bash /data26s/comparat/simulations/erosim/eRASS1_UNIT_fA1i_2021_10_12_SKYMAP/092084/092084_pipeline_det2.sh > /data26s/comparat/simulations/erosim/eRASS1_UNIT_fA1i_2021_10_12_SKYMAP/092084/092084_pipeline_det2_RS.log &
"""
# !/usr/bin/env python
import sys, os, glob
import numpy as n
from astropy.table import Table, vstack
import astropy.io.fits as fits
sky_map_hdu = Table.read(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS.fits') )

GE_name = sys.argv[1] # 'sim_evt_e4_merge'

for sky_tile in sky_map_hdu[(sky_map_hdu['OWNER']==2)|(sky_map_hdu['OWNER']==0)][:1]:

    sky_tile_id = str(sky_tile['SRVMAP'])
    str_field = str(sky_tile['SRVMAP']).zfill(6)

    # directory
    #eRASSn = sys.argv[1]
    field_id = str_field # sys.argv[2]
    indir = os.path.join("/home/idies/workspace/erosim/Uchuu/LCerass/", field_id, GE_name)
    outdir = os.path.join(indir, 'eSASS')
    os.system('mkdir -p ' + outdir)
    # sb05_030108_020_FlareGTI_c030.fits  sb05_030108_020_Image_c030.fits.gz
    outprefix = field_id + "_"  # ""
    print(outdir)
    print(sys.argv)
    # input files
    EvtFiles = os.path.join(indir, 'evt_'+field_id+'.fits')
    print(EvtFiles)
    # SINGLE BAND
    # single band image
    out_im1 = os.path.join(outdir, outprefix + 'pipeline_img1.sh')
    # single band detection
    out_det1 = os.path.join(outdir, outprefix + 'pipeline_det1.sh')
    # ermldet + sensmap + source tool single band
    out_src1 = os.path.join(outdir, outprefix + 'pipeline_Src1.sh')
    print('sh '+out_im1)
    print('sh '+out_det1)
    print('sh '+out_src1)
    #
    Sixte = False
    # version
    # keep same naming convention as https://wiki.mpe.mpg.de/eRosita/PipeProc947
    #Three band detection:
    #1: 0.2 - 0.6 keV
    #2: 0.6 - 2.3 keV
    #3: 2.3 - 5.0 keV
    #Single band detection:
    #4: 0.2 - 2.3

    VerBand1 = str(1)
    VerBand2 = str(2)
    VerBand3 = str(3)
    VerBand4 = str(4)

    VerBands = n.array([VerBand1, VerBand2, VerBand3])#, VerBand4])

    # energy bands
    eminstr1, emaxstr1 = str(0.2), str(0.6)
    eminstr2, emaxstr2 = str(0.6), str(2.3)
    eminstr3, emaxstr3 = str(2.3), str(5.0)
    eminstr4, emaxstr4 = str(0.2), str(2.3)

    energy_band_mins = n.array([ eminstr1, eminstr2, eminstr3])#, eminstr4 ])
    energy_band_maxs = n.array([ emaxstr1, emaxstr2, emaxstr3])#, emaxstr4 ])


    # single band image
    f_out = open(out_im1, 'w')
    f_out.write("""#!/bin/bash/ \n""")

    #f_out.write(" \n ")
    #f_out.write("source /home/erosita/sw/eSASSusers_240410/bin/esass-init.sh \n")

    f_out.write(" \n ")
    f_out.write(" \n ")

    VerBand = VerBand4
    eminstr = eminstr4
    emaxstr = emaxstr4
    # output files
    f_out.write(f"# Create image, expoure maps for band {VerBand}\n")

    EvtImgFiles = os.path.join(outdir, f"{outprefix}02{VerBand}_EvtImg.fits")
    ExpMapFiles = os.path.join(outdir, f"{outprefix}02{VerBand}_ExpMap.fits")
    UnVigExpMap = os.path.join(outdir, f"{outprefix}02{VerBand}_UnvExpMap.fits")
    DetMask = os.path.join(outdir, f"{outprefix}02{VerBand}_DetMsk.fits")

    f_out.write("\n")
    f_out.write("# radec2xy  \n")
    f_out.write("\n")
    # radec2xy file=tmp_${hid}.fits ra0=${ra_cen} dec0=${dec_cen}
    cmd_radec = "radec2xy file="+EvtFiles+" ra0="+str(sky_tile['RA_CEN'])+" dec0="+str(sky_tile['DE_CEN'])
    f_out.write(cmd_radec)

    f_out.write("\n")
    f_out.write("# EVTOOL  \n")
    f_out.write("\n")
    f_out.write("# EVTOOL creates event image " + os.path.basename(EvtImgFiles) + " \n")
    f_out.write("\n")

    f_out.write("evtool \\\n")
    f_out.write("eventfiles=\"" + EvtFiles + "\" \\\n")
    f_out.write("outfile=\"" + EvtImgFiles + "\" \\\n")
    f_out.write("events=yes \\\n")
    f_out.write("image=yes \\\n")
    f_out.write("rebin=80 \\\n")
    f_out.write("size='3240 3240' \\\n")
    f_out.write("pattern=15 \\\n")
    f_out.write("center_position='0 0' \\\n")
    f_out.write("flag=0 \\\n")
    f_out.write("emin=" + """\"""" + eminstr + """\" \\\n""")
    f_out.write("emax=" + """\"""" + emaxstr + """\" \n """)

    f_out.write(" \n ")
    f_out.write(" \n ")


    f_out.write(" \n ")
    f_out.write(" \n ")

    f_out.write("\n")
    f_out.write("# EXPMAP c030 \n")
    f_out.write("\n")
    f_out.write("# EXPMAP creates exposure map " + os.path.basename(ExpMapFiles) + " \n")
    f_out.write("\n")

    f_out.write("expmap \\\n")
    f_out.write("inputdatasets=\"" + EvtFiles + "\" \\\n")
    f_out.write("templateimage=\"" + EvtImgFiles + "\"  \\\n")
    f_out.write("mergedmaps=\"" + ExpMapFiles + "\" \\\n")
    f_out.write("withvignetting=yes  \\\n")
    f_out.write("withmergedmaps=yes  \\\n")
    f_out.write("withsinglemaps=no   \\\n")
    f_out.write("withinputmaps=no   \\\n")
    f_out.write("withweights=yes     \\\n")
    f_out.write("plindex=-1.7     \\\n")
    f_out.write("emin=" + """\"""" + eminstr + """\" \\\n""")
    f_out.write("emax=" + """\"""" + emaxstr + """\" \n """)

    f_out.write("\n")
    f_out.write("# ERMASK c030 \n")
    f_out.write("\n")
    f_out.write("# ERMASK creates a mask " + os.path.basename(DetMask) + " \n")
    f_out.write("\n")

    f_out.write("ermask \\\n")
    f_out.write("expimage=\"" + ExpMapFiles + "\" \\\n")
    f_out.write("detmask=\"" + DetMask + "\"  \\\n")
    f_out.write("threshold1=0.0005 \\\n")
    f_out.write("threshold2=1.0 \\\n")
    f_out.write("regionfile_flag=no \n")

    f_out.write("\n")

    f_out.close()

    ##
    ##
    # DETECTION 1 bands
    ##
    ##

    f_out = open(out_det1, 'w')
    f_out.write("""#!/bin/bash/ \n""")

    #f_out.write("# source the new esass for detection \n")
    #f_out.write("source /home/erosita/sw/eSASSusers_240410/bin/esass-init.sh \n")

    MinDetLike = 5
    ERML_CUTRAD = str(15)
    ERML_MULTRAD = str(15)
    ERML_likemin = str(5)
    ERML_extlikemin = str(3)
    ERML_extmin = str(2)
    ERML_extmax = str(15)

    # energy band
    Nmlin = 3

    ecfstr1 = '1.028e+12'
    ecfstr2 = '1.087e+12'
    ecfstr3 = "1.147e+11"
    ecfstr4 = "1.074e+12"

    f_out.write("# Now do single band detection")
    VerBand = VerBand4
    eminstr, emaxstr = str(int(1000*float(eminstr4))), str(int(1000*float(emaxstr4)))
    ecfstr = ecfstr4
    #input files
    EvtImgFiles = os.path.join(outdir, f"{outprefix}02{VerBand}_EvtImg.fits")
    ExpMapFiles = os.path.join(outdir, f"{outprefix}02{VerBand}_ExpMap.fits")
    UnVigExpMap = os.path.join(outdir, f"{outprefix}02{VerBand}_UnvExpMap.fits")
    DetMask = os.path.join(outdir, f"{outprefix}02{VerBand}_DetMsk.fits")

    # output files

    SrcMapFiles = {n: os.path.join(outdir, f"{outprefix}02{VerBand}_Src{n:d}Map.fits") for n in (1, 2, 3)}
    BkgMapFiles = {n: os.path.join(outdir, f"{outprefix}02{VerBand}_Bg{n:d}Map.fits") for n in (1, 2, 3)}
    CheMskFiles = os.path.join(outdir, f"{outprefix}02{VerBand}_CheMsk.fits")

    BoxCats = {n: os.path.join(outdir, f"{outprefix}02{VerBand}_Bo{n:d}Cat.fits") for n in (1, 2, 3, 4, 5)}
    MLCats = {n: os.path.join(outdir, f"{outprefix}02{VerBand}_ML{n}Cat.fits") for n in (1, 2, 'A')} #uses only 1
    ExtCat = os.path.join(outdir, f"{outprefix}02{VerBand}_ExtCat.fits")
    MLinCat = {1: BoxCats[Nmlin], 2: os.path.join(outdir, f"{outprefix}02{VerBand}_ML1Cat.fits"),
            'A': os.path.join(outdir, f"{outprefix}02{VerBand}_ML2Cat.fits")}
    ApeCat = lambda eef: os.path.join(outdir, f"{outprefix}02{VerBand}_ApeCat_eef" + eef + ".fits")
    ApeSenFiles = lambda eef: os.path.join(outdir, f"{outprefix}02{VerBand}_ApeSen_eef" + eef + ".fits")
    PSFMapFiles = lambda eef: os.path.join(outdir, f"{outprefix}02{VerBand}_PSFMap_eef" + eef + ".fits")
    SrcCats = {n: os.path.join(outdir, f"{outprefix}02{VerBand}_Sc{n:d}Cat.fits") for n in (1, 2, 3)} #only use 1

    f_out.write("# BACKGROUND and BOX detection c030 \n")

    f_out.write("\n")
    f_out.write("# BACKGROUND and BOX detection. Iteration 1 c030 \n")
    f_out.write("\n")
    f_out.write("# ERBOX creates  " + os.path.basename(BoxCats[1]) + " \n")

    f_out.write("erbox \\\n")
    f_out.write("boxlist=\"" + BoxCats[1] + "\" \\\n")
    f_out.write("images=\"" + EvtImgFiles + "\" \\\n")
    f_out.write("expimages=\"" + ExpMapFiles + "\" \\\n")
    f_out.write("detmasks=\"" + DetMask + "\" \\\n")
    f_out.write("hrdef= \\\n"            )
    f_out.write("nruns=3 \\\n")
    f_out.write("likemin=6 \\\n")
    f_out.write("boxsize=4 \\\n")
    f_out.write("compress_flag=N \\\n")
    f_out.write("bkgima_flag=N \\\n")
    f_out.write("expima_flag=Y \\\n")
    f_out.write("detmask_flag=Y \\\n")
    f_out.write("ecf=" + """\"""" + ecfstr + """\" \\\n""")
    f_out.write("emin=" + """\"""" + eminstr + """\" \\\n""")
    f_out.write("emax=" + """\"""" + emaxstr + """\" \n """)

    f_out.write("\n")
    f_out.write("\n")
    f_out.write(
        "# erbackmap creates  " + os.path.basename(BkgMapFiles[1]) + " and " + os.path.basename(CheMskFiles) + " \n")
    f_out.write("\n")

    f_out.write("erbackmap \\\n")
    f_out.write("boxlist=\"" + BoxCats[1] + "\" \\\n")
    f_out.write("bkgimage=\"" + BkgMapFiles[1] + "\" \\\n")
    f_out.write("image=\"" + EvtImgFiles + "\" \\\n")
    f_out.write("expimage=\"" + ExpMapFiles + "\" \\\n")
    f_out.write("detmask=\"" + DetMask + "\" \\\n")
    f_out.write("cheesemask=\"" + CheMskFiles + "\" \\\n")
    f_out.write("usermask_flag=N \\\n")
    f_out.write("usermask=\\\n")
    f_out.write("idband=1 \\\n")
    f_out.write("scut=0.00005 \\\n")
    f_out.write("mlmin=6 \\\n")
    f_out.write("maxcut=0.5 \\\n")
    f_out.write("fitmethod=smooth \\\n")
    f_out.write("nsplinenodes=36 \\\n")
    f_out.write("degree=2 \\\n")
    f_out.write("smoothflag=Y \\\n")
    f_out.write("smoothval=15.0 \\\n")
    f_out.write("smoothmax=360 \\\n")
    f_out.write("snr=40 \\\n")
    f_out.write("excesssigma=1000. \\\n")
    f_out.write("nfitrun=3 cheesemask_flag=Y \\\n")
    f_out.write("detmask_flag=Y expima_flag=Y \\\n")
    f_out.write("""expima2_flag=N expimage2="" \\\n""")
    f_out.write("emin=" + """\"""" + eminstr + """\" \\\n""")
    f_out.write("emax=" + """\"""" + emaxstr + """\" \n """)

    f_out.write("\n")
    f_out.write("# BACKGROUND and BOX detection. Iteration 2 c030 \n")
    f_out.write("\n")
    f_out.write("# ERBOX creates  " + os.path.basename(BoxCats[2]) + ", it now uses " + os.path.basename(
        BkgMapFiles[1]) + " \n")
    f_out.write("\n")

    f_out.write("erbox \\\n")
    f_out.write("boxlist=\"" + BoxCats[2] + "\" \\\n")
    f_out.write("images=\"" + EvtImgFiles + "\" \\\n")
    f_out.write("expimages=\"" + ExpMapFiles + "\" \\\n")
    f_out.write("bkgimages=\"" + BkgMapFiles[1] + "\" \\\n")
    f_out.write("detmasks=\"" + DetMask + "\" \\\n")
    f_out.write("hrdef= \\\n"            )
    f_out.write("nruns=3 \\\n")
    f_out.write("likemin=4 \\\n")
    f_out.write("boxsize=4 \\\n")
    f_out.write("compress_flag=N \\\n")
    f_out.write("bkgima_flag=Y \\\n")
    f_out.write("expima_flag=Y \\\n")
    f_out.write("detmask_flag=Y \\\n")
    f_out.write("ecf=" + """\"""" + ecfstr + """\" \\\n""")
    f_out.write("emin=" + """\"""" + eminstr + """\" \\\n""")
    f_out.write("emax=" + """\"""" + emaxstr + """\" \n """)

    f_out.write("\n")
    f_out.write("\n")
    f_out.write(
        "# erbackmap creates  " + os.path.basename(BkgMapFiles[2]) + " and " + os.path.basename(CheMskFiles) + " \n")
    f_out.write("# first remove previous cheesemask, then run erbackmap command \n")
    f_out.write("\n")

    f_out.write("rm " + CheMskFiles + " \n")

    f_out.write("\n")

    f_out.write("erbackmap \\\n")
    f_out.write("boxlist=\"" + BoxCats[2] + "\" \\\n")
    f_out.write("bkgimage=\"" + BkgMapFiles[2] + "\" \\\n")
    f_out.write("image=\"" + EvtImgFiles + "\" \\\n")
    f_out.write("expimage=\"" + ExpMapFiles + "\" \\\n")
    f_out.write("detmask=\"" + DetMask + "\" \\\n")
    f_out.write("cheesemask=\"" + CheMskFiles + "\" \\\n")
    f_out.write("usermask_flag=N \\\n")
    f_out.write("usermask=\\\n")
    f_out.write("idband=1 \\\n")
    f_out.write("scut=0.00005 \\\n")
    f_out.write("mlmin=6 \\\n")
    f_out.write("maxcut=0.5 \\\n")
    f_out.write("fitmethod=smooth \\\n")
    f_out.write("nsplinenodes=36 \\\n")
    f_out.write("degree=2 \\\n")
    f_out.write("smoothflag=Y \\\n")
    f_out.write("smoothval=15.0 \\\n")
    f_out.write("smoothmax=360 \\\n")
    f_out.write("snr=40 \\\n")
    f_out.write("excesssigma=1000. \\\n")
    f_out.write("nfitrun=3 cheesemask_flag=Y \\\n")
    f_out.write("detmask_flag=Y expima_flag=Y \\\n")
    f_out.write("""expima2_flag=N expimage2="" \\\n""")
    f_out.write("emin=" + """\"""" + eminstr + """\" \\\n""")
    f_out.write("emax=" + """\"""" + emaxstr + """\" \n """)

    f_out.write("\n")
    f_out.write("# BACKGROUND and BOX detection. Iteration 3 c030 \n")
    f_out.write("\n")
    f_out.write("# ERBOX creates  " + os.path.basename(BoxCats[3]) + ", it now uses " + os.path.basename(
        BkgMapFiles[2]) + " \n")
    f_out.write("\n")

    f_out.write("erbox \\\n")
    f_out.write("boxlist=\"" + BoxCats[3] + "\" \\\n")
    f_out.write("images=\"" + EvtImgFiles + "\" \\\n")
    f_out.write("expimages=\"" + ExpMapFiles + "\" \\\n")
    f_out.write("bkgimages=\"" + BkgMapFiles[2] + "\" \\\n")
    f_out.write("detmasks=\"" + DetMask + "\" \\\n")
    f_out.write("hrdef= \\\n"            )
    f_out.write("nruns=3 \\\n")
    f_out.write("likemin=4 \\\n")
    f_out.write("boxsize=4 \\\n")
    f_out.write("compress_flag=N \\\n")
    f_out.write("bkgima_flag=Y \\\n")
    f_out.write("expima_flag=Y \\\n")
    f_out.write("detmask_flag=Y \\\n")
    f_out.write("ecf=" + """\"""" + ecfstr + """\" \\\n""")
    f_out.write("emin=" + """\"""" + eminstr + """\" \\\n""")
    f_out.write("emax=" + """\"""" + emaxstr + """\" \n """)

    f_out.write("\n")
    f_out.write("\n")
    f_out.write(
        "# erbackmap creates  " + os.path.basename(BkgMapFiles[3]) + " and " + os.path.basename(CheMskFiles) + " \n")
    f_out.write("# first remove previous cheesemask, then run erbackmap command \n")
    f_out.write("\n")

    f_out.write("rm " + CheMskFiles + " \n")
    f_out.write("\n")

    f_out.write("erbackmap \\\n")
    f_out.write("boxlist=\"" + BoxCats[3] + "\" \\\n")
    f_out.write("bkgimage=\"" + BkgMapFiles[3] + "\" \\\n")
    f_out.write("image=\"" + EvtImgFiles + "\" \\\n")
    f_out.write("expimage=\"" + ExpMapFiles + "\" \\\n")
    f_out.write("detmask=\"" + DetMask + "\" \\\n")
    f_out.write("cheesemask=\"" + CheMskFiles + "\" \\\n")
    f_out.write("usermask_flag=N \\\n")
    f_out.write("usermask=\\\n")
    f_out.write("idband=1 \\\n")
    f_out.write("scut=0.00005 \\\n")
    f_out.write("mlmin=6 \\\n")
    f_out.write("maxcut=0.5 \\\n")
    f_out.write("fitmethod=smooth \\\n")
    f_out.write("nsplinenodes=36 \\\n")
    f_out.write("degree=2 \\\n")
    f_out.write("smoothflag=Y \\\n")
    f_out.write("smoothval=15.0 \\\n")
    f_out.write("smoothmax=360 \\\n")
    f_out.write("snr=40 \\\n")
    f_out.write("excesssigma=1000. \\\n")
    f_out.write("nfitrun=3 cheesemask_flag=Y \\\n")
    f_out.write("detmask_flag=Y expima_flag=Y \\\n")
    f_out.write("""expima2_flag=N expimage2="" \\\n""")
    f_out.write("emin=" + """\"""" + eminstr + """\" \\\n""")
    f_out.write("emax=" + """\"""" + emaxstr + """\" \n """)

    f_out.write("\n")
    f_out.write("\n")

    #f_out.write("cd " + outdir + "\n")
    #f_out.write("\n")
    #f_out.write("chgrp erosim * \n")

    #f_out.write("\n")

    f_out.close()





    f_out = open(out_src1, 'w')
    f_out.write("""#!/bin/bash/ \n""")

    #f_out.write("# source the new esass for detection \n")
    #f_out.write("source /home/erosita/sw/eSASSusers_240410/bin/esass-init.sh \n")

    f_out.write("# DETECTION, ermldet c030 \n")

    f_out.write("\n")
    f_out.write(
        "# ermldet creates  " + os.path.basename(MLCats[1]) + " and " + os.path.basename(SrcMapFiles[1]) + " \n")
    f_out.write("# it starts from  " + os.path.basename(MLinCat[1]) + " and " + os.path.basename(
        BkgMapFiles[3]) + " and all other image products \n")

    f_out.write("\n")

    f_out.write("ermldet \\\n")
    f_out.write("boxlist=\"" + MLinCat[1] + "\" \\\n")
    f_out.write("mllist=\"" + MLCats[1] + "\" \\\n")
    f_out.write("images=\"" + EvtImgFiles + "\" \\\n")
    f_out.write("expimages=\"" + ExpMapFiles + "\" \\\n")
    f_out.write("bkgimages=\"" + BkgMapFiles[3] + "\" \\\n")
    f_out.write("detmasks=\"" + DetMask + "\" \\\n")
    f_out.write("srcimages=\"" + SrcMapFiles[1] + "\" \\\n")
    f_out.write("likemin=" + ERML_likemin + " \\\n")
    f_out.write("extlikemin=" + ERML_extlikemin + " \\\n")
    f_out.write("cutrad=" + ERML_CUTRAD + " \\\n")
    f_out.write("multrad=" + ERML_MULTRAD + " \\\n")
    f_out.write("compress_flag=N \\\n")
    f_out.write("extmin=" + ERML_extmin + " \\\n")
    f_out.write("extmax=" + ERML_extmax + " \\\n")
    f_out.write("expima_flag=Y  \\\n")
    f_out.write("detmask_flag=Y  \\\n")
    f_out.write("extentmodel=beta  \\\n")
    f_out.write("thres_flag=Y  \\\n")
    f_out.write("thres_col=scts  \\\n")
    f_out.write("thres_val=25.  \\\n")
    f_out.write("nmaxfit=4  \\\n")
    f_out.write("nmulsou=2  \\\n")
    f_out.write("fitext_flag=yes  \\\n")
    f_out.write("srcima_flag=Y  \\\n")
    f_out.write("shapelet_flag=Y \\\n")
    f_out.write("photon_flag=Y \\\n")
    if Sixte:
        f_out.write("sixte_flag=Y \\\n")
    else:
        f_out.write("sixte_flag=N \\\n")
    f_out.write("hrdef= \\\n")
    f_out.write("fitpos_flag=Y \\\n")
    f_out.write("twostage_flag=Y \\\n")
    f_out.write("""extlike_slope="0.0 0.0"  \\\n""")
    f_out.write("ecf=" + """\"""" + ecfstr + """\" \\\n""")
    f_out.write("emin=" + """\"""" + eminstr + """\" \\\n""")
    f_out.write("emax=" + """\"""" + emaxstr + """\" \n """)

    f_out.write("\n")
    f_out.write("\n")

    f_out.write("catprep infile=" + MLCats[1] + " outfile=" + SrcCats[1] + " \\\n")
    f_out.write("clobber=N \\\n")
    f_out.write("det_algo=ML \\\n")
    f_out.write("p_srcdensity=10 \\\n")
    f_out.write("owner=2 \\\n")
    f_out.write("flux_threshold=0.0 \\\n")
    f_out.write("detuid_prefix=eRO \n")

    f_out.write("\n")
    f_out.write("\n")

    f_out.write("# APETOOL c030 \n")
    f_out.write("# performs aperture photometry within 0.75 EEF of the PSF size \n")


    def cmd_apetool(eef):
        f_out.write("\n")
        f_out.write("# APETOOL creates  " + os.path.basename(ApeCat(eef)) + " and " + os.path.basename(
            ApeSenFiles(eef)) + " and " + os.path.basename(PSFMapFiles(eef)) + " \n")
        f_out.write("\n")

        f_out.write("apetool \\\n")
        f_out.write("mllist=\"" + MLCats[1] + "\" \\\n")
        f_out.write("apelist=\"" + MLCats[1] + "\" \\\n")
        f_out.write("apelistout=\"" + ApeCat(eef) + "\" \\\n")
        f_out.write("images=\"" + EvtImgFiles + "\" \\\n")
        f_out.write("expimages=\"" + ExpMapFiles + "\" \\\n")
        f_out.write("bkgimages=\"" + BkgMapFiles[3] + "\" \\\n")
        f_out.write("psfmaps=\"" + PSFMapFiles(eef) + "\" \\\n")
        f_out.write("detmasks=\"" + DetMask + "\" \\\n")
        f_out.write("apesenseimages=\"" + ApeSenFiles(eef) + "\" \\\n")
        f_out.write("srcimages=\"" + SrcMapFiles[1] + "\" \\\n")
        f_out.write("eindex=\"\" " + " \\\n")
        f_out.write("eefextract=" + eef + " \\\n")
        f_out.write("pthresh=4e-6" + " \\\n")
        f_out.write("cutrad=" + ERML_CUTRAD + " \\\n")
        f_out.write("psfmapsampling=11" + " \\\n")
        f_out.write("apexflag=yes" + " \\\n")
        f_out.write("stackflag=no" + " \\\n")

        f_out.write("psfmapflag=yes" + " \\\n")
        f_out.write("apesenseflag=yes" + " \\\n")
        f_out.write("shapepsf=yes" + " \\\n")
        f_out.write("emin=" + eminstr + " \\\n")
        f_out.write("emax=" + emaxstr + " \n")
        f_out.write("\n")
        #f_out.write("mv " + MLCats[1] + " " + ApeCat(eef) + "\n")

    # cmd_apetool("0.5")
    #cmd_apetool("0.6")
    cmd_apetool("0.75")
    # cmd_apetool("0.8")
    # cmd_apetool("1.5")
    # cmd_apetool("2.0")

    f_out.write("\n")


    f_out.write("# run srctool\n")
    f_out.write("export SRCTOOL_AUTOREG_BACK_TO_SRC_AREA_RATIO=200 \n "   )

    f_out.write('# do single band\n')

    SrcCat = os.path.join(outdir, f"{outprefix}02{VerBand4}_Sc1Cat.fits") # 192075_0216_Sc1Cat.fits.gz
    # output files
    SrcReg = os.path.join(outdir, f"{outprefix}_02_23_srcAUTO.reg")
    BkgReg = os.path.join(outdir, f"{outprefix}_02_23_bkgAUTO.reg")

    f_out.write("# source region extractions \n"   )

    f_out.write("\n"   )
    f_out.write("# srctool creates  "+os.path.basename(SrcReg)+" "+ os.path.basename(BkgReg) +" \n"   )

    f_out.write("srctool \\\n")
    f_out.write("eventfiles=\""+EvtFiles+"\" \\\n")
    f_out.write("srccoord=\""+SrcCat+"\" \\\n")
    f_out.write("todo=NOSRCGTI \\\n"          )
    f_out.write("""insts="1 2 3 4 5 6 7"  \\\n""")
    f_out.write("flagsel=0 \\\n"      )
    f_out.write("srcreg=AUTO \\\n"        )
    f_out.write("backreg=AUTO \\\n"        )
    f_out.write("psftype=2D_PSF \\\n"        )
    f_out.write("outsrcreg=\""+SrcReg+"\" \\\n")
    f_out.write("outbackreg=\""+BkgReg+"\" \n")


    f_out.write("\n"   )

    #f_out.write("cd " + outdir + "\n")
    #f_out.write("\n")
    #f_out.write("chgrp erosim * \n")

    #f_out.write("\n")

    f_out.close()


