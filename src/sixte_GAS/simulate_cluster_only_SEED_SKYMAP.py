"""
Simulate eROSITA event data from sky model using SIXTE.
Uses sixte version 2.7

To use sixte version 3.0, this needs a significant upgrade

"""
import subprocess as sp
from multiprocessing import Pool
from functools import partial
import os
import sys, glob
import numpy as np
from astropy.io import fits
from astropy.table import Table

#One iteration for subprocess
def one_iter_func(sky_tile, other_elements):
    """
    Single iteration of loop function over healpix pixels. Writes the files to path_2_eRO_catalog.

    Input:

    - sky_tile
    - other_elements: list containing basename, seed, LC_dir, erass_option, sixte_version
    """

    #Identify sky tile
    sky_tile_id = str(sky_tile['SRVMAP'])
    str_field = str(sky_tile['SRVMAP']).zfill(6)

    #Splits other_elements into relevant quantities
    basename, seed, LC_dir, erass_option, sixte_version, log_fn_for_check = other_elements
    print('\nTile {0} - seed: {1} LC_dir: {2} eRASS option: {3}'.format(str_field, seed, LC_dir, erass_option))

    #Loop as usual
    ra_cen = sky_tile['RA_CEN']
    dec_cen = sky_tile['DE_CEN']
    print('\nTile {0} - sky tile ID: {1}'.format(str_field, sky_tile_id))
    #Next line was problematic on 07.11.2025, rseppi added split_simput to handle simput files larger than 1000
    #simput_files = np.array([os.path.join(os.environ['UCHUU'], LC_dir, str_field, basename+'_simput_final.fits'), ])
    simput_file = os.path.join(os.environ['UCHUU'], LC_dir, str_field, basename+'_simput_final.fits')
    simput_files = split_simput(simput_file)
    data_dir = os.path.join(os.environ['UCHUU'], LC_dir, str_field, erass_option + "_sixte_"+ sixte_version+ "_SEED_"+str(seed).zfill(3) +"_events_cluster_"+basename )
    print('\nTile {0} - output path: {1}'.format(str_field, data_dir))

    #Select attitude file and time interval
    if erass_option=='Att_eRASS8':
        p_2_att_file = "/home/idies/workspace/erosim/erosita_attitude/eRASS_4yr_epc85_att.fits"
        t_starts = np.array([617943605])
        t_stops = np.array([744174005])

    #Iterate over steps
    for jj, (t0, t1) in enumerate(zip(t_starts, t_stops)):
        print('\nTile {0} - Step: {1} t_starts: {2} t_stops: {3}'.format(str_field, jj, t0, t1))
        bkg = 0
        t_start = t0 # 617943605.0
        exposure = t1 - t0 # 31536000 * 4 # = 4 years  # 31536000 = 1year # 15750000 = 1/2 year
        try:
            Simulator(
            jj,
            bkg,
            t_start,
            exposure,
            int(seed),
            simput_files,
            data_dir,
            ra_cen,
            dec_cen,
            p_2_att_file,
            str_field).run_all()
        except(FileNotFoundError):
            print('\nTile {0} - Step {1} - file is missing'.format(str_field, jj))

    #Save the tile_id to the logfile to mark that it is done
    with open(log_fn_for_check, 'a') as lffc:
        lffc.write('{0}\n'.format(str_field))    

# =============================================================================

#Splits SIMPUT file into multiple parts
def split_simput(input_simput, max_per=1000, flux_min=8e-16):
    """
    Split a SIMPUT file into multiple parts with at most `max_per` sources each.
    Returns a list of output file paths.
    """

    clu_simput = Table.read(input_simput)

    # Apply flux cut
    mask = clu_simput["FLUX"] > flux_min
    clu_simput = clu_simput[mask]

    n = len(clu_simput)

    # Output directory and base name
    base_dir = os.path.dirname(input_simput)
    base_name = os.path.basename(input_simput).replace(".fits", "")

    def make_hdul(subtab):
        """Build a SIMPUT HDUList from a subset table."""
        hdu_cols = fits.ColDefs([
            fits.Column(name="SRC_ID", format='K', unit='',
                        array=subtab["SRC_ID"].astype(int)),
            fits.Column(name="RA", format='D', unit='deg',
                        array=subtab["RA"]),
            fits.Column(name="DEC", format='D', unit='deg',
                        array=subtab["DEC"]),
            fits.Column(name="E_MIN", format='D', unit='keV',
                        array=np.ones(len(subtab)) * 0.5),
            fits.Column(name="E_MAX", format='D', unit='keV',
                        array=np.ones(len(subtab)) * 2.0),
            fits.Column(name="FLUX", format='D', unit='erg/s/cm**2',
                        array=subtab["FLUX"]),
            fits.Column(name="IMAGE", format='100A', unit='',
                        array=subtab["IMAGE"]),
            fits.Column(name="SPECTRUM", format='100A', unit='',
                        array=subtab["SPECTRUM"]),
            fits.Column(name="IMGROTA", format='D', unit='deg',
                        array=subtab["IMGROTA"]),
            fits.Column(name="IMGSCAL", format='D', unit='',
                        array=subtab["IMGSCAL"]),
        ])
        hdu = fits.BinTableHDU.from_columns(hdu_cols)
        hdu.name = 'SRC_CAT'
        hdu.header['HDUCLASS'] = 'HEASARC/SIMPUT'
        hdu.header['HDUCLAS1'] = 'SRC_CAT'
        hdu.header['HDUVERS'] = '1.1.0'
        hdu.header['RADESYS'] = 'FK5'
        hdu.header['EQUINOX'] = 2000.0
        return fits.HDUList([fits.PrimaryHDU(), hdu])

    outfiles = []

    # Case 1: everything fits into a single SIMPUT file
    if n <= max_per:
        hdul = make_hdul(clu_simput)
        out_path = os.path.join(
            base_dir,
            f"{base_name}_fluxcut.fits"
        )
        hdul.writeto(out_path, overwrite=True)
        outfiles.append(out_path)
        return outfiles

    # Case 2: need to split into multiple parts
    n_parts = int(np.ceil(n / float(max_per)))
    edges = (np.arange(n_parts + 1) * max_per).astype(int)
    edges[-1] = n  # ensure last edge matches exactly

    for i in range(n_parts):
        lo, hi = edges[i], edges[i + 1]
        subtab = clu_simput[lo:hi]

        hdul = make_hdul(subtab)
        out_path = os.path.join(
            base_dir,
            f"{base_name}_fluxcut_part{i}.fits"
        )
        hdul.writeto(out_path, overwrite=True)
        outfiles.append(out_path)

    return outfiles

# =============================================================================

#Deals with processes
def proprocess(command):

    process = sp.Popen(
        command,
        shell=True,
        stdout=sp.PIPE,
        stderr=sp.STDOUT,
        text=True,
    )

    for line in process.stdout:
        sys.stdout.write(line)

    process.wait()

# =============================================================================

#Simulator class to manage SIXTE
class Simulator:
    """
    SIXTE simulator for eROSITA observations.
    1. Compute GTI file for given simput
    2. Simulate eROSITA observations of simput, using GTI to speed things up.
    """

    def __init__(self, block_id, with_bkg_par, t_start, exposure, seed, simput, data_dir, ra_cen, dec_cen, p_2_att_file, str_field):
        """
        :param with_bkg_par: Simulate with particle background.
        :param t_start: Start time of simulation. Input units of [s]
        :param exposure: Length of time to simulate for after t_start
        :param seed: Seed for random number generator.
        :param simput: Simput files (ie. the sky model)
        """
        self._block_id = 't'+str(block_id)
        self._with_bkg_par = bool(with_bkg_par)
        self._t_start = float(t_start)  # secs
        self._exposure = float(exposure)
        self._seed = int(seed)
        self._simput = simput
        self._N_simputs = len(simput)
        self._data_dir = data_dir
        self._ra_cen = ra_cen
        self._dec_cen = dec_cen
        self._p_2_att_file = p_2_att_file
        self._str_field = str_field

    def make_event_directory(self):
        """
        Check for whether directory exists and create if not.
        """
        try:
            os.makedirs(self._data_dir)
        except OSError as e:
            print('\nTile {0} - Directory already exists'.format(self._str_field))

    def compute_gti(self, simput_ero, gti_file):
        """
        Compute the GTI (good time interval) file given the simput file.
        Use this as input to SIXTE call for reducing computational expense.
        """
        cmd = ["ero_vis",
                "GTIfile={0}".format(gti_file),
                "Simput={0}".format(simput_ero),
                "Exposure={0:f}".format(self._exposure),
                "Attitude={0}".format(self._p_2_att_file),
                "TSTART={0:f}".format(self._t_start),
                "dt=0.5",
                "visibility_range=1.0",
                "clobber=yes"
                ]

        command = " ".join(cmd)
        print('\nTile {0} - Compute GTI with command:\n{1}'.format(self._str_field, command))
#        os.system(command) #Offending line in NEW SciServer?
        proprocess(command)
        hd=fits.open(gti_file)
        number_of_lines = hd[1].header['NAXIS2']
        print('\nTile {0} - GTI file has {1} lines'.format(self._str_field, number_of_lines))
        if number_of_lines==0:
            cmd_del = ['rm', gti_file]
            command_del = " ".join(cmd_del)
            print('\nTile {0} - Removing GTI file with command:\n{1}'.format(self._str_field, command_del))
            proprocess(command_del)
            
    def run_sixte(self):
        """
        Launch erosim from python.
        """

        print('\nTile {0} - Number of simput files: {1}'.format(self._str_field, self._N_simputs))
        prefix = os.path.join(self._data_dir, self._block_id+"erass_" )
        if self._N_simputs==1:
            cmd = ["erosim",
                    "Simput={0}".format(self._simput[0]),
                    "Prefix={0}".format(prefix),
                    "Attitude={0}".format(self._p_2_att_file),
                    "RA={0}".format(self._ra_cen),
                    "Dec={0}".format(self._dec_cen),
                    "GTIFile={0}/erass.gti".format(self._data_dir),
                    "TSTART={0}".format(self._t_start),
                    "Exposure={0}".format(self._exposure),
                    "MJDREF=51543.875",
                    "dt=0.5",
                    "Seed={0}".format(self._seed),
                    "clobber=yes",
                    "chatter=3"
                    ]
            if self._with_bkg_par:
                cmd.append("Background=yes")
            else:
                cmd.append("Background=no")

            command = " ".join(cmd)
            print('\nTile {0} - Running SIXTE with command:\n{1}'.format(self._str_field, command))
            proprocess(command)
            return  # we're done

        elif (self._N_simputs >= 2) and (self._N_simputs <= 6):
            cmd = ["erosim", "Simput={0}".format(self._simput[0])]
            for jj in np.arange(len(self._simput))[1:]:
                cmd.append("Simput"+str(int(jj+1))+"={0}".format(self._simput[jj]))
            cmd_end = ["Prefix={0}".format(prefix),
                    "Attitude={0}".format(self._p_2_att_file),
                    "RA={0}".format(self._ra_cen),
                    "Dec={0}".format(self._dec_cen),
                    "GTIFile={0}/erass.gti".format(self._data_dir),
                    "TSTART={0}".format(self._t_start),
                    "Exposure={0}".format(self._exposure),
                    "MJDREF=51543.875",
                    "dt=0.5",
                    "Seed={0}".format(self._seed),
                    "clobber=yes",
                    "chatter=3"
                    ]
            for el in cmd_end:
                cmd.append(el)

            if self._with_bkg_par is True:
                cmd.append("Background=yes")
            else:
                cmd.append("Background=no")
            command = " ".join(cmd)
            print('\nTile {0} - Running SIXTE with command:\n{1}'.format(self._str_field, command))
            proprocess(command)

        # ---------- >6 SIMPUTS: two runs ----------
        elif self._N_simputs > 6:
            # ---- First run: first 6 SIMPUTs ----
            cmd = ["erosim", "Simput={0}".format(self._simput[0])]
            # add Simput2..Simput6 with indices 1..5
            for jj in np.arange(1, 6):
                cmd.append("Simput" + str(int(jj + 1)) + "={0}".format(self._simput[jj]))

            cmd_end = [
                "Prefix={0}".format(prefix),
                "Attitude={0}".format(self._p_2_att_file),
                "RA={0}".format(self._ra_cen),
                "Dec={0}".format(self._dec_cen),
                "GTIFile={0}/erass.gti".format(self._data_dir),
                "TSTART={0}".format(self._t_start),
                "Exposure={0}".format(self._exposure),
                "MJDREF=51543.875",
                "dt=0.5",
                "Seed={0}".format(self._seed),
                "clobber=yes",
                "chatter=3"
            ]
            cmd.extend(cmd_end)

            if self._with_bkg_par:
                cmd.append("Background=yes")
            else:
                cmd.append("Background=no")

            command = " ".join(cmd)
            print('\nTile {0} - FIRST SIXTE RUN (SIMPUT 1–6) - Running SIXTE with command:\n{1}'.format(self._str_field, command))
            proprocess(command)

            # ---- Second run: remaining SIMPUTs (6,7,...) ----
            remaining = self._simput[6:]
            prefix2 = os.path.join(self._data_dir, self._block_id + "erass_run2_")

            # first file in this run is "Simput="
            cmd1 = ["erosim", "Simput={0}".format(remaining[0])]
            # then Simput2, Simput3, ... for the rest
            for k, simput_path in enumerate(remaining[1:], start=2):
                cmd1.append("Simput{0:d}={1}".format(k, simput_path))

            cmd_end1 = [
                "Prefix={0}".format(prefix2),
                "Attitude={0}".format(self._p_2_att_file),
                "RA={0}".format(self._ra_cen),
                "Dec={0}".format(self._dec_cen),
                "GTIFile={0}/erass.gti".format(self._data_dir),
                "TSTART={0}".format(self._t_start),
                "Exposure={0}".format(self._exposure),
                "MJDREF=51543.875",
                "dt=0.5",
                "Seed={0}".format(self._seed),
                "clobber=yes",
                "chatter=3"
            ]
            cmd1.extend(cmd_end1)

            if self._with_bkg_par:
                cmd1.append("Background=yes")
            else:
                cmd1.append("Background=no")

            command1 = " ".join(cmd1)
            print('\nTile {0} - SECOND SIXTE RUN (SIMPUT 7+) - Running SIXTE with command:\n{1}'.format(self._str_field, command1))
            proprocess(command1)

    #This is what is called
    def run_all(self):
        """
        Run SIXTE simulation of eRASS 8
        """

        #Create directory for output
        print('\nTile {0} - Make event directory'.format(self._str_field))
        self.make_event_directory()
        
        #Compute GTI file
        print('\nTile {0} - Compute gti with ero_vis'.format(self._str_field))
        if self._N_simputs==1:
            path_to_gti = os.path.join(self._data_dir, "erass.gti")
            self.compute_gti(self._simput[0], path_to_gti)

        if self._N_simputs>=2:
            for el in self._simput:
                gti_file = os.path.join(self._data_dir, "erass_" + os.path.basename(el)[:-5] + ".gti")
                self.compute_gti(el, gti_file)
            print('\nTile {0} - Merging gti files'.format(self._str_field))
            path_to_gti_list = os.path.join(self._data_dir, "gti.list")
            path_to_gti = os.path.join(self._data_dir, "erass.gti")
            gti_files = glob.glob(self._data_dir + "/*.gti")
            if len(gti_files)>0:
                command_list = "ls " + self._data_dir + "/*.gti > " + path_to_gti_list
                print('\nTile {0} - Listing GTI files:\n{1}'.format(self._str_field, command_list))
                proprocess(command_list)
                command_merge = "mgtime ingtis=@ " + path_to_gti_list + " outgti=" + path_to_gti + " merge=OR"
                print('\nTile {0} - Merging GTI files with command:\n{1}'.format(self._str_field, command_merge))
                proprocess(command_merge)

        if self._N_simputs==0:
            print('\nTile {0} - No simput'.format(self._str_field))
            return 0.

        #If a GTI file is present then run SIXTE        
        print('\nTile {0} - Using GTI file:\n{1}'.format(self._str_field, os.path.isfile(path_to_gti)))
        if os.path.isfile(path_to_gti):
            print('\nTile {0} - Running SIXTE with erosim'.format(self._str_field))
            self.run_sixte()

            list_2_delete = np.array(glob.glob(os.path.join(self._data_dir,  "*.gti" )))
            print('\nTile {0} - Deleting files:\n{1}'.format(self._str_field, list_2_delete))
            del_list = np.array([os.remove(el) for el in  list_2_delete])

if __name__ == '__main__':

    #Read skymap
    sky_map_hdu = Table.read(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS.fits'))
    sky_map_hdu['ID'] = np.arange(len(sky_map_hdu))

    #Load basename and seed
    basename = sys.argv[1]
    seed = int(sys.argv[2])

    LC_dir = 'LCerass'
    erass_option = "Att_eRASS8"
    sixte_version = 'v27'
    print('This run is:\n Seed {0}\n Lightcone {1}\n eRASS option {2}'.format(seed, LC_dir, erass_option))
    
    #Number of cores to run task and joblabel
    ncores = int(sys.argv[3])
    joblabel = int(sys.argv[4])

    #Log file - stores tile IDs of tiles already processed.
    log_fn_for_check = sys.argv[5]

    #Open logfile and read what's on it
    if os.path.exists(log_fn_for_check):
        av_tile_id_list_int = np.loadtxt(log_fn_for_check, dtype = int, unpack = True)
        av_tile_id_list = [str(x).zfill(6) for x in av_tile_id_list_int]
        #Cycle sky_tiles as usual
        tiles_to_consider = []
        for sky_tile in sky_map_hdu[(sky_map_hdu['OWNER']==2)|(sky_map_hdu['OWNER']==0)]:
            str_field = str(sky_tile['SRVMAP']).zfill(6)
            if str_field in av_tile_id_list:
                print('Sky tile ID {0} already done: skipping it.'.format(str_field))    
            else:
                tiles_to_consider.append(sky_tile)
    else:
        lfcf = open(log_fn_for_check, 'w')
        lfcf.close
        #Cycle sky_tiles as usual
        tiles_to_consider = []
        for sky_tile in sky_map_hdu[(sky_map_hdu['OWNER']==2)|(sky_map_hdu['OWNER']==0)]:
            str_field = str(sky_tile['SRVMAP']).zfill(6)
            tiles_to_consider.append(sky_tile)
    
    #Define function for multiprocessing
    onepool_func = partial(one_iter_func, other_elements = [basename, seed, LC_dir, erass_option, sixte_version, log_fn_for_check])

#    #Map to cores    
#    with Pool(ncores) as p:
#        if joblabel == 1:
#            p.map(onepool_func, tiles_to_consider[:int(round((len(tiles_to_consider)-1)/2))])
#        elif joblabel == 2:
#            p.map(onepool_func, tiles_to_consider[int(round((len(tiles_to_consider)-1)/2)):])
#        elif joblabel == 3:
#            p.map(onepool_func, tiles_to_consider)
    if joblabel == 3:
        for ti in tiles_to_consider:
            onepool_func[ti]




    
    
    
    
    
    
    
    
    
    

        
