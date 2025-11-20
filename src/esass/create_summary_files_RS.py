import sys, os, glob
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table, join, vstack, Column
from astropy.coordinates import SkyCoord
from sklearn.neighbors import BallTree
from tqdm import tqdm
from astropy.cosmology import FlatLambdaCDM
import healpy as hp
import xspec
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoUCHUU = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
h = 0.6774
cosmo = cosmoUCHUU

'''
Creates summmary files in /home/idies/workspace/erosim/Uchuu/LCerass/SummaryFiles
takes the desired experiment directory as input, for example:

python create_summary_files_RS.py GE_e4_merge_AGNseed001_SimBKG_CLUseed001

R. Seppi (19.11.2025)
'''

deg_to_rad = np.pi/180.
cosmo = FlatLambdaCDM(H0=67.77, Om0=0.29)

basedir = '/home/idies/workspace/erosim/Uchuu/LCerass'

subdir = sys.argv[1]

print('Directory:', subdir)

p_2_match_cats = np.sort(glob.glob(os.path.join(basedir, '??????', subdir, 'Sc1*Uniq*fits')))
p_2_match_cats_Lext0 = np.sort(glob.glob(os.path.join(basedir, '??????', subdir, 'Sc_Lext0*Uniq*fits')))
p_2_clu_matched = os.path.join(basedir, 'SummaryFiles', 'LC_eRASS_Clusters_Input_'+subdir+'.fits')
p_2_eSASS_matched = os.path.join(basedir, 'SummaryFiles', 'LC_eRASS_eSASS_Output_'+subdir+'.fits')

p_2_clu_matched_Lext0 = os.path.join(basedir, 'SummaryFiles', 'LC_eRASS_Clusters_Input_'+subdir+'_Lext0.fits')
p_2_eSASS_matched_Lext0 = os.path.join(basedir, 'SummaryFiles', 'LC_eRASS_eSASS_Output_'+subdir+'_Lext0.fits')

coolfunc=np.loadtxt('/home/idies/workspace/erosim/software/st_mod_data/data/models/model_GAS/coolfunc.dat')
xgrid_ext = np.loadtxt('/home/idies/workspace/erosim/software/st_mod_data/data/models/model_GAS/radial_binning.txt')
profiles = Table.read('/home/idies/workspace/erosim/software/st_mod_data/data/models/model_GAS/profiles_010z015_1e14M2e14.fits')

def calc_lx(prof, kt, m5, z, fraction=1):
    """
    Compute the X-ray luminosity in the profile to be extended to 3x r500c.
    .. math::
        r_{500c} = \left(\\frac{3 M_{500c}}{ 4. \pi 500 \\rho_c(z)  }\\right)^{1/3} [ cm ]
    * profile\_emission = profile x rescale_factor
    * rescale\_factor = $\sqrt(kT/10.0) E^3(z)$
    * CF(kT) = cooling function, show the curve
    * L$_X(r)$ = $\Sigma_{<r}$( profile_emission $r_{500c}^2 2 \pi x CF(kT)$ Mpc=3.0856776e+24 dx )
    * L$_{500c}$ = L$_X$(1)
    """
    Mpc = 3.0856776e+24
    msun = 1.98892e33
    ez2 = cosmo.efunc(z) ** 2
    rhoc = cosmo.critical_density(z).value
    r500 = np.power(m5 * msun / 4. * 3. / np.pi / 500. / rhoc, 1. / 3.)
    resfact = np.sqrt(kt / 10.0) * np.power(ez2, 3. / 2.)
    prof_em = prof * resfact  # emission integral
    tlambda = np.interp(kt, coolfunc[:, 0], coolfunc[:, 1])  # cooling function
    dx = np.empty(len(xgrid_ext))
    dx[0] = xgrid_ext[0]
    dx[1:len(xgrid_ext)] = (np.roll(xgrid_ext, -1) - xgrid_ext)[:len(xgrid_ext) - 1]
    # print(prof_em*self.xgrid_ext*r500**2*2.*np.pi*tlambda*Mpc*dx)
    lxcum = np.cumsum(prof_em * xgrid_ext * r500 ** 2 * 2. * np.pi * tlambda * Mpc * dx)  # riemann integral
    lx_500 = np.interp(fraction, xgrid_ext, lxcum)  # evaluated at R500
    return lx_500




rsp_path="/home/idies/workspace/erosim/software/st_mod_data/data/onaxis_tm0_rmf_2023-01-17.fits"
arf_path="/home/idies/workspace/erosim/software/st_mod_data/data/survey_tm0_arf_filter_2023-01-17.fits"


def predict_ctr_xspec(Lx_erg_s, z, Z_solar, kT_keV,
                      rsp_path, arf_path,
                      band_rest_keV=(0.5, 2.0),
                      abund='aspl'):
    """
    Predict observed count rate (cts/s) for a tbabs*apec model normalized to
    intrinsic source-frame luminosity Lx using a single-run XSPEC
    """

    E1_rest, E2_rest = band_rest_keV

    # --- STEP 1: Calculate L_intrinsic for norm=1.0 ---
    # We use a temporary model file and log file
    xcm_file_1 = 'temp_Lx_1.xcm'
    log_file_1 = 'temp_Lx_1.log'

    with open(xcm_file_1, 'w') as f:
        # XSPEC Setup
        print('abund %s' % abund, file=f)

        # Define Model with norm=1.0 (This must be done first)
        print('model apec', file=f)
        print('%.4f' % kT_keV, file=f)  # kT
        print('%.4f' % Z_solar, file=f)  # Abundanc
        print('%.4f' % z, file=f)  # Redshift
        print('1.0', file=f)  # Norm = 1.0

        # Fakeit to load the response (RMF/ARF)
        print('fakeit none', file=f)
        print(rsp_path, file=f)
        print(arf_path, file=f)
        print('', file=f)  # No background file
        print('', file=f)  # No back exposure
        print('temp.pi', file=f)  # Temporary PHA file
        print('1.0', file=f)  # Exposure (1.0s)
        print('\n', file=f)

        # Calculate L_intrinsic for the current model (norm=1.0) and log it
        print('log %s' % log_file_1, file=f)
        # lumin E1 E2 z (This correctly gives the intrinsic Lx)
        print('lumin %.2f %.2f %.4f' % (E1_rest, E2_rest, z), file=f)
        print('log none', file=f)
        print('quit', file=f)
        print('y', file=f)

    # Execute XSPEC
    os.system(f'xspec - {xcm_file_1} > /dev/null 2>&1')

    # Read L_at_norm_1 from the log file
    L_at_norm_1 = 0.0
    try:
        with open(log_file_1, 'r') as f:
            for line in f:
                if 'Model Luminosity' in line:
                    L_at_norm_1 = float(line.split()[2])
                    break
    except FileNotFoundError:
        print("XSPEC log file not found. Check XSPEC installation/paths.")
        return 0.0
    # Clean up first run files
    os.system(f'rm {xcm_file_1} temp.pi {log_file_1}')

    if L_at_norm_1 <= 0.0:
        raise ValueError(f"Intrinsic luminosity calculated for norm=1.0 is zero ({L_at_norm_1:.2e}). Check kT/Z/band.")

    # Calculate final norm in Python
    final_norm = Lx_erg_s / L_at_norm_1

    # --- STEP 2: Rerun XSPEC with final norm and read rate ---
    xcm_file_2 = 'temp_rate_2.xcm'
    log_file_2 = 'temp_rate_2.log'

    with open(xcm_file_2, 'w') as f:
        # XSPEC Setup (same as before)
        print('abund %s' % abund, file=f)

        # Define Model with the new, scaled norm
        print('model apec', file=f)
        print('newpar 1 %.4f' % kT_keV, file=f)
        print('newpar 2 %.4f' % Z_solar, file=f)
        print('newpar 3 %.4f' % z, file=f)
        print('newpar 4 %.8e' % final_norm, file=f)  # Set the final, correct norm

        # Fakeit to load the response (RMF/ARF)
        print('fakeit none', file=f)
        print(rsp_path, file=f)
        print(arf_path, file=f)
        print('', file=f)
        print('', file=f)
        print('temp.pi', file=f)
        print('1.0', file=f)
        print('\n', file=f)

        # Show rate and log it
        print('log %s' % log_file_2, file=f)
        # This rate is now the final, absorbed count rate
        print('show rate', file=f)
        print('log none', file=f)

        print('quit', file=f)
        print('y', file=f)

    # Execute XSPEC
    os.system(f'xspec - {xcm_file_2} > /dev/null 2>&1')

    # Read the final rate
    rate_cts_per_s = 0.0
    try:
        with open(log_file_2, 'r') as f:
            for line in f:
                if 'Model predicted rate:' in line:
                    rate_cts_per_s = float(line.split()[4])
                    break
    except FileNotFoundError:
        print("Final XSPEC log file not found.")

    # Clean up all temporary files
    os.system(f'rm {xcm_file_2} temp.pi {log_file_2}')

    return rate_cts_per_s

#prepare to select unique area in esass catalogues
#sky_map_hdu = Table.read('SKYMAPS.fits')
try:
    sky_map_hdu = Table.read(os.path.join(os.environ['GIT_STMOD_DATA'],
                                          'data/models/eROSITA', 'SKYMAPS.fits'))
except:
    sky_map_hdu = Table.read('/home/idies/workspace/erosim/software/st_mod_data/data/models/eROSITA/SKYMAPS.fits')

def get_srvmap(ra, dec):
    return sky_map_hdu['SRVMAP'].value[(sky_map_hdu['RA_MIN'] < ra) & (
            sky_map_hdu['RA_MAX'] >= ra) & (sky_map_hdu['DE_MIN'] < dec) & (
                                               sky_map_hdu['DE_MAX'] >= dec)]


p2bkg_simput = os.path.join(os.environ['HOME'], 'workspace/erosim/simput/bkg_erosita_simput_full_sky', 'catalogue.fits')
bkg_simput = Table.read(p2bkg_simput, memmap=True)
sel_BG_feature = bkg_simput['FLUX']>2e-11
high_BKG_ids = bkg_simput['SRC_ID'][sel_BG_feature]
nside = 128

clu_all = []
eSASS_all = []
for p2uniq in tqdm(p_2_match_cats):
    tile = p2uniq.split('/')[-3]
    print(tile)
    Any = Table.read(os.path.join(basedir, tile, subdir, 'Sc1_'+tile+'_IDMatch_Any_Tot2.0.fits'), memmap=True)
    Uniq = Table.read(p2uniq, memmap=True)
    Uniq['ID_simput'] = Uniq['ID_simput'].astype(np.int32)
    Any['ID_simput'] = Any['ID_simput'].astype(np.int32)
    Uniq.rename_column('ID_simput', 'ID_uniq')
    Uniq['ID_Any'] = Any['ID_simput']

    clu = Table.read(os.path.join(basedir, tile, 'Xgas_bHS0.8_simput_final.fits'), memmap=True)
    XGAS = Table.read(os.path.join(basedir, tile, 'XGAS.fits'), memmap=True)

    # Boolean mask where XGAS['id'] is in clu['SRC_ID']
    mask = np.isin(XGAS['id'].astype(np.int64), clu['SRC_ID'].astype(np.int64))

    # Filtered XGAS table
    XGAS_filtered = XGAS[mask]
    R500c_arcmin = XGAS_filtered['R500c']/cosmo.kpc_proper_per_arcmin(z=XGAS_filtered['redshift_S'])
    XGAS_filtered['R500c_arcmin'] = R500c_arcmin.value

    #collect photon statistic
    apertures = np.array([0.1, 0.3, 0.5, 1., 2.]) # this will be times R500c
    p2evt_clu = os.path.join(basedir, tile, subdir, 'simCLUevt_'+tile+'.fits')
    p2evt_agn = os.path.join(basedir, tile, subdir, 'simAGNevt_'+tile+'.fits')
    p2evt_bkg = os.path.join(basedir, tile, subdir, 'simBKGevt_'+tile+'.fits')

    if os.path.isfile(p2evt_clu):
        clu_evt = Table.read(p2evt_clu, memmap=True)
        sel_evt = (clu_evt['SIGNAL'] > 0.2) & (clu_evt['SIGNAL'] < 2.3)
        clu_evt = clu_evt[sel_evt]

        coord_cat_CLU = deg_to_rad * np.transpose([XGAS_filtered['DEC'], XGAS_filtered['RA']])
        coord_evt_CLU = deg_to_rad * np.transpose([clu_evt['DEC'], clu_evt['RA']])
        Tree_CLU = BallTree(coord_evt_CLU, metric='haversine')

        Counts_02_23_CLU_CLU = np.empty((len(XGAS_filtered), len(apertures)), dtype=int)
        for j, ap in enumerate(apertures):
            r = deg_to_rad * ap * XGAS_filtered['R500c_arcmin'] /60. # deg
            Counts_02_23_CLU_CLU[:, j] = Tree_CLU.query_radius(coord_cat_CLU, r=r, count_only=True)

        XGAS_filtered.add_column(Column(name='COUNTS_02_23_CLU_CLU', data=Counts_02_23_CLU_CLU))

        if os.path.isfile(p2evt_agn):
            agn_evt = Table.read(p2evt_agn, memmap=True)
            sel_evt = (agn_evt['SIGNAL'] > 0.2) & (agn_evt['SIGNAL'] < 2.3)
            agn_evt = agn_evt[sel_evt]

            coord_evt_AGN = deg_to_rad * np.transpose([agn_evt['DEC'], agn_evt['RA']])
            Tree_AGN = BallTree(coord_evt_AGN, metric='haversine')

            Counts_02_23_CLU_AGN = np.empty((len(XGAS_filtered), len(apertures)), dtype=int)
            for j, ap in enumerate(apertures):
                r = deg_to_rad * ap
                Counts_02_23_CLU_AGN[:, j] = Tree_AGN.query_radius(coord_cat_CLU, r=r, count_only=True)
            XGAS_filtered.add_column(Column(name='COUNTS_02_23_CLU_AGN', data=Counts_02_23_CLU_AGN))

        if os.path.isfile(p2evt_bkg):
            bkg_evt = Table.read(p2evt_bkg, memmap=True)
            sel_evt = (bkg_evt['SIGNAL'] > 0.2) & (bkg_evt['SIGNAL'] < 2.3)
            bkg_evt = bkg_evt[sel_evt]

            coord_evt_BKG = deg_to_rad * np.transpose([bkg_evt['DEC'], bkg_evt['RA']])
            Tree_BKG = BallTree(coord_evt_BKG, metric='haversine')

            Counts_02_23_CLU_BKG = np.empty((len(XGAS_filtered), len(apertures)), dtype=int)
            for j, ap in enumerate(apertures):
                r = deg_to_rad * ap
                Counts_02_23_CLU_BKG[:, j] = Tree_BKG.query_radius(coord_cat_CLU, r=r, count_only=True)
            XGAS_filtered.add_column(Column(name='COUNTS_02_23_CLU_BKG', data=Counts_02_23_CLU_BKG))

    hdul_exp = fits.open(os.path.join(basedir, tile, subdir, 'eSASS', tile+'_024_ExpMap.fits'))
    ExpMap = hdul_exp[0].data
    hdul_Bg3 = fits.open(os.path.join(basedir, tile, subdir, 'eSASS', tile + '_024_Bg3Map.fits'))
    Bg3Map = hdul_Bg3[0].data
    wcs = WCS(hdul_exp[0].header)
    # Convert sky coordinates (RA, DEC) to pixel coordinates
    x_pix, y_pix = wcs.wcs_world2pix(XGAS_filtered['RA'], XGAS_filtered['DEC'], 0)
    x_pix = np.round(x_pix).astype(int)
    y_pix = np.round(y_pix).astype(int)

    # Clip to stay within image bounds
    ny, nx = ExpMap.shape
    x_pix = np.clip(x_pix, 0, nx - 1)
    y_pix = np.clip(y_pix, 0, ny - 1)

    # Extract exposure time at those pixel positions
    Texp = ExpMap[y_pix, x_pix]
    BgMap = Bg3Map[y_pix, x_pix]
    XGAS_filtered['ExpMap_eSASS'] = Texp
    XGAS_filtered['Bg3Map_eSASS'] = BgMap

    #Add eSASS detections
    eSASS_names = ['RA', 'DEC', 'DET_LIKE_0', 'EXT', 'EXT_LIKE', 'ML_CTS_0', 'ML_RATE_0', 'ML_FLUX_0']
    for col in eSASS_names:
        if col in ['RA', 'DEC']:
            XGAS_filtered[col+'_eSASS'] = np.ones(len(XGAS_filtered))*-99.
        else:
            XGAS_filtered[col] = np.ones(len(XGAS_filtered))*-99.

    # Create a mapping from ID to row index in Uniq
    id_to_index = {id_val: i for i, id_val in enumerate(Uniq['ID_uniq'])}

    # For each row in glist, check if id is in Uniq
    for i, gid in enumerate(XGAS_filtered['id'].astype('int32')):
        idx = id_to_index.get(gid)
        if idx is not None:
            for col in eSASS_names:
                # Rename RA and DEC
                target_col = col + '_eSASS' if col in ['RA', 'DEC'] else col
                XGAS_filtered[target_col][i] = Uniq[col][idx]

    eSASS_offset = np.ones(len(XGAS_filtered))*-99.
    seldet = XGAS_filtered['DET_LIKE_0']>0
    coord1 = SkyCoord(XGAS_filtered[seldet]['RA'], XGAS_filtered[seldet]['DEC'], unit='deg')
    coord2 = SkyCoord(XGAS_filtered[seldet]['RA_eSASS'], XGAS_filtered[seldet]['DEC_eSASS'], unit='deg')
    eSASS_offset[seldet] = coord1.separation(coord2).arcsecond
    XGAS_filtered['eSASS_offset'] = eSASS_offset
    XGAS_filtered['eSASS_offset'].units = 'arcsec'
    ipix = hp.ang2pix(nside, XGAS_filtered['RA'], XGAS_filtered['DEC'], nest=True, lonlat=True)
    in_good_region = ~np.isin(ipix, high_BKG_ids)
    XGAS_filtered['in_good_region'] = in_good_region
    clu_all.append(XGAS_filtered)



    # Now do eSASS
    SRV_value = np.array([get_srvmap(e0, e1) for e0, e1 in zip(Uniq['RA'], Uniq['DEC'])])
    sel_area = np.isin(SRV_value.T[0], int(tile))
    Uniq = Uniq[sel_area]
    # unique counterpart of AGN or Star
    PNT_sel =  (Uniq['ID_uniq'] >= 0) & (Uniq['ID_uniq'] < 1e5)
    # unique counterpart of cluster
    EXT_sel = (Uniq['ID_uniq'] >= 1e6)
    # no unique counterpart and linked to agn or star
    PNT2_sel = ((Uniq['ID_uniq'] < 0) & ((Uniq['ID_Any'] >= 0) & (Uniq['ID_Any'] < 1e5)))
    # no unique counterpart and linked to cluster
    EXT2_sel = (Uniq['ID_uniq'] < 0) & (Uniq['ID_Any'] >= 1e6)
    # BKG/FG feature
    # spurious, not linked to anything
    BKG_sel = (Uniq['ID_Any'] < 0)  # & (~sel_BG_feature)

    #Uniq['SRVMAP'] = SRV_value[sel_area]
    Uniq['class_PNT'] = PNT_sel
    Uniq['class_PNT2'] = PNT2_sel
    Uniq['class_EXT'] = EXT_sel
    Uniq['class_EXT2'] = EXT2_sel
    Uniq['class_BKG'] = BKG_sel
    ipix = hp.ang2pix(nside, Uniq['RA'], Uniq['DEC'], nest=True, lonlat=True)
    in_good_region = ~np.isin(ipix, high_BKG_ids)
    Uniq['in_good_region'] = in_good_region

    eSASS_all.append(Uniq)

if len(clu_all)>0:
    clu_all_table = vstack(clu_all)
    clu_all_table.write(p_2_clu_matched, overwrite=True)

    eSASS_all_table = vstack(eSASS_all)
    eSASS_all_table.write(p_2_eSASS_matched, overwrite=True)

    print('Written', p_2_eSASS_matched, 'and the basic version of',p_2_clu_matched)
    print('Now adding Lx<4r500...')
    CLUSTER_LX_soft_RF_4R500c = np.zeros(len(clu_all_table))
    for jj, cl in enumerate(tqdm(clu_all_table)):
        prof = profiles[cl['idx_profile']]['profiles']
        L500_prof = calc_lx(prof, cl['CLUSTER_kT'], cl['M500c'], cl['redshift_S'], fraction=1)
        resc_fact = np.power(10, cl['CLUSTER_LX_soft_RF_R500c']) / L500_prof
        CLUSTER_LX_soft_RF_4R500c[jj] = np.log10(
            resc_fact * calc_lx(prof, cl['CLUSTER_kT'], cl['M500c'], cl['redshift_S'], fraction=4))
    clu_all_table['CLUSTER_LX_soft_RF_4R500c'] = CLUSTER_LX_soft_RF_4R500c

    print('Now adding the unabsorbed count rate...')
    CTR_R500c = np.zeros(len(clu_all_table))
    CTR_4R500c = np.zeros(len(clu_all_table))
    for jj, cl in enumerate(tqdm(clu_all_table)):
        CTR_R500c[jj] = predict_ctr_xspec(
            Lx_erg_s=np.power(10, cl['CLUSTER_LX_soft_RF_R500c']),  # intrinsic Lx in rest 0.5–2 keV
            z=cl['redshift_S'],
            Z_solar=0.3,
            kT_keV=cl['CLUSTER_kT'],
            rsp_path=rsp_path,
            arf_path=arf_path,
            band_rest_keV=(0.5, 2.0)
        )

        CTR_4R500c[jj] = predict_ctr_xspec(
            Lx_erg_s=np.power(10, cl['CLUSTER_LX_soft_RF_4R500c']),  # intrinsic Lx in rest 0.5–2 keV
            z=cl['redshift_S'],
            Z_solar=0.3,
            kT_keV=cl['CLUSTER_kT'],
            rsp_path=rsp_path,
            arf_path=arf_path,
            band_rest_keV=(0.5, 2.0)
        )

    clu_all_table['CTR_R500c_unabs'] = CTR_R500c
    clu_all_table['CTR_4R500c_unabs'] = CTR_4R500c

    print('done!')

    print('Written', p_2_clu_matched, 'and', p_2_eSASS_matched, 'using', len(clu_all), 'tiles')
else:
    print('No tiles found to process!')






print('Now doing summary catalogues for the extlikemin=0 run...')
clu_all = []
eSASS_all = []
for p2uniq in tqdm(p_2_match_cats_Lext0):
    tile = p2uniq.split('/')[-3]
    print(tile)
    Any = Table.read(os.path.join(basedir, tile, subdir, 'Sc_Lext0'+tile+'_IDMatch_Any_Tot2.0.fits'), memmap=True)
    Uniq = Table.read(p2uniq, memmap=True)
    Uniq['ID_simput'] = Uniq['ID_simput'].astype(np.int32)
    Any['ID_simput'] = Any['ID_simput'].astype(np.int32)
    Uniq.rename_column('ID_simput', 'ID_uniq')
    Uniq['ID_Any'] = Any['ID_simput']

    clu = Table.read(os.path.join(basedir, tile, 'Xgas_bHS0.8_simput_final.fits'), memmap=True)
    XGAS = Table.read(os.path.join(basedir, tile, 'XGAS.fits'), memmap=True)

    # Boolean mask where XGAS['id'] is in clu['SRC_ID']
    mask = np.isin(XGAS['id'].astype(np.int64), clu['SRC_ID'].astype(np.int64))

    # Filtered XGAS table
    XGAS_filtered = XGAS[mask]
    R500c_arcmin = XGAS_filtered['R500c']/cosmo.kpc_proper_per_arcmin(z=XGAS_filtered['redshift_S'])
    XGAS_filtered['R500c_arcmin'] = R500c_arcmin.value

    #collect photon statistic
    apertures = np.array([0.1, 0.3, 0.5, 1., 2.]) # this will be times R500c
    p2evt_clu = os.path.join(basedir, tile, subdir, 'simCLUevt_'+tile+'.fits')
    p2evt_agn = os.path.join(basedir, tile, subdir, 'simAGNevt_'+tile+'.fits')
    p2evt_bkg = os.path.join(basedir, tile, subdir, 'simBKGevt_'+tile+'.fits')

    if os.path.isfile(p2evt_clu):
        clu_evt = Table.read(p2evt_clu, memmap=True)
        sel_evt = (clu_evt['SIGNAL'] > 0.2) & (clu_evt['SIGNAL'] < 2.3)
        clu_evt = clu_evt[sel_evt]

        coord_cat_CLU = deg_to_rad * np.transpose([XGAS_filtered['DEC'], XGAS_filtered['RA']])
        coord_evt_CLU = deg_to_rad * np.transpose([clu_evt['DEC'], clu_evt['RA']])
        Tree_CLU = BallTree(coord_evt_CLU, metric='haversine')

        Counts_02_23_CLU_CLU = np.empty((len(XGAS_filtered), len(apertures)), dtype=int)
        for j, ap in enumerate(apertures):
            r = deg_to_rad * ap * XGAS_filtered['R500c_arcmin'] /60. # deg
            Counts_02_23_CLU_CLU[:, j] = Tree_CLU.query_radius(coord_cat_CLU, r=r, count_only=True)

        XGAS_filtered.add_column(Column(name='COUNTS_02_23_CLU_CLU', data=Counts_02_23_CLU_CLU))

        if os.path.isfile(p2evt_agn):
            agn_evt = Table.read(p2evt_agn, memmap=True)
            sel_evt = (agn_evt['SIGNAL'] > 0.2) & (agn_evt['SIGNAL'] < 2.3)
            agn_evt = agn_evt[sel_evt]

            coord_evt_AGN = deg_to_rad * np.transpose([agn_evt['DEC'], agn_evt['RA']])
            Tree_AGN = BallTree(coord_evt_AGN, metric='haversine')

            Counts_02_23_CLU_AGN = np.empty((len(XGAS_filtered), len(apertures)), dtype=int)
            for j, ap in enumerate(apertures):
                r = deg_to_rad * ap
                Counts_02_23_CLU_AGN[:, j] = Tree_AGN.query_radius(coord_cat_CLU, r=r, count_only=True)
            XGAS_filtered.add_column(Column(name='COUNTS_02_23_CLU_AGN', data=Counts_02_23_CLU_AGN))

        if os.path.isfile(p2evt_bkg):
            bkg_evt = Table.read(p2evt_bkg, memmap=True)
            sel_evt = (bkg_evt['SIGNAL'] > 0.2) & (bkg_evt['SIGNAL'] < 2.3)
            bkg_evt = bkg_evt[sel_evt]

            coord_evt_BKG = deg_to_rad * np.transpose([bkg_evt['DEC'], bkg_evt['RA']])
            Tree_BKG = BallTree(coord_evt_BKG, metric='haversine')

            Counts_02_23_CLU_BKG = np.empty((len(XGAS_filtered), len(apertures)), dtype=int)
            for j, ap in enumerate(apertures):
                r = deg_to_rad * ap
                Counts_02_23_CLU_BKG[:, j] = Tree_BKG.query_radius(coord_cat_CLU, r=r, count_only=True)
            XGAS_filtered.add_column(Column(name='COUNTS_02_23_CLU_BKG', data=Counts_02_23_CLU_BKG))

    hdul_exp = fits.open(os.path.join(basedir, tile, subdir, 'eSASS', tile+'_024_ExpMap.fits'))
    ExpMap = hdul_exp[0].data
    hdul_Bg3 = fits.open(os.path.join(basedir, tile, subdir, 'eSASS', tile + '_024_Bg3Map.fits'))
    Bg3Map = hdul_Bg3[0].data
    wcs = WCS(hdul_exp[0].header)
    # Convert sky coordinates (RA, DEC) to pixel coordinates
    x_pix, y_pix = wcs.wcs_world2pix(XGAS_filtered['RA'], XGAS_filtered['DEC'], 0)
    x_pix = np.round(x_pix).astype(int)
    y_pix = np.round(y_pix).astype(int)

    # Clip to stay within image bounds
    ny, nx = ExpMap.shape
    x_pix = np.clip(x_pix, 0, nx - 1)
    y_pix = np.clip(y_pix, 0, ny - 1)

    # Extract exposure time at those pixel positions
    Texp = ExpMap[y_pix, x_pix]
    BgMap = Bg3Map[y_pix, x_pix]
    XGAS_filtered['ExpMap_eSASS'] = Texp
    XGAS_filtered['Bg3Map_eSASS'] = BgMap

    #Add eSASS detections
    eSASS_names = ['RA', 'DEC', 'DET_LIKE_0', 'EXT', 'EXT_LIKE', 'ML_CTS_0', 'ML_RATE_0', 'ML_FLUX_0']
    for col in eSASS_names:
        if col in ['RA', 'DEC']:
            XGAS_filtered[col+'_eSASS'] = np.ones(len(XGAS_filtered))*-99.
        else:
            XGAS_filtered[col] = np.ones(len(XGAS_filtered))*-99.

    # Create a mapping from ID to row index in Uniq
    id_to_index = {id_val: i for i, id_val in enumerate(Uniq['ID_uniq'])}

    # For each row in glist, check if id is in Uniq
    for i, gid in enumerate(XGAS_filtered['id'].astype('int32')):
        idx = id_to_index.get(gid)
        if idx is not None:
            for col in eSASS_names:
                # Rename RA and DEC
                target_col = col + '_eSASS' if col in ['RA', 'DEC'] else col
                XGAS_filtered[target_col][i] = Uniq[col][idx]

    eSASS_offset = np.ones(len(XGAS_filtered))*-99.
    seldet = XGAS_filtered['DET_LIKE_0']>0
    coord1 = SkyCoord(XGAS_filtered[seldet]['RA'], XGAS_filtered[seldet]['DEC'], unit='deg')
    coord2 = SkyCoord(XGAS_filtered[seldet]['RA_eSASS'], XGAS_filtered[seldet]['DEC_eSASS'], unit='deg')
    eSASS_offset[seldet] = coord1.separation(coord2).arcsecond
    XGAS_filtered['eSASS_offset'] = eSASS_offset
    XGAS_filtered['eSASS_offset'].units = 'arcsec'
    ipix = hp.ang2pix(nside, XGAS_filtered['RA'], XGAS_filtered['DEC'], nest=True, lonlat=True)
    in_good_region = ~np.isin(ipix, high_BKG_ids)
    XGAS_filtered['in_good_region'] = in_good_region
    clu_all.append(XGAS_filtered)



    # Now do eSASS
    SRV_value = np.array([get_srvmap(e0, e1) for e0, e1 in zip(Uniq['RA'], Uniq['DEC'])])
    sel_area = np.isin(SRV_value.T[0], int(tile))
    Uniq = Uniq[sel_area]
    # unique counterpart of AGN or Star
    PNT_sel =  (Uniq['ID_uniq'] >= 0) & (Uniq['ID_uniq'] < 1e5)
    # unique counterpart of cluster
    EXT_sel = (Uniq['ID_uniq'] >= 1e6)
    # no unique counterpart and linked to agn or star
    PNT2_sel = ((Uniq['ID_uniq'] < 0) & ((Uniq['ID_Any'] >= 0) & (Uniq['ID_Any'] < 1e5)))
    # no unique counterpart and linked to cluster
    EXT2_sel = (Uniq['ID_uniq'] < 0) & (Uniq['ID_Any'] >= 1e6)
    # BKG/FG feature
    # spurious, not linked to anything
    BKG_sel = (Uniq['ID_Any'] < 0)  # & (~sel_BG_feature)

    #Uniq['SRVMAP'] = SRV_value[sel_area]
    Uniq['class_PNT'] = PNT_sel
    Uniq['class_PNT2'] = PNT2_sel
    Uniq['class_EXT'] = EXT_sel
    Uniq['class_EXT2'] = EXT2_sel
    Uniq['class_BKG'] = BKG_sel
    ipix = hp.ang2pix(nside, Uniq['RA'], Uniq['DEC'], nest=True, lonlat=True)
    in_good_region = ~np.isin(ipix, high_BKG_ids)
    Uniq['in_good_region'] = in_good_region

    eSASS_all.append(Uniq)

if len(clu_all)>0:
    clu_all_table_Lext0 = vstack(clu_all)
    clu_all_table_Lext0.write(p_2_clu_matched_Lext0, overwrite=True)

    eSASS_all_table = vstack(eSASS_all)
    eSASS_all_table.write(p_2_eSASS_matched_Lext0, overwrite=True)

    print('Written', p_2_eSASS_matched_Lext0, 'and the basic version of',p_2_clu_matched_Lext0)
    print('Now adding Lx<4r500 and CTR reading the table from the standard run...')

    # --- Read original summary (with Lx<4R500 and CTR columns) ---
    clu_all_orig = Table.read(p_2_clu_matched)  # LC_eRASS_Clusters_Input_<subdir>.fits

    # --- Sanity checks on IDs ---
    # 1. Every Lext0 cluster id must exist in the original summary
    assert np.all(np.isin(clu_all_table_Lext0['id'], clu_all_orig['id'])), \
        "Some Lext0 IDs are not present in the original summary!"

    # --- Prepare columns to copy from the original summary ---
    cols_to_copy = [
        'CLUSTER_LX_soft_RF_4R500c',
        'CTR_R500c_unabs',
        'CTR_4R500c_unabs'
    ]

    missing_cols = [c for c in cols_to_copy if c not in clu_all_orig.colnames]
    if len(missing_cols) > 0:
        raise RuntimeError(f"Original summary missing columns: {missing_cols}")

    # Build a small table with only id + the needed columns
    orig_small = clu_all_orig['id']  # id column
    for c in cols_to_copy:
        orig_small[c] = clu_all_orig[c]

    # --- Join by 'id' to bring those columns into the Lext0 table ---
    clu_all_table_Lext0 = join(
        clu_all_table_Lext0,
        orig_small,
        keys='id',
        join_type='left',
        table_names=['', '_orig'],
        uniq_col_name='{col_name}{table_name}'
    )

    # Now clu_all_Lext0 has the three copied columns with exactly the same values
    # as in the original run, matched by id.

    # Save the Lext0 cluster summary with the reused columns
    clu_all_table_Lext0.write(p_2_clu_matched_Lext0, overwrite=True)

    print('Written', p_2_clu_matched_Lext0, 'and', p_2_eSASS_matched_Lext0, 'using', len(clu_all), 'tiles')
else:
    print('No tiles found to process!')