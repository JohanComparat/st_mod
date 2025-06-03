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

'''
Creates summmary files in /home/idies/workspace/erosim/Uchuu/LCerass/SummaryFiles
takes the desired experiment directory as input, for example:

python create_summary_files_RS.py GE_e4_merge_AGNseed001_SimBKG_CLUseed001

R. Seppi (19.05.2025)
'''

deg_to_rad = np.pi/180.
cosmo = FlatLambdaCDM(H0=67.77, Om0=0.29)

basedir = '/home/idies/workspace/erosim/Uchuu/LCerass'

subdir = sys.argv[1]

print('Directory:', subdir)

p_2_match_cats = np.sort(glob.glob(os.path.join(basedir, '??????', subdir, 'Sc1*Uniq*fits')))
p_2_clu_matched = os.path.join(basedir, 'SummaryFiles', 'LC_eRASS_Clusters_Input_'+subdir+'.fits')
p_2_eSASS_matched = os.path.join(basedir, 'SummaryFiles', 'LC_eRASS_eSASS_Output_'+subdir+'.fits')


#prepare to select unique area in esass catalogues
#sky_map_hdu = Table.read(os.path.join(os.environ['GIT_STMOD_DATA'],
#                                      'data/models/eROSITA', 'SKYMAPS.fits'))
sky_map_hdu = Table.read('SKYMAPS.fits')
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

    print('Written', p_2_clu_matched, 'and', p_2_eSASS_matched, 'using', len(clu_all), 'tiles')
else:
    print('No tiles found to process!')