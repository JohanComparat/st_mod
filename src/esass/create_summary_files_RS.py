#!/usr/bin/env python3
import sys
import os
import glob
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed

from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table, join, vstack, Column
from astropy.coordinates import SkyCoord
from sklearn.neighbors import BallTree
from scipy.interpolate import RegularGridInterpolator
from tqdm import tqdm
from astropy.cosmology import FlatLambdaCDM
import healpy as hp
import xspec
import astropy.units as u

# ----------------------------------------------------------------------
# Cosmology / constants
# ----------------------------------------------------------------------
cosmoUCHUU = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
h = 0.6774
cosmo = cosmoUCHUU
deg_to_rad = np.pi / 180.0

"""
Creates summary files in /home/idies/workspace/erosim/Uchuu/LCerass/SummaryFiles
Example:
    python create_summary_files_RS.py GE_e4_merge_AGNseed001_SimBKG_CLUseed001

R. Seppi (31.03.2026)
Speed-up + fixed ID handling.
"""

# ----------------------------------------------------------------------
# Global paths
# ----------------------------------------------------------------------
basedir = '/home/idies/workspace/erosim/Uchuu/LCerass'
subdir = sys.argv[1]

print('Directory:', subdir)

print('Reading in p2cats...')
p_2_match_cats = np.sort(glob.glob(os.path.join(basedir, '??????', subdir, 'Sc1*Uniq*fits')))
p_2_match_cats_Lext0 = np.sort(glob.glob(os.path.join(basedir, '??????', subdir, 'Sc_Lext0*Uniq*fits')))

p_2_clu_matched = os.path.join(basedir, 'SummaryFiles', f'LC_eRASS_Clusters_Input_{subdir}.fits')
p_2_eSASS_matched = os.path.join(basedir, 'SummaryFiles', f'LC_eRASS_eSASS_Output_{subdir}.fits')

p_2_clu_matched_Lext0 = os.path.join(basedir, 'SummaryFiles', f'LC_eRASS_Clusters_Input_{subdir}_Lext0.fits')
p_2_eSASS_matched_Lext0 = os.path.join(basedir, 'SummaryFiles', f'LC_eRASS_eSASS_Output_{subdir}_Lext0.fits')

coolfunc = np.loadtxt('/home/idies/workspace/erosim/software/st_mod_data/data/models/model_GAS/coolfunc.dat')
xgrid_ext = np.loadtxt('/home/idies/workspace/erosim/software/st_mod_data/data/models/model_GAS/radial_binning.txt')
profiles = Table.read(
    '/home/idies/workspace/erosim/software/st_mod_data/data/models/model_GAS/profiles_010z015_1e14M2e14.fits')

rsp_path = "/home/idies/workspace/erosim/software/st_mod_data/data/onaxis_tm0_rmf_2023-01-17.fits"
arf_path = "/home/idies/workspace/erosim/software/st_mod_data/data/survey_tm0_arf_filter_2023-01-17.fits"
p2ecf_grid = "/home/idies/workspace/erosim/software/st_mod_data/data/ecf_grid.npz"


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------
def calc_lx(prof, kt, m5, z, fraction=1):
    Mpc = 3.0856776e+24
    msun = 1.98892e33
    ez2 = cosmo.efunc(z) ** 2
    rhoc = cosmo.critical_density(z).value
    r500 = np.power(m5 * msun / 4. * 3. / np.pi / 500. / rhoc, 1. / 3.)
    resfact = np.sqrt(kt / 10.0) * np.power(ez2, 3. / 2.)
    prof_em = prof * resfact
    tlambda = np.interp(kt, coolfunc[:, 0], coolfunc[:, 1])
    dx = np.empty(len(xgrid_ext))
    dx[0] = xgrid_ext[0]
    dx[1:] = (np.roll(xgrid_ext, -1) - xgrid_ext)[:-1]
    lxcum = np.cumsum(prof_em * xgrid_ext * r500 ** 2 * 2. * np.pi * tlambda * Mpc * dx)
    return np.interp(fraction, xgrid_ext, lxcum)


def ctr_factor_xspec(z, Z_solar, kT_keV,
                     rsp_path, arf_path,
                     band_rest_keV=(0.5, 2.0),
                     abund='aspl'):
    import os

    E1_rest, E2_rest = band_rest_keV
    xcm_file = 'temp_Lx_rate.xcm'
    log_file = 'temp_Lx_rate.log'

    with open(xcm_file, 'w') as f:
        print(f"abund {abund}", file=f)
        print("model apec", file=f)
        print(f"{kT_keV:.4f}", file=f)
        print(f"{Z_solar:.4f}", file=f)
        print(f"{z:.4f}", file=f)
        print("1.0", file=f)

        print("fakeit none", file=f)
        print(rsp_path, file=f)
        print(arf_path, file=f)
        print("", file=f)
        print("", file=f)
        print("temp.pi", file=f)
        print("1.0", file=f)
        print("", file=f)

        print(f"log {log_file}", file=f)
        print(f"lumin {E1_rest:.2f} {E2_rest:.2f} {z:.4f}", file=f)
        print("show rate", file=f)
        print("log none", file=f)
        print("quit", file=f)
        print("y", file=f)

    os.system(f"xspec - {xcm_file} > /dev/null 2>&1")
    os.system('rm -f temp.pi')

    L_at_norm_1 = None
    rate_at_norm_1 = None
    try:
        with open(log_file, "r") as f:
            for line in f:
                if "Model Luminosity" in line:
                    L_at_norm_1 = float(line.split()[2])
                if "Model predicted rate:" in line:
                    rate_at_norm_1 = float(line.split()[4])
    finally:
        os.system(f"rm -f {xcm_file} {log_file} temp.pi")

    if (L_at_norm_1 is None) or (L_at_norm_1 <= 0.0):
        raise RuntimeError(f"Bad L_at_norm_1: {L_at_norm_1}")
    if (rate_at_norm_1 is None) or (rate_at_norm_1 < 0.0):
        raise RuntimeError(f"Bad rate_at_norm_1: {rate_at_norm_1}")

    return rate_at_norm_1 / L_at_norm_1


def build_ecf_grid(kT_grid, z_grid, rsp_path, arf_path,
                   band_rest_keV=(0.5, 2.0), Z_solar=0.3,
                   abund='aspl', out_file='ecf_grid.npz'):
    print('Building ECF grid in', out_file)
    ecf = np.zeros((len(kT_grid), len(z_grid)), dtype=float)

    for i, kT in enumerate(tqdm(kT_grid)):
        for j, z in enumerate(z_grid):
            ecf[i, j] = ctr_factor_xspec(
                z=z,
                Z_solar=Z_solar,
                kT_keV=kT,
                rsp_path=rsp_path,
                arf_path=arf_path,
                band_rest_keV=band_rest_keV,
                abund=abund
            )

    np.savez(out_file, kT_grid=kT_grid, z_grid=z_grid, ecf=ecf)
    print(f"Saved ECF grid to {out_file}")
    return ecf


# ----------------------------------------------------------------------
# Build/load ECF grid once in parent process
# ----------------------------------------------------------------------
kT_grid = np.arange(0.1, 13.1, 0.2)
z_grid = np.arange(0.01, 1.51, 0.02)

if not os.path.isfile(p2ecf_grid):
    build_ecf_grid(kT_grid, z_grid, rsp_path, arf_path, out_file=p2ecf_grid)

# ----------------------------------------------------------------------
# Read sky map and precompute bounds once
# ----------------------------------------------------------------------
try:
    sky_map_hdu = Table.read(os.path.join(os.environ['GIT_STMOD_DATA'],
                                          'data/models/eROSITA', 'SKYMAPS.fits'))
except Exception:
    sky_map_hdu = Table.read('/home/idies/workspace/erosim/software/st_mod_data/data/models/eROSITA/SKYMAPS.fits')

tile_bounds = {}
for row in sky_map_hdu:
    tile_id = int(row['SRVMAP'])
    tile_bounds[tile_id] = (
        row['RA_MIN'],
        row['RA_MAX'],
        row['DE_MIN'],
        row['DE_MAX']
    )


def get_inner_tile_bounds(tile, border_arcmin=20.0):
    ra_min, ra_max, dec_min, dec_max = tile_bounds[int(tile)]
    delta = border_arcmin / 60.0
    return (
        ra_min + delta,
        ra_max - delta,
        dec_min + delta,
        dec_max - delta
    )


# ----------------------------------------------------------------------
# Background map info
# ----------------------------------------------------------------------
p2bkg_simput = os.path.join(os.environ['HOME'], 'workspace/erosim/simput/bkg_erosita_simput_full_sky', 'catalogue.fits')
bkg_simput = Table.read(p2bkg_simput, memmap=True)
sel_BG_feature = bkg_simput['FLUX'] > 2e-11
high_BKG_ids = bkg_simput['SRC_ID'][sel_BG_feature]
nside = 128


# ----------------------------------------------------------------------
# Tile processing
# ----------------------------------------------------------------------
def process_tile(tile, mode='standard'):
    """
    Process one tile.
    mode = 'standard' or 'lext0'
    Returns:
        clu_tile_table, esass_tile_table
    """
    try:
        tile_int = int(tile)

        if mode == 'standard':
            uniq_path = os.path.join(basedir, tile, subdir, f'Sc1_{tile}_IDMatch_Uniq_Tot2.0.fits')
            any_path = os.path.join(basedir, tile, subdir, f'Sc1_{tile}_IDMatch_Any_Tot2.0.fits')
        else:
            uniq_path = os.path.join(basedir, tile, subdir, f'Sc_Lext0_{tile}_IDMatch_Uniq_Tot2.0.fits')
            any_path = os.path.join(basedir, tile, subdir, f'Sc_Lext0_{tile}_IDMatch_Any_Tot2.0.fits')

        if (not os.path.isfile(uniq_path)) or (not os.path.isfile(any_path)):
            return None, None

        xgas_path = os.path.join(basedir, tile, 'XGAS.fits')
        clu_simput_path = os.path.join(basedir, tile, 'Xgas_bHS0.8_simput_final.fits')
        if (not os.path.isfile(xgas_path)) or (not os.path.isfile(clu_simput_path)):
            return None, None

        XGAS = Table.read(xgas_path, memmap=True)

        # --- cluster table restricted to cluster simput IDs
        clu = Table.read(clu_simput_path, memmap=True)
        mask = np.isin(np.asarray(XGAS['id'], dtype=np.int64), np.asarray(clu['SRC_ID'], dtype=np.int64))
        XGAS_filtered = XGAS[mask]

        if len(XGAS_filtered) == 0:
            return None, None

        Any = Table.read(any_path, memmap=True)
        Uniq = Table.read(uniq_path, memmap=True)

        # --- fixed ID semantics
        Uniq['ID_simput'] = np.asarray(Uniq['ID_simput'], dtype=np.int64)
        Any['ID_simput'] = np.asarray(Any['ID_simput'], dtype=np.int64)

        if 'ID_uniq' not in Uniq.colnames:
            Uniq.rename_column('ID_simput', 'ID_uniq')
        if 'ID_simput32' in Uniq.colnames and 'ID_uniq32' not in Uniq.colnames:
            Uniq.rename_column('ID_simput32', 'ID_uniq32')

        # safer Any -> Uniq mapping by ID_cat, not row order
        any_map = {int(row['ID_cat']): int(row['ID_simput']) for row in Any}
        Uniq['ID_Any'] = np.array([any_map.get(int(x), -99) for x in Uniq['ID_cat']], dtype=np.int64)

        if 'ID_simput32' in Any.colnames:
            any32_map = {int(row['ID_cat']): int(row['ID_simput32']) for row in Any}
            Uniq['ID_Any32'] = np.array([any32_map.get(int(x), -99) for x in Uniq['ID_cat']], dtype=np.int32)

        # --- restrict Uniq to the tile area using precomputed box, not get_srvmap per source
        if tile_int not in tile_bounds:
            return None, None

        ra_min, ra_max, dec_min, dec_max = tile_bounds[tile_int]
        sel_area = (
                (Uniq['RA'] >= ra_min) &
                (Uniq['RA'] < ra_max) &
                (Uniq['DEC'] >= dec_min) &
                (Uniq['DEC'] < dec_max)
        )
        Uniq = Uniq[sel_area]

        # add flag within border
        ra_min_in, ra_max_in, dec_min_in, dec_max_in = get_inner_tile_bounds(tile, border_arcmin=30.0)
        within_border = (
                (Uniq['RA'] >= ra_min_in) &
                (Uniq['RA'] <= ra_max_in) &
                (Uniq['DEC'] >= dec_min_in) &
                (Uniq['DEC'] <= dec_max_in)
        )

        Uniq['within_border'] = within_border

        R500c_arcmin = XGAS_filtered['R500c'] / cosmo.kpc_proper_per_arcmin(z=XGAS_filtered['redshift_S'])
        XGAS_filtered['R500c_arcmin'] = R500c_arcmin.value

        # --- photon statistics
        apertures = np.array([0.1, 0.3, 0.5, 1.0, 2.0])

        p2evt_clu = os.path.join(basedir, tile, subdir, f'simCLUevt_{tile}.fits')
        p2evt_agn = os.path.join(basedir, tile, subdir, f'simAGNevt_{tile}.fits')
        p2evt_bkg = os.path.join(basedir, tile, subdir, f'simBKGevt_{tile}.fits')
        coord_cat_CLU = deg_to_rad * np.transpose([XGAS_filtered['DEC'], XGAS_filtered['RA']])

        if os.path.isfile(p2evt_clu):
            clu_evt = Table.read(p2evt_clu, memmap=True)
            sel_evt = (clu_evt['SIGNAL'] > 0.2) & (clu_evt['SIGNAL'] < 2.3)
            clu_evt = clu_evt[sel_evt]

            coord_evt_CLU = deg_to_rad * np.transpose([clu_evt['DEC'], clu_evt['RA']])
            Tree_CLU = BallTree(coord_evt_CLU, metric='haversine')

            counts = np.empty((len(XGAS_filtered), len(apertures)), dtype=int)
            for j, ap in enumerate(apertures):
                r = deg_to_rad * ap * XGAS_filtered['R500c_arcmin'] / 60.0
                counts[:, j] = Tree_CLU.query_radius(coord_cat_CLU, r=r, count_only=True)
            XGAS_filtered.add_column(Column(name='COUNTS_02_23_CLU_CLU', data=counts))

        if os.path.isfile(p2evt_agn):
            agn_evt = Table.read(p2evt_agn, memmap=True)
            sel_evt = (agn_evt['SIGNAL'] > 0.2) & (agn_evt['SIGNAL'] < 2.3)
            agn_evt = agn_evt[sel_evt]

            coord_evt_AGN = deg_to_rad * np.transpose([agn_evt['DEC'], agn_evt['RA']])
            Tree_AGN = BallTree(coord_evt_AGN, metric='haversine')

            counts = np.empty((len(XGAS_filtered), len(apertures)), dtype=int)
            for j, ap in enumerate(apertures):
                r = deg_to_rad * ap
                counts[:, j] = Tree_AGN.query_radius(coord_cat_CLU, r=r, count_only=True)
            XGAS_filtered.add_column(Column(name='COUNTS_02_23_CLU_AGN', data=counts))

        if os.path.isfile(p2evt_bkg):
            bkg_evt = Table.read(p2evt_bkg, memmap=True)
            sel_evt = (bkg_evt['SIGNAL'] > 0.2) & (bkg_evt['SIGNAL'] < 2.3)
            bkg_evt = bkg_evt[sel_evt]

            coord_evt_BKG = deg_to_rad * np.transpose([bkg_evt['DEC'], bkg_evt['RA']])
            Tree_BKG = BallTree(coord_evt_BKG, metric='haversine')

            counts = np.empty((len(XGAS_filtered), len(apertures)), dtype=int)
            for j, ap in enumerate(apertures):
                r = deg_to_rad * ap
                counts[:, j] = Tree_BKG.query_radius(coord_cat_CLU, r=r, count_only=True)
            XGAS_filtered.add_column(Column(name='COUNTS_02_23_CLU_BKG', data=counts))

        # --- exposure / background maps
        exp_path = os.path.join(basedir, tile, subdir, 'eSASS', f'{tile}_024_ExpMap.fits')
        bg_path = os.path.join(basedir, tile, subdir, 'eSASS', f'{tile}_024_Bg3Map.fits')
        if os.path.isfile(exp_path) and os.path.isfile(bg_path):
            with fits.open(exp_path) as hdul_exp:
                ExpMap = hdul_exp[0].data
                wcs = WCS(hdul_exp[0].header)

            with fits.open(bg_path) as hdul_bg:
                Bg3Map = hdul_bg[0].data

            x_pix, y_pix = wcs.wcs_world2pix(XGAS_filtered['RA'], XGAS_filtered['DEC'], 0)
            x_pix = np.round(x_pix).astype(int)
            y_pix = np.round(y_pix).astype(int)

            ny, nx = ExpMap.shape
            x_pix = np.clip(x_pix, 0, nx - 1)
            y_pix = np.clip(y_pix, 0, ny - 1)

            XGAS_filtered['ExpMap_eSASS'] = ExpMap[y_pix, x_pix]
            XGAS_filtered['Bg3Map_eSASS'] = Bg3Map[y_pix, x_pix]

        # --- fast XGAS <- Uniq join instead of Python loop
        eSASS_names = ['ID_cat', 'RA', 'DEC', 'RADEC_ERR', 'DET_LIKE_0', 'EXT',
                       'EXT_LIKE', 'DIST_NN', 'ML_CTS_0', 'ML_CTS_ERR_0',
                       'ML_RATE_0', 'ML_RATE_ERR_0', 'ML_FLUX_0', 'ML_FLUX_ERR_0',
                       'ML_BKG_0', 'ML_EXP_1']

        uniq_small = Table()
        uniq_small['id'] = np.asarray(Uniq['ID_uniq'], dtype=np.int64)

        for col in eSASS_names:
            if col in ['RA', 'DEC']:
                uniq_small[col + '_eSASS'] = Uniq[col]
            else:
                uniq_small[col] = Uniq[col]

        XGAS_filtered['id'] = np.asarray(XGAS_filtered['id'], dtype=np.int64)
        XGAS_filtered = join(
            XGAS_filtered,
            uniq_small,
            keys='id',
            join_type='left',
            table_names=['', '_uniq'],
            uniq_col_name='{col_name}{table_name}'
        )

        within_border_clu = (
                (XGAS_filtered['RA'] >= ra_min_in) &
                (XGAS_filtered['RA'] <= ra_max_in) &
                (XGAS_filtered['DEC'] >= dec_min_in) &
                (XGAS_filtered['DEC'] <= dec_max_in)
        )

        XGAS_filtered['within_border'] = within_border_clu

        # fill missing values from join, including masked columns
        for col in eSASS_names:
            target_col = col + '_eSASS' if col in ['RA', 'DEC'] else col

            if target_col not in XGAS_filtered.colnames:
                XGAS_filtered[target_col] = np.full(len(XGAS_filtered), -99.0)
                continue

            coldata = XGAS_filtered[target_col]

            if hasattr(coldata, 'filled'):
                # MaskedColumn -> replace masked entries explicitly
                filled = coldata.filled(-99)
                XGAS_filtered[target_col] = filled
            else:
                arr = np.asarray(coldata)
                if np.issubdtype(arr.dtype, np.floating):
                    bad = ~np.isfinite(arr)
                    if np.any(bad):
                        arr = arr.copy()
                        arr[bad] = -99.0
                        XGAS_filtered[target_col] = arr

        # --- offsets
        eSASS_offset = np.full(len(XGAS_filtered), -99.0)
        eSASS_offset_kpc = np.full(len(XGAS_filtered), -99.0)

        matched = np.asarray(XGAS_filtered['ID_cat']) >= 0
        if np.any(matched):
            coord1 = SkyCoord(
                XGAS_filtered['RA'][matched],
                XGAS_filtered['DEC'][matched],
                unit='deg'
            )
            coord2 = SkyCoord(
                XGAS_filtered['RA_eSASS'][matched],
                XGAS_filtered['DEC_eSASS'][matched],
                unit='deg'
            )

            sep_angle = coord1.separation(coord2)
            eSASS_offset[matched] = sep_angle.arcsecond

            DA = cosmo.angular_diameter_distance(XGAS_filtered['redshift_S'][matched])
            offset_kpc = (sep_angle.to(u.rad).value * DA).to(u.kpc)
            eSASS_offset_kpc[matched] = offset_kpc.value

        XGAS_filtered['eSASS_offset'] = eSASS_offset
        XGAS_filtered['eSASS_offset'].unit = 'arcsec'
        XGAS_filtered['eSASS_offset_kpc'] = eSASS_offset_kpc
        XGAS_filtered['eSASS_offset_kpc'].unit = 'kpc'

        # add N split sources
        # Primary Any counts by cluster ID
        any_primary = np.asarray(Any['ID_simput'], dtype=np.int64)
        valid_primary = any_primary >= 1e6
        u_primary, c_primary = np.unique(any_primary[valid_primary], return_counts=True)
        n_primary_map = dict(zip(u_primary, c_primary))

        # From Uniq, build:
        # cluster_id -> final matched row info
        uniq_ids = np.asarray(Uniq['ID_uniq'], dtype=np.int64)
        uniq_any = np.asarray(Uniq['ID_Any'], dtype=np.int64)

        # keep only matched cluster rows
        sel_uniq_cluster = uniq_ids >= 1e6
        uniq_ids = uniq_ids[sel_uniq_cluster]
        uniq_any = uniq_any[sel_uniq_cluster]

        uniq_any_map = dict(zip(uniq_ids, uniq_any))

        cluster_ids = np.asarray(XGAS_filtered['id'], dtype=np.int64)
        N_eSASS_det = np.zeros(len(XGAS_filtered), dtype=np.int32)

        for i, cid in enumerate(cluster_ids):
            n_det = n_primary_map.get(int(cid), 0)

            # if this cluster has a final unique match that came from promoted secondary,
            # add one because that source is not counted in Any primary == cid
            if cid in uniq_any_map:
                if int(uniq_any_map[cid]) != int(cid):
                    n_det += 1

            N_eSASS_det[i] = n_det

        XGAS_filtered['N_eSASS_det'] = N_eSASS_det

        ipix = hp.ang2pix(nside, XGAS_filtered['RA'], XGAS_filtered['DEC'], nest=True, lonlat=True)
        XGAS_filtered['in_good_region'] = ~np.isin(ipix, high_BKG_ids)

        # --- classify eSASS detections
        PNT_sel = (Uniq['ID_uniq'] >= 0) & (Uniq['ID_uniq'] < 1e5)
        EXT_sel = (Uniq['ID_uniq'] >= 1e6)
        PNT2_sel = ((Uniq['ID_uniq'] < 0) & ((Uniq['ID_Any'] >= 0) & (Uniq['ID_Any'] < 1e5)))
        EXT2_sel = ((Uniq['ID_uniq'] < 0) & (Uniq['ID_Any'] >= 1e6))
        BKG_sel = (Uniq['ID_Any'] < 0)

        Uniq['class_PNT'] = PNT_sel
        Uniq['class_PNT2'] = PNT2_sel
        Uniq['class_EXT'] = EXT_sel
        Uniq['class_EXT2'] = EXT2_sel
        Uniq['class_BKG'] = BKG_sel

        ipix = hp.ang2pix(nside, Uniq['RA'], Uniq['DEC'], nest=True, lonlat=True)
        Uniq['in_good_region'] = ~np.isin(ipix, high_BKG_ids)

        return XGAS_filtered, Uniq

    except Exception as e:
        print(f"[{mode}] tile {tile} failed: {e}")
        return None, None


# ----------------------------------------------------------------------
# Parallel runner
# ----------------------------------------------------------------------
def run_parallel(tiles, mode='standard', max_workers=None):
    clu_all = []
    esass_all = []

    if max_workers is None:
        ncpu = os.cpu_count() or 4
        max_workers = max(1, min(16, ncpu - 1))

    print(f'Running mode={mode} with {max_workers} workers over {len(tiles)} tiles')

    with ProcessPoolExecutor(max_workers=max_workers) as exe:
        futures = {exe.submit(process_tile, tile, mode): tile for tile in tiles}

        for fut in tqdm(as_completed(futures), total=len(futures), desc=f'Processing {mode}'):
            tile = futures[fut]
            clu_tile, esass_tile = fut.result()

            if clu_tile is not None and len(clu_tile) > 0:
                clu_all.append(clu_tile)
            if esass_tile is not None and len(esass_tile) > 0:
                esass_all.append(esass_tile)

    return clu_all, esass_all


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
if __name__ == '__main__':
    # standard
    tiles_standard = [p.split('/')[-3] for p in p_2_match_cats]
    clu_all, eSASS_all = run_parallel(tiles_standard, mode='standard')

    if len(clu_all) > 0:
        clu_all_table = vstack(clu_all)
        clu_all_table.write(p_2_clu_matched, overwrite=True)

        eSASS_all_table = vstack(eSASS_all)
        eSASS_all_table.write(p_2_eSASS_matched, overwrite=True)

        print('Written', p_2_eSASS_matched, 'and the basic version of', p_2_clu_matched)
        print('Now adding Lx<4r500...')

        CLUSTER_LX_soft_RF_4R500c = np.zeros(len(clu_all_table))
        for jj, cl in enumerate(tqdm(clu_all_table)):
            prof = profiles[cl['idx_profile']]['profiles']
            L500_prof = calc_lx(prof, cl['CLUSTER_kT'], cl['M500c'], cl['redshift_S'], fraction=1)
            resc_fact = np.power(10, cl['CLUSTER_LX_soft_RF_R500c']) / L500_prof
            CLUSTER_LX_soft_RF_4R500c[jj] = np.log10(
                resc_fact * calc_lx(prof, cl['CLUSTER_kT'], cl['M500c'], cl['redshift_S'], fraction=4)
            )
        clu_all_table['CLUSTER_LX_soft_RF_4R500c'] = CLUSTER_LX_soft_RF_4R500c

        print('Now adding the unabsorbed count rate...')
        ecf_data = np.load(p2ecf_grid)
        kT_grid = ecf_data['kT_grid']
        z_grid = ecf_data['z_grid']
        ecf_grid = ecf_data['ecf']

        ecf_interp = RegularGridInterpolator(
            (kT_grid, z_grid),
            ecf_grid,
            bounds_error=False,
            fill_value=None
        )

        points = np.column_stack([
            clu_all_table['CLUSTER_kT'].astype(float),
            clu_all_table['redshift_S'].astype(float)
        ])
        ecf_arr = ecf_interp(points)

        Lx_R500c = 10 ** clu_all_table['CLUSTER_LX_soft_RF_R500c']
        Lx_4R500c = 10 ** clu_all_table['CLUSTER_LX_soft_RF_4R500c']

        clu_all_table['CTR_R500c_unabs'] = Lx_R500c * ecf_arr
        clu_all_table['CTR_4R500c_unabs'] = Lx_4R500c * ecf_arr
        clu_all_table.write(p_2_clu_matched, overwrite=True)

        print('Written', p_2_clu_matched, 'and', p_2_eSASS_matched, 'using', len(clu_all), 'tiles')
    else:
        print('No tiles found to process for standard run!')
    sys.exit()
    # Lext0
    print('Now doing summary catalogues for the extlikemin=0 run...')
    tiles_lext0 = [p.split('/')[-3] for p in p_2_match_cats_Lext0]
    clu_all, eSASS_all = run_parallel(tiles_lext0, mode='lext0')

    if len(clu_all) > 0:
        clu_all_table_Lext0 = vstack(clu_all)
        clu_all_table_Lext0.write(p_2_clu_matched_Lext0, overwrite=True)

        eSASS_all_table = vstack(eSASS_all)
        eSASS_all_table.write(p_2_eSASS_matched_Lext0, overwrite=True)

        print('Written', p_2_eSASS_matched_Lext0, 'and the basic version of', p_2_clu_matched_Lext0)
        print('Now adding Lx<4r500 and CTR reading the table from the standard run...')

        clu_all_orig = Table.read(p_2_clu_matched)

        assert np.all(np.isin(clu_all_table_Lext0['id'], clu_all_orig['id'])), \
            "Some Lext0 IDs are not present in the original summary!"

        cols_to_copy = [
            'CLUSTER_LX_soft_RF_4R500c',
            'CTR_R500c_unabs',
            'CTR_4R500c_unabs'
        ]

        missing_cols = [c for c in cols_to_copy if c not in clu_all_orig.colnames]
        if len(missing_cols) > 0:
            raise RuntimeError(f"Original summary missing columns: {missing_cols}")

        orig_small = clu_all_orig['id']
        for c in cols_to_copy:
            orig_small[c] = clu_all_orig[c]

        clu_all_table_Lext0 = join(
            clu_all_table_Lext0,
            orig_small,
            keys='id',
            join_type='left',
            table_names=['', '_orig'],
            uniq_col_name='{col_name}{table_name}'
        )

        clu_all_table_Lext0.write(p_2_clu_matched_Lext0, overwrite=True)
        print('Written', p_2_clu_matched_Lext0, 'and', p_2_eSASS_matched_Lext0, 'using', len(clu_all), 'tiles')
    else:
        print('No tiles found to process for Lext0 run!')