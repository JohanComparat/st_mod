#!/usr/bin/env python3
"""
Photon based input-output matches for eROSITA simulations.

Example:
    python photon_matching_RS.py GE_e4_merge_AGNseed001_SimBKG_CLUseed001 121048

R. Seppi (27.03.2026)
Patched to fix int32 SRC_ID collisions a posteriori using positional unwrapping.
"""

import sys
import os
import numpy as np
from astropy.table import Table, vstack
from astropy.io import fits
from sklearn.neighbors import BallTree
from scipy import stats


MinTotalCts = 2.0
poiss_ppf = 0.8
MinAperCts = 1   # There has to be at least one to know the SimputID

basedir = '/home/idies/workspace/erosim/Uchuu/LCerass'


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------
def angsep_arcsec(ra1, dec1, ra2, dec2):
    """
    Approximate angular separation in arcsec.
    Good enough here because we only use it to choose among collided IDs.
    Inputs can be scalars or arrays.
    """
    dra = (ra1 - ra2) * np.cos(np.deg2rad(0.5 * (dec1 + dec2)))
    ddec = dec1 - dec2
    return 3600.0 * np.sqrt(dra**2 + ddec**2)


def build_wrap_lookup(simput_id64):
    """
    Build wrapped int32 IDs and map wrapped id -> list of indices in simput arrays.
    """
    simput_id32 = simput_id64.astype(np.int32)
    wrap_to_idx = {}
    for i, wid in enumerate(simput_id32):
        wrap_to_idx.setdefault(int(wid), []).append(i)
    return simput_id32, wrap_to_idx


def unwrap_srcid_by_position(wrapped_id, ra_det, dec_det,
                             wrap_to_idx, simput_id64, simput_ra, simput_dec):
    """
    Resolve one wrapped int32 source id to the most likely original int64 source id,
    using the detection position.
    Returns:
        true_id64, n_candidates_for_this_wrapped_id, sep_arcsec_to_chosen_candidate
    """
    cand_idx = wrap_to_idx.get(int(np.int32(wrapped_id)), [])
    if len(cand_idx) == 0:
        return -99, 0, np.nan

    if len(cand_idx) == 1:
        j = cand_idx[0]
        sep = angsep_arcsec(ra_det, dec_det, simput_ra[j], simput_dec[j])
        return int(simput_id64[j]), 1, float(sep)

    cand_ra = simput_ra[cand_idx]
    cand_dec = simput_dec[cand_idx]
    sep = angsep_arcsec(ra_det, dec_det, cand_ra, cand_dec)
    k = np.argmin(sep)
    j = cand_idx[k]
    return int(simput_id64[j]), len(cand_idx), float(sep[k])


def matchlines(list1, list2):
    """
    Return indices in list2 corresponding to values in list1.
    If not found, return -1.
    """
    mapper = {r: i for i, r in enumerate(list2)}
    return np.array([mapper.get(x, -1) for x in list1], dtype=int)


def read_table_if_exists(path, memmap=False):
    if os.path.isfile(path):
        return Table.read(path, memmap=memmap)
    return None


def ensure_1d_srcid(col):
    arr = np.asarray(col)
    if arr.ndim == 1:
        return arr
    if arr.ndim == 2:
        return arr[:, 0]
    raise ValueError(f"SRC_ID has unexpected shape {arr.shape}")


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
subdir = sys.argv[1]
tile = sys.argv[2]

print('Directory:', subdir)
print('Tile:', tile)

workdir = os.path.join(basedir, tile, subdir)
if not os.path.isdir(workdir):
    print(workdir, 'does not exist!')
    sys.exit(1)

# Paths
p_2_esass       = os.path.join(workdir, 'eSASS', tile + '_024_Sc1Cat.fits')
p_2_esass_Lext0 = os.path.join(workdir, 'eSASS', tile + '_024_ScCat_Lext0.fits')
p_2_evt_clu     = os.path.join(workdir, 'simCLUevt_' + tile + '.fits')
p_2_evt_agn     = os.path.join(workdir, 'simAGNevt_' + tile + '.fits')
p_2_evt_sta     = os.path.join(workdir, 'simSTAevt_' + tile + '.fits')
p_2_evt_bkg     = os.path.join(workdir, 'simBKGevt_' + tile + '.fits')

p_2_clu         = os.path.join(basedir, tile, 'Xgas_bHS0.8_simput_final.fits')
p_2_agn         = os.path.join(basedir, tile, 'AGN_list_sigma_0.8_fsat_8.0_simput.fits')
p_2_star        = os.path.join(basedir, tile, 'SIMPUT_STA.fits')

# ----------------------------------------------------------------------
# Read simput source tables, KEEP ORIGINAL int64 IDs
# ----------------------------------------------------------------------
Simput_IDs = []
Simput_RA = []
Simput_DEC = []
Simput_FLUX = []

clu = read_table_if_exists(p_2_clu, memmap=True)
if clu is not None:
    Simput_IDs.append(np.asarray(clu['SRC_ID'], dtype=np.int64))
    Simput_RA.append(np.asarray(clu['RA'], dtype=float))
    Simput_DEC.append(np.asarray(clu['DEC'], dtype=float))
    Simput_FLUX.append(np.asarray(clu['FLUX'], dtype=float))

agn = read_table_if_exists(p_2_agn, memmap=True)
if agn is not None:
    Simput_IDs.append(np.asarray(agn['SRC_ID'], dtype=np.int64))
    Simput_RA.append(np.asarray(agn['RA'], dtype=float))
    Simput_DEC.append(np.asarray(agn['DEC'], dtype=float))
    Simput_FLUX.append(np.asarray(agn['FLUX'], dtype=float))

star = read_table_if_exists(p_2_star, memmap=True)
if star is not None:
    Simput_IDs.append(np.asarray(star['SRC_ID'], dtype=np.int64))
    Simput_RA.append(np.asarray(star['RA'], dtype=float))
    Simput_DEC.append(np.asarray(star['DEC'], dtype=float))
    Simput_FLUX.append(np.asarray(star['FLUX'], dtype=float))

if len(Simput_IDs) == 0:
    raise RuntimeError("No simput source catalogues found.")

Simput_SrcID64 = np.hstack(Simput_IDs).astype(np.int64)
Simput_RA = np.hstack(Simput_RA).astype(float)
Simput_DEC = np.hstack(Simput_DEC).astype(float)
Simput_FLUX = np.hstack(Simput_FLUX).astype(float)

Simput_SrcID32, wrap_to_simput_idx = build_wrap_lookup(Simput_SrcID64)

# ----------------------------------------------------------------------
# Read event files
# ----------------------------------------------------------------------
Evt_Src_list = []

Evt_Clu = read_table_if_exists(p_2_evt_clu)
if Evt_Clu is not None:
    Evt_Src_list.append(Evt_Clu)

Evt_Agn = read_table_if_exists(p_2_evt_agn)
if Evt_Agn is not None:
    Evt_Src_list.append(Evt_Agn)

Evt_Sta = read_table_if_exists(p_2_evt_sta)
if Evt_Sta is not None:
    Evt_Src_list.append(Evt_Sta)

if Evt_Src_list:
    Evt_Src = vstack(Evt_Src_list)
    sel_src = (Evt_Src['SIGNAL'] > 0.2) & (Evt_Src['SIGNAL'] < 2.3)
    Evt_Src = Evt_Src[sel_src]
else:
    Evt_Src = None

Evt_Bkg = read_table_if_exists(p_2_evt_bkg)
if Evt_Bkg is not None:
    sel_bkg = (Evt_Bkg['SIGNAL'] > 0.2) & (Evt_Bkg['SIGNAL'] < 2.3)
    Evt_Bkg = Evt_Bkg[sel_bkg]

if Evt_Src is None and Evt_Bkg is None:
    raise RuntimeError("No source or background event files found.")

# Combine all events and assign SRC_ID = -1 to background
if Evt_Src is not None and Evt_Bkg is not None:
    src_ids_1d = ensure_1d_srcid(Evt_Src['SRC_ID']).astype(np.int32)
    bkg_ids_1d = -1 * np.ones(len(Evt_Bkg), dtype=np.int32)

    event_energy = np.hstack((
        np.asarray(Evt_Src['SIGNAL']),
        np.asarray(Evt_Bkg['SIGNAL'])
    ))
    event_src_id = np.hstack((src_ids_1d, bkg_ids_1d))
    event_ra = np.hstack((
        np.asarray(Evt_Src['RA']),
        np.asarray(Evt_Bkg['RA'])
    ))
    event_dec = np.hstack((
        np.asarray(Evt_Src['DEC']),
        np.asarray(Evt_Bkg['DEC'])
    ))

elif Evt_Src is not None:
    event_energy = np.asarray(Evt_Src['SIGNAL'])
    event_src_id = ensure_1d_srcid(Evt_Src['SRC_ID']).astype(np.int32)
    event_ra = np.asarray(Evt_Src['RA'])
    event_dec = np.asarray(Evt_Src['DEC'])

else:
    event_energy = np.asarray(Evt_Bkg['SIGNAL'])
    event_src_id = -1 * np.ones(len(Evt_Bkg), dtype=np.int32)
    event_ra = np.asarray(Evt_Bkg['RA'])
    event_dec = np.asarray(Evt_Bkg['DEC'])

# ----------------------------------------------------------------------
# Build event trees
# ----------------------------------------------------------------------
if Evt_Src is not None:
    EvtSrcID32 = ensure_1d_srcid(Evt_Src['SRC_ID']).astype(np.int32)
    EvtRA_Src = np.asarray(Evt_Src['RA']) * np.pi / 180.0
    EvtDEC_Src = np.asarray(Evt_Src['DEC']) * np.pi / 180.0
    event_Npix_Src = np.asarray(Evt_Src['NPIXELS'])
    EvtCoord = np.transpose([EvtDEC_Src, EvtRA_Src])
    EvtTree = BallTree(EvtCoord, metric='haversine')
else:
    EvtSrcID32 = np.array([], dtype=np.int32)
    event_Npix_Src = np.array([], dtype=float)
    EvtTree = None

if Evt_Bkg is not None:
    EvtCoordBkg = np.transpose([
        np.asarray(Evt_Bkg['DEC']) * np.pi / 180.0,
        np.asarray(Evt_Bkg['RA']) * np.pi / 180.0
    ])
    EvtTreeBkg = BallTree(EvtCoordBkg, metric='haversine')
else:
    EvtTreeBkg = None

# ----------------------------------------------------------------------
# Counts per WRAPPED source id
# ----------------------------------------------------------------------
if len(EvtSrcID32) > 0:
    u, inv, photoncts = np.unique(EvtSrcID32, return_inverse=True, return_counts=True)
    EvtSrcTotCts = {int(u[i]): int(photoncts[i]) for i in range(len(u))}
else:
    EvtSrcTotCts = {}

# For each original simput source, count photons using the WRAPPED id.
# If two original sources collide to same int32, they will share these wrapped counts.
SimputSrcTotCts_wrap = np.array([EvtSrcTotCts.get(int(i), 0) for i in Simput_SrcID32], dtype=np.int32)

# ----------------------------------------------------------------------
# Save simput counts
# ----------------------------------------------------------------------
simput_srccts_path = os.path.join(workdir, 'simput_srccts_' + tile + '.fits')
#if not os.path.isfile(simput_srccts_path):
if len(EvtSrcID32) > 0:
    u, inv, photoncts = np.unique(EvtSrcID32, return_inverse=True, return_counts=True)
    EvtSrcCtsB1 = {int(u[i]): int(photoncts[i]) for i in range(len(u))}
else:
    EvtSrcCtsB1 = {}

SimputSrcCtsB1_wrap = np.array([EvtSrcCtsB1.get(int(i), 0) for i in Simput_SrcID32], dtype=np.int32)
ok = SimputSrcTotCts_wrap > 0

out = Table()
out['SRC_ID64'] = Simput_SrcID64[ok]
out['SRC_ID32'] = Simput_SrcID32[ok]
out['CtsTot_wrapid'] = SimputSrcTotCts_wrap[ok]
out['CtsB1_wrapid'] = SimputSrcCtsB1_wrap[ok]
out.write(simput_srccts_path, overwrite=True)
print('>>', simput_srccts_path)

# ----------------------------------------------------------------------
# Contamination file
# Keep the logic based on WRAPPED IDs, but write the primary source as int64.
# ----------------------------------------------------------------------
id_contam_path = os.path.join(workdir, 'ID_contam_' + tile + '.fits')
#if (EvtTree is not None) and (not os.path.isfile(id_contam_path)):
simOK = (SimputSrcTotCts_wrap > 0)

SimputSrcID64_sel = Simput_SrcID64[simOK]
SimputSrcID32_sel = Simput_SrcID32[simOK]
Simput_RA_sel = Simput_RA[simOK]
Simput_DEC_sel = Simput_DEC[simOK]
Simput_FLUX_sel = Simput_FLUX[simOK]
SimputSrcTotCts_sel = SimputSrcTotCts_wrap[simOK]

SimputCoord = np.transpose([Simput_DEC_sel * np.pi / 180.0, Simput_RA_sel * np.pi / 180.0])

FittingRad = 60.0
Contamination_MinCts = 3

EvtInd4Sim = EvtTree.query_radius(SimputCoord, r=FittingRad * np.pi / 180.0 / 3600.0)

blend_ID1 = []
blend_FLUX = []
blend_ID2_wrap = []
blend_cts1 = []
blend_cts2 = []

for n in range(len(SimputSrcID32_sel)):
    if len(EvtInd4Sim[n]) == 0:
        continue

    u, inv, photoncts = np.unique(EvtSrcID32[EvtInd4Sim[n]], return_inverse=True, return_counts=True)

    if SimputSrcID32_sel[n] not in u:
        continue

    eventcts = photoncts
    if len(u) <= 1:
        continue

    thisone = (u == SimputSrcID32_sel[n])
    if np.sum(thisone) != 1:
        continue

    ctsthisone = int(eventcts[thisone][0])
    other_cts = eventcts[~thisone]
    other_u = u[~thisone]

    second = np.argmax(other_cts)
    if other_cts[second] >= max(Contamination_MinCts, np.sqrt(ctsthisone)):
        second_id = int(other_u[second])
    else:
        second_id = int(-1 * other_u[second])

    blend_ID1.append(int(SimputSrcID64_sel[n]))
    blend_FLUX.append(float(Simput_FLUX_sel[n]))
    blend_ID2_wrap.append(second_id)
    blend_cts1.append(ctsthisone)
    blend_cts2.append(int(other_cts[second]))

out = Table()
out['ID_1'] = np.array(blend_ID1, dtype=np.int64)
out['FLUX'] = np.array(blend_FLUX, dtype=float)
out['ID_2_wrap32'] = np.array(blend_ID2_wrap, dtype=np.int32)
out['counts_1'] = np.array(blend_cts1, dtype=np.int32)
out['counts_2'] = np.array(blend_cts2, dtype=np.int32)
out.write(id_contam_path, overwrite=True)
print('>>', id_contam_path)


# ----------------------------------------------------------------------
# Matching routine for one source catalogue
# ----------------------------------------------------------------------
def run_match_for_catalogue(src_cat_path, prefix):
    if not os.path.isfile(src_cat_path):
        print('Missing source catalogue:', src_cat_path)
        return

    print(f'Now matching catalogue: {os.path.basename(src_cat_path)}')

    hdu = fits.open(src_cat_path)
    SrcCat = Table(hdu[1].data)

    AperRad = 20.0 * np.ones(len(SrcCat))
    AperRad[np.where(SrcCat['EXT_LIKE'] > 0)] = 60.0

    SrcCoord = np.transpose([SrcCat['DEC'] * np.pi / 180.0, SrcCat['RA'] * np.pi / 180.0])

    if EvtTree is not None:
        EvtInd4Cat = EvtTree.query_radius(SrcCoord, r=AperRad * np.pi / 180.0 / 3600.0)
    else:
        EvtInd4Cat = [np.array([], dtype=int) for _ in range(len(SrcCat))]

    if EvtTreeBkg is not None:
        EvtInd4Cat_bkg = EvtTreeBkg.query_radius(SrcCoord, r=AperRad * np.pi / 180.0 / 3600.0)
        SrcCat_N_Bkg = np.array([len(arr) for arr in EvtInd4Cat_bkg], dtype=np.int32)
    else:
        SrcCat_N_Bkg = np.zeros(len(SrcCat), dtype=np.int32)

    # Any match
    any_path = os.path.join(workdir, f'{prefix}_{tile}_IDMatch_Any_Tot{MinTotalCts}.fits')

    SrcID32 = np.zeros(len(SrcCoord), dtype=np.int32) - 99
    SrcApeCts = np.zeros(len(SrcCoord), dtype=np.int32) - 99
    SrcID32_2 = np.zeros(len(SrcCoord), dtype=np.int32) - 99
    SrcApeCts2 = np.zeros(len(SrcCoord), dtype=np.int32) - 99

    for n in range(len(SrcCat)):
        if len(EvtInd4Cat[n]) == 0:
            continue

        this_evt_ids = EvtSrcID32[EvtInd4Cat[n]]
        ok = np.array([EvtSrcTotCts.get(int(i), 0) for i in this_evt_ids]) >= MinTotalCts

        if np.sum(ok) == 0:
            continue

        u, inv, photoncts = np.unique(this_evt_ids[ok], return_inverse=True, return_counts=True)
        uc = np.array([[u[i], photoncts[i]] for i in range(len(u)) if photoncts[i] >= MinAperCts], dtype=np.int64)

        apebkg = int(SrcCat_N_Bkg[n])

        if len(uc) == 0:
            continue

        if len(uc) == 1:
            if uc[0, 1] + apebkg > stats.poisson.ppf(0.97725, apebkg):
                SrcID32[n] = np.int32(uc[0, 0])
                SrcApeCts[n] = np.int32(uc[0, 1])

        else:
            uc = uc[np.argsort(uc[:, 1])]
            if uc[-1, 1] + apebkg > stats.poisson.ppf(0.97725, apebkg):
                SrcID32[n] = np.int32(uc[-1, 0])
                SrcApeCts[n] = np.int32(uc[-1, 1])

                if uc[-2, 1] + apebkg > stats.poisson.ppf(0.97725, apebkg):
                    SrcID32_2[n] = np.int32(uc[-2, 0])
                    SrcApeCts2[n] = np.int32(uc[-2, 1])

    # Resolve wrapped IDs to original int64 IDs
    SrcID64 = np.zeros(len(SrcCoord), dtype=np.int64) - 99
    SrcID64_2 = np.zeros(len(SrcCoord), dtype=np.int64) - 99
    NCand = np.zeros(len(SrcCoord), dtype=np.int16)
    NCand2 = np.zeros(len(SrcCoord), dtype=np.int16)
    SepMatch = np.zeros(len(SrcCoord), dtype=float) + np.nan
    SepMatch2 = np.zeros(len(SrcCoord), dtype=float) + np.nan

    for n in range(len(SrcCat)):
        if SrcID32[n] >= 0:
            SrcID64[n], NCand[n], SepMatch[n] = unwrap_srcid_by_position(
                SrcID32[n], float(SrcCat['RA'][n]), float(SrcCat['DEC'][n]),
                wrap_to_simput_idx, Simput_SrcID64, Simput_RA, Simput_DEC
            )
        if SrcID32_2[n] >= 0:
            SrcID64_2[n], NCand2[n], SepMatch2[n] = unwrap_srcid_by_position(
                SrcID32_2[n], float(SrcCat['RA'][n]), float(SrcCat['DEC'][n]),
                wrap_to_simput_idx, Simput_SrcID64, Simput_RA, Simput_DEC
            )

    out_any = Table()
    out_any['ID_cat'] = SrcCat['ID_SRC']
    out_any['RA'] = SrcCat['RA']
    out_any['DEC'] = SrcCat['DEC']
    out_any['RADEC_ERR'] = SrcCat['RADEC_ERR']
    out_any['DET_LIKE_0'] = SrcCat['DET_LIKE_0']
    out_any['EXT'] = SrcCat['EXT']
    out_any['EXT_LIKE'] = SrcCat['EXT_LIKE']

    out_any['DIST_NN'] = SrcCat['DIST_NN']
    out_any['ML_CTS_0'] = SrcCat['ML_CTS_0']
    out_any['ML_CTS_ERR_0'] = SrcCat['ML_CTS_ERR_0']
    out_any['ML_RATE_0'] = SrcCat['ML_RATE_0']
    out_any['ML_RATE_ERR_0'] = SrcCat['ML_RATE_ERR_0']
    out_any['ML_FLUX_0'] = SrcCat['ML_FLUX_0']
    out_any['ML_FLUX_ERR_0'] = SrcCat['ML_FLUX_ERR_0']
    out_any['ML_BKG_0'] = SrcCat['ML_BKG_0']
    out_any['ML_EXP_1'] = SrcCat['ML_EXP_1']

    out_any['ID_simput'] = SrcID64
    out_any['aperture_counts'] = SrcApeCts
    out_any['ID_simput_2'] = SrcID64_2
    out_any['aperture_counts_2'] = SrcApeCts2

    out_any['ID_simput32'] = SrcID32
    out_any['ID_simput32_2'] = SrcID32_2
    out_any['N_wrap_candidates'] = NCand
    out_any['N_wrap_candidates_2'] = NCand2
    out_any['sep_match_arcsec'] = SepMatch
    out_any['sep_match_arcsec_2'] = SepMatch2

    out_any.write(any_path, overwrite=True)
    print('>>', any_path)

    # --------------------------------------------------------------
    # Unique match stage
    # Use the TRUE int64 IDs for uniqueness.
    # total_counts is tracked on WRAPPED ids because that is what events know.
    # --------------------------------------------------------------
    uniq_path = os.path.join(workdir, f'{prefix}_{tile}_IDMatch_Uniq_Tot{MinTotalCts}.fits')

    selected = np.zeros(len(SrcCoord), dtype=bool)

    valid_primary = SrcApeCts > 0
    s = np.argsort(SrcApeCts)[::-1]

    # unique on true int64 IDs
    _uID, ind = np.unique(SrcID64[s], return_index=True)
    selected[s[ind]] = True
    selected[~valid_primary] = False
    selected[SrcID64 < 0] = False

    for n in np.where(~selected & (SrcApeCts2 > 0))[0]:
        if SrcApeCts[n] <= 0:
            continue
        if np.sum(selected) == 0:
            continue
        if SrcID64[n] not in SrcID64[selected]:
            continue

        if np.any((SrcID64 == SrcID64[n]) & selected):
            if np.max(SrcApeCts[(SrcID64 == SrcID64[n]) & selected]) < SrcApeCts[n]:
                continue

        if SrcID64_2[n] not in SrcID64[selected] and SrcID64_2[n] >= 0:
            SrcID64[n] = SrcID64_2[n]
            SrcID32[n] = SrcID32_2[n]
            SrcApeCts[n] = SrcApeCts2[n]

            SrcID64_2[n] = -99
            SrcID32_2[n] = -99
            SrcApeCts2[n] = -99
            selected[n] = True

    SrcTotCts_wrap = np.zeros(len(SrcCoord), dtype=np.int32)

    for n in np.where(selected & (SrcID32 >= 0))[0]:
        SrcTotCts_wrap[n] = np.int32(EvtSrcTotCts.get(int(SrcID32[n]), 0))

    SrcID64_out = SrcID64.copy()
    SrcApeCts_out = SrcApeCts.copy()
    SrcID32_out = SrcID32.copy()

    SrcID64_out[~selected] = -99
    SrcID32_out[~selected] = -99
    SrcApeCts_out[~selected] = -99

    out_uniq = Table()
    out_uniq['ID_cat'] = SrcCat['ID_SRC']
    out_uniq['RA'] = SrcCat['RA']
    out_uniq['DEC'] = SrcCat['DEC']
    out_uniq['RADEC_ERR'] = SrcCat['RADEC_ERR']
    out_uniq['DET_LIKE_0'] = SrcCat['DET_LIKE_0']
    out_uniq['EXT'] = SrcCat['EXT']
    out_uniq['EXT_LIKE'] = SrcCat['EXT_LIKE']

    out_uniq['DIST_NN'] = SrcCat['DIST_NN']
    out_uniq['ML_CTS_0'] = SrcCat['ML_CTS_0']
    out_uniq['ML_CTS_ERR_0'] = SrcCat['ML_CTS_ERR_0']
    out_uniq['ML_RATE_0'] = SrcCat['ML_RATE_0']
    out_uniq['ML_RATE_ERR_0'] = SrcCat['ML_RATE_ERR_0']
    out_uniq['ML_FLUX_0'] = SrcCat['ML_FLUX_0']
    out_uniq['ML_FLUX_ERR_0'] = SrcCat['ML_FLUX_ERR_0']
    out_uniq['ML_BKG_0'] = SrcCat['ML_BKG_0']
    out_uniq['ML_EXP_1'] = SrcCat['ML_EXP_1']

    out_uniq['ID_simput'] = SrcID64_out
    out_uniq['ID_simput32'] = SrcID32_out
    out_uniq['aperture_counts'] = SrcApeCts_out
    out_uniq['total_counts_wrapid'] = SrcTotCts_wrap
    out_uniq['N_wrap_candidates'] = NCand
    out_uniq['sep_match_arcsec'] = SepMatch

    out_uniq.write(uniq_path, overwrite=True)
    print('>>', uniq_path)


# ----------------------------------------------------------------------
# Run both catalogues
# ----------------------------------------------------------------------
run_match_for_catalogue(p_2_esass, 'Sc1')
run_match_for_catalogue(p_2_esass_Lext0, 'Sc_Lext0')
