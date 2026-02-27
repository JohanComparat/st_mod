import os, glob
from astropy.wcs import WCS
import astropy.io.fits as fits
import time
from astropy.table import Table, vstack
import numpy as np
import sys

#Set paths and environmental variables
basedir = '/home/idies/workspace/erosim/Uchuu/LCerass'
base_st_mod_data = '/home/idies/workspace/erosim/software/st_mod_data'
os.environ['UCHUU']='/home/idies/workspace/erosim/Uchuu'
os.environ['GIT_STMOD']='/home/idies/workspace/erosim/software/st_mod'
os.environ['GIT_STMOD_DATA']='/home/idies/workspace/erosim/software/st_mod_data'
LC_dir = 'LCerass'
top_dir = os.path.join(os.environ['UCHUU'], LC_dir)

#Will be used to get local Bg and ExpMap to rescale pdet
p2erass1clu = '/home/idies/workspace/erosim/eRASS1_mock_input_CLUSTERS_extlike3.fits'
erass1_clu = Table.read(p2erass1clu, memmap=True)

agn_seed = sys.argv[1] # Seed for the AGN realization e.g. 1
clu_seed = sys.argv[2] # Seed for the CLU realization e.g. 1
exp_name = sys.argv[3] # Type of the experiment, whether e4 or e5
if exp_name == 'e4':
    real_data_name = 's4'
elif exp_name == 'e5':
    real_data_name = 's5'
mergeType = 'GE'

nl = lambda sel: len(sel.nonzero()[0])

sky_map_hdu = Table.read(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS.fits'))

def get_srvmap(ra, dec):
    return sky_map_hdu['SRVMAP'].value[
        (sky_map_hdu['RA_MIN'] < ra) & (sky_map_hdu['RA_MAX'] >= ra) & (sky_map_hdu['DE_MIN'] < dec) & (
                    sky_map_hdu['DE_MAX'] >= dec)]

def remove_events(evlist, blacklist):
    '''Events from evlist will be removed following the energy statistics
    in blacklist (at least roughly)
    '''
    e_blacklist = blacklist['SIGNAL']
    for bl_e in e_blacklist:
        e_evlist = evlist['SIGNAL']
        e_evlist_argsort = np.argsort(e_evlist)
        e_evlist_sorted = e_evlist[e_evlist_argsort]
        pos = np.searchsorted(e_evlist_sorted, bl_e)
        print("Index to remove: ", e_evlist_argsort[pos])
        tokeep = np.delete(e_evlist_argsort, pos)
        evlist = evlist[tokeep]
    # Reshuffle list otherwise it outputs sorted by energy
    id_shuffle = np.arange(len(evlist))
    np.random.shuffle(id_shuffle)
    return evlist[id_shuffle]

def remove_events_binned(evlist, blacklist, nbins=400, emin=None, emax=None, rng=None):
    """
    Remove events from evlist so that the removed set has approximately the
    same SIGNAL distribution as blacklist (in nbins).
    """
    if rng is None:
        rng = np.random.default_rng()

    e_ev = np.asarray(evlist['SIGNAL'])
    e_bl = np.asarray(blacklist['SIGNAL'])

    if emin is None:
        emin = min(e_ev.min(), e_bl.min())
    if emax is None:
        emax = max(e_ev.max(), e_bl.max())
    if not np.isfinite(emin) or not np.isfinite(emax) or emin >= emax:
        return evlist  # nothing sensible to do

    edges = np.linspace(emin, emax, nbins + 1)

    # how many to remove per bin
    bl_counts, _ = np.histogram(e_bl, bins=edges)

    # indices of evlist in each bin
    bin_id = np.searchsorted(edges, e_ev, side="right") - 1
    valid = (bin_id >= 0) & (bin_id < nbins)
    idx = np.nonzero(valid)[0]
    bin_id = bin_id[valid]

    # collect indices to remove
    to_remove = []
    # group indices by bin (fast-ish)
    for b in range(nbins):
        k = bl_counts[b]
        if k <= 0:
            continue
        in_bin = idx[bin_id == b]
        if in_bin.size == 0:
            continue
        k = min(k, in_bin.size)
        # sample k indices uniformly within bin
        to_remove.append(rng.choice(in_bin, size=k, replace=False))

    if not to_remove:
        return evlist

    to_remove = np.concatenate(to_remove)
    keep_mask = np.ones(len(evlist), dtype=bool)
    keep_mask[to_remove] = False

    # return shuffled (to avoid energy ordering)
    kept_idx = np.nonzero(keep_mask)[0]
    rng.shuffle(kept_idx)
    return evlist[kept_idx]

def ctr_BKG_percent_detection(bg_cts_per_pix, b, a):
    '''
    from erass1 mock
    AGN 10% -0.86810338809034, -4.02348244025905
    CLU 10% -0.74810338809034, -3.29348244025905
    AGN 50% -0.76810338809034,-3.62348244025905
    CLU 50% -0.59810338809034,-2.39348244025905
    '''
    return 10 ** a * bg_cts_per_pix ** b

print('{0} tiles to process'.format(len(sky_map_hdu[(sky_map_hdu['OWNER'] == 2) | (sky_map_hdu['OWNER'] == 0)])))

fails = []
good = []
for sky_tile in sky_map_hdu[(sky_map_hdu['OWNER'] == 2) | (sky_map_hdu['OWNER'] == 0)]:#[:1]:

    #Set reference time
    t0 = time.time()

    sky_tile_id = str(sky_tile['SRVMAP'])
    str_field = str(sky_tile['SRVMAP']).zfill(6)
    print('\nTile {0} - Started processing'.format(str_field))
    evt_list = np.array(glob.glob(
        os.path.join(os.environ['UCHUU'], LC_dir, str_field, '{0}_c030'.format(real_data_name),
                     '*_Image_c030.fits.gz')))

    p2_real_BgMap = os.path.join(os.environ['UCHUU'], LC_dir, str_field, '{0}_eSASS'.format(real_data_name),
                     str_field+'_024_Bg3Map.fits')
    p2_real_ExpMap = os.path.join(os.environ['UCHUU'], LC_dir, str_field, '{0}_eSASS'.format(real_data_name),
                     str_field+'_024_ExpMap.fits')

    real_BgMap = fits.open(p2_real_BgMap)[0].data
    real_ExpMap = fits.open(p2_real_ExpMap)[0].data

    esass_dir = os.path.join(os.environ['UCHUU'], LC_dir, str_field,
                             mergeType + '_{0}_merge_AGNseed'.format(exp_name) + agn_seed.zfill(
                                 3) + '_SimBKG_CLUseed' + clu_seed.zfill(3))
    if not os.path.isdir(esass_dir):
        os.system('mkdir -p ' + esass_dir)
    # esass_GE_dir = os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'GE_{0}_merge_AGNseed'.format(exp_name)+agn_seed.zfill(3)+'_SimBKG_CLUseed'+clu_seed.zfill(3))

    path_2_event_file = os.path.join(esass_dir, 'evt_' + str_field + '.fits')
    path_2_simeventAGN_file = os.path.join(esass_dir, 'simAGNevt_' + str_field + '.fits')
    path_2_simeventCLU_file = os.path.join(esass_dir, 'simCLUevt_' + str_field + '.fits')
    path_2_simeventBKG_file = os.path.join(esass_dir, 'simBKGevt_' + str_field + '.fits')
    if len(evt_list) == 0 or os.path.isfile(path_2_event_file):
        print(
            'Tile {0} -  Not processed. Reason:\n no real data event file: {1}\n merged simulated event file exists: {2}'.format(
                str_field, len(evt_list) == 0, os.path.isfile(path_2_event_file)))
        fails.append(1)
        continue
    bg_dir = os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'pBG2')  # 'evt_particle_???.fits' )
    BG_evt_files = np.array(glob.glob(os.path.join(bg_dir, '*.fits')))
    if len(BG_evt_files) == 0:
        print('\nTile {0} - Not processed. Reason:\n no background file: {1}'.format(str_field, len(BG_evt_files) == 0))
        fails.append(2)
        continue
    hdul_raw = fits.open(evt_list[0])
    try:
        gtis = np.array([np.sum(hdul_raw['GTI1'].data['STOP'] - hdul_raw['GTI1'].data['START'])
                            , np.sum(hdul_raw['GTI2'].data['STOP'] - hdul_raw['GTI2'].data['START'])
                            , np.sum(hdul_raw['GTI3'].data['STOP'] - hdul_raw['GTI3'].data['START'])
                            , np.sum(hdul_raw['GTI4'].data['STOP'] - hdul_raw['GTI4'].data['START'])
                            , np.sum(hdul_raw['GTI5'].data['STOP'] - hdul_raw['GTI5'].data['START'])
                            , np.sum(hdul_raw['GTI6'].data['STOP'] - hdul_raw['GTI6'].data['START'])
                            , np.sum(hdul_raw['GTI7'].data['STOP'] - hdul_raw['GTI7'].data['START'])])
    except KeyError:
        print('\nTile {0} - Has a problem: KeyError'.format(str_field))
        fails.append(3)
        continue
    # for eee in hdul_raw[1:]:
    # 	print(eee.header['EXTNAME'])
    hdul = fits.open(evt_list[0])
    hdul['EVENTS'].data['RA'][hdul['EVENTS'].data['RA'] == 0] = 1e-6
    #no need to recomupte skytile, we already have ra dec borders in sky_tile table
    #SRV_ev = np.array([get_srvmap(e0, e1) for e0, e1 in zip(hdul['EVENTS'].data['RA'], hdul['EVENTS'].data['DEC'])])
    #to_replace = (np.hstack((SRV_ev)) == sky_tile['SRVMAP'])
    #ids_to_replace = np.arange(len(to_replace))[to_replace]
    #N_ev_OBS = len(ids_to_replace)
    to_replace = (
            (hdul['EVENTS'].data['RA'] > sky_tile['RA_MIN']) & (hdul['EVENTS'].data['RA'] <= sky_tile['RA_MAX']) &
            (hdul['EVENTS'].data['DEC'] > sky_tile['DE_MIN']) & (hdul['EVENTS'].data['DEC'] <= sky_tile['DE_MAX'])
    )
    ids_to_replace = np.flatnonzero(to_replace)
    N_ev_OBS = ids_to_replace.size

    agn_dir = os.path.join(os.environ['UCHUU'], LC_dir, str_field,
                           "eRASS8_SEED_" + str(agn_seed).zfill(3) + "_events_AGN_2025_04")
    cluster_dir = os.path.join(os.environ['UCHUU'], LC_dir, str_field,
                               "Att_eRASS8_sixte_v27_SEED_" + str(clu_seed).zfill(3) + "_events_cluster_Xgas_bHS0.8")
    stars_dir = os.path.join(os.environ['UCHUU'], LC_dir, str_field, 'stars')  # , 'simulated_photons_ccd?.fits' )

    # Get the AGN and CLU event lists for the exact survey GTI
    data_A = []
    data_C = []
    # data_B = []
    t_max_A = []
    t_max_C = []
    frac_A = []
    frac_C = []
    for NCCD, t_obs in zip(np.arange(7) + 1, gtis):
        agn_evt_files = np.array(glob.glob(os.path.join(agn_dir, 't0erass_ccd' + str(NCCD) + '_evt.fits')))
        hdu_A = fits.open(agn_evt_files[0])
        texp_cum = np.cumsum(hdu_A['STDGTI'].data['STOP'] - hdu_A['STDGTI'].data['START'])
        t_max = hdu_A['STDGTI'].data['STOP'][np.searchsorted(texp_cum, t_obs) + 1]
        t_max_A.append(t_max)
        t_a = Table(hdu_A['EVENTS'].data)
        frac_A.append(nl(t_a['TIME'] < t_max) / len(t_a))
        t_a['TM_NR'] = NCCD
        t_a = t_a[t_a['TIME'] < t_max]
        data_A.append(t_a)

        CL_evt_files = np.array(glob.glob(os.path.join(cluster_dir, 't0erass*ccd' + str(NCCD) + '_evt.fits')))
        hdu_C = fits.open(CL_evt_files[0])
        texp_cum = np.cumsum(hdu_C['STDGTI'].data['STOP'] - hdu_C['STDGTI'].data['START'])
        # t_max = hdu_C['STDGTI'].data['STOP'][np.searchsorted(texp_cum, t_obs)+1]
        # t_max_C.append(t_max)
        t_c = Table(hdu_C['EVENTS'].data)
        frac_C.append(nl(t_c['TIME'] < t_max) / len(t_c))
        t_c['TM_NR'] = NCCD
        t_c = t_c[t_c['TIME'] < t_max]
        data_C.append(t_c)

    t_max_A = np.array(t_max_A)
    # t_max_C = np.array(t_max_C)
    frac_A = np.array(frac_A)
    frac_C = np.array(frac_C)
    data_A = vstack((data_A))
    data_C = vstack((data_C))
    data_A['is_in_unique_area'] = (data_A['RA'] >= sky_tile['RA_MIN']) & (data_A['RA'] <= sky_tile['RA_MAX']) & (
                data_A['DEC'] >= sky_tile['DE_MIN']) & (data_A['DEC'] <= sky_tile['DE_MAX'])
    data_C['is_in_unique_area'] = (data_C['RA'] >= sky_tile['RA_MIN']) & (data_C['RA'] <= sky_tile['RA_MAX']) & (
                data_C['DEC'] >= sky_tile['DE_MIN']) & (data_C['DEC'] <= sky_tile['DE_MAX'])
    # now in data_A, data_C I have agn and clu events reteained after the gti filtering as in the real data plus a flag about the unique tile area
    # need to determine how many of these events are sources and how many are background

    # Oversample the BKG [from eRASS1], CLU, AGN event lists
    oversampling = 1.5
    bg_all = []
    frac_B = []
    for el in BG_evt_files:
        hh = fits.open(el)
        exp_bg = hh[1].header['EXPOSURE']
        tt0 = Table.read(el)
        NCCD = int(os.path.basename(el).split('_')[-2][-1])
        tt0['TM_NR'] = NCCD
        tt0.keep_columns(['TIME', 'RA', 'DEC', 'RAWX', 'RAWY', 'PHA', 'SIGNAL', 'TM_NR'])
        tt0 = tt0[tt0['TIME'] < np.max(t_max_A) * oversampling]  # to have too many BG events
        # time_sel_bg = ( np.random.uniform(0,1,len(tt0)) < 1/4 ) # (gtis[NCCD-1]/exp_bg) )
        # frac_B.append( nl( time_sel_bg ) / len(tt0) )
        # tt0=tt0[time_sel_bg]
        bg_all.append(tt0)

    frac_B = np.array(frac_B)
    data_B_oversampled = vstack((bg_all))
    id_B_shuffle = np.arange(len(data_B_oversampled))
    np.random.shuffle(id_B_shuffle)
    data_B_oversampled = data_B_oversampled[id_B_shuffle]
    data_B_oversampled['is_in_unique_area'] = (data_B_oversampled['RA'] >= sky_tile['RA_MIN']) & (
                data_B_oversampled['RA'] <= sky_tile['RA_MAX']) & (data_B_oversampled['DEC'] >= sky_tile['DE_MIN']) & (
                                                          data_B_oversampled['DEC'] <= sky_tile['DE_MAX'])
    # already trimmed to the unique area
    data_B_oversampled.write(os.path.join(esass_dir, 'all_BG_evts_oversampled.fits'), overwrite=True)

    #do the same oversampling for AGN and CLU
    clu_oversampled = []
    agn_oversampled = []
    for NCCD in np.arange(7) + 1:
        agn_evt_files = np.array(glob.glob(os.path.join(agn_dir, 't0erass_ccd' + str(NCCD) + '_evt.fits')))
        hdu_A = Table.read(agn_evt_files[0])
        tt0 = hdu_A[hdu_A['TIME'] < np.max(t_max_A) * oversampling]
        agn_oversampled.append(tt0)
        #
        clu_evt_files = np.array(glob.glob(os.path.join(cluster_dir, 't0erass_ccd' + str(NCCD) + '_evt.fits')))
        hdu_C = Table.read(clu_evt_files[0])
        tt0 = hdu_C[hdu_C['TIME'] < np.max(t_max_A) * oversampling]
        clu_oversampled.append(tt0)
    data_A_oversampled = vstack(agn_oversampled)
    id_A_shuffle = np.arange(len(data_A_oversampled))
    np.random.shuffle(id_A_shuffle)
    data_A_oversampled = data_A_oversampled[id_A_shuffle]
    data_A_oversampled['is_in_unique_area'] = (data_A_oversampled['RA'] >= sky_tile['RA_MIN']) & (
                data_A_oversampled['RA'] <= sky_tile['RA_MAX']) & (data_A_oversampled['DEC'] >= sky_tile['DE_MIN']) & (
                                                          data_A_oversampled['DEC'] <= sky_tile['DE_MAX'])
    data_C_oversampled = vstack((clu_oversampled))
    id_C_shuffle = np.arange(len(data_C_oversampled))
    np.random.shuffle(id_C_shuffle)
    data_C_oversampled = data_C_oversampled[id_C_shuffle]
    data_C_oversampled['is_in_unique_area'] = (data_C_oversampled['RA'] >= sky_tile['RA_MIN']) & (
                data_C_oversampled['RA'] <= sky_tile['RA_MAX']) & (data_C_oversampled['DEC'] >= sky_tile['DE_MIN']) & (
                                                          data_C_oversampled['DEC'] <= sky_tile['DE_MAX'])

    print('\nTile {0} - Number of:\n AGNs {1}\n Clusters {2}\n Oversampled background {3}\n Total simulated data {4}\n Number of observed events {5}'.format(str_field, len(data_A), len(data_C), len(data_B_oversampled), len(data_A) + len(data_C) + len(data_B_oversampled), N_ev_OBS))

    # verify the background levels between the simulated and real data
    # exposure and background maps (from the real data)
    expmap = fits.open(os.path.join(os.environ['UCHUU'], LC_dir, str_field,
                                    real_data_name + '_eSASS/' + str_field + '_024_ExpMap.fits'))
    bg3mapD = fits.open(os.path.join(os.environ['UCHUU'], LC_dir, str_field,
                                     real_data_name + '_eSASS/' + str_field + '_024_Bg3Map.fits'))
    local_BKG = np.median(bg3mapD[0].data)
    local_EXP = np.median(expmap[0].data)

    sel_local = erass1_clu['tile'] == str_field
    local_BKG_erass1 = np.median(erass1_clu['BG3Model'][sel_local])
    local_EXP_erass1 = np.median(erass1_clu['TexpModel'][sel_local])

    # select events that will belong top detected sources
    # deduce events that will belong to the background
    texp_factor_to_e1 = float(local_EXP/local_EXP_erass1)
    print('\nTile {0} - texp_factor_to_e1: {1}'.format(str_field, texp_factor_to_e1))

    A_DetLimit_eRASS1 = ctr_BKG_percent_detection(local_BKG / texp_factor_to_e1, -0.76810338809034,-3.62348244025905) * local_EXP / texp_factor_to_e1
    C_DetLimit_eRASS1 = ctr_BKG_percent_detection(local_BKG / texp_factor_to_e1, -0.59810338809034,-2.39348244025905) * local_EXP / texp_factor_to_e1
    print('\nTile {0} - A_DetLimit_eRASS1: {1}'.format(str_field, A_DetLimit_eRASS1))
    print('\nTile {0} - C_DetLimit_eRASS1: {1}'.format(str_field, C_DetLimit_eRASS1))

    # To erass4 the SNR of count rate goes down by sqrt(4): you get higher counts by a factor of 2, but BG increased by 4, so effectively SNR went down.
    # A source at the limit producing 1 ct in erass1 will produce for counts in erass4, where the limit is now 2 counts, so you need more counts in erass4 than erass1 but the same source has higher significance in erass4 than erass1
    A_DetLimit_eRASSn = ctr_BKG_percent_detection(local_BKG / texp_factor_to_e1, -0.76810338809034,-3.62348244025905) / np.sqrt(texp_factor_to_e1) * local_EXP
    C_DetLimit_eRASSn = ctr_BKG_percent_detection(local_BKG / texp_factor_to_e1, -0.59810338809034,-2.39348244025905) / np.sqrt(texp_factor_to_e1) * local_EXP
    print('\nTile {0} - A_DetLimit_eRASSn: {1}'.format(str_field, A_DetLimit_eRASSn))
    print('\nTile {0} - C_DetLimit_eRASSn: {1}'.format(str_field, C_DetLimit_eRASSn))

    #A_DetLimit_eRASS1 = 4
    #C_DetLimit_eRASS1 = 8
    #from the oversampled pool of agn and clu events identifiy the ones that should fall as part of the BKG
    src_id_A, n_evt_A = np.unique(data_A_oversampled['SRC_ID'].T[0], return_counts=True)
    #not_detected_src_id_A = src_id_A[n_evt_A < A_DetLimit_eRASS1]
    not_detected_src_id_A = src_id_A[n_evt_A < A_DetLimit_eRASSn]
    selection_bg_A = np.isin(data_A_oversampled['SRC_ID'].T[0], not_detected_src_id_A)
    BG_A = data_A_oversampled[selection_bg_A]
    src_id_C, n_evt_C = np.unique(data_C_oversampled['SRC_ID'].T[0], return_counts=True)
    #not_detected_src_id_C = src_id_C[n_evt_C < C_DetLimit_eRASS1]
    not_detected_src_id_C = src_id_C[n_evt_C < C_DetLimit_eRASSn]
    selection_bg_C = np.isin(data_C_oversampled['SRC_ID'].T[0], not_detected_src_id_C)
    BG_C = data_C_oversampled[selection_bg_C]
    data_A_oversampled['is_BG'] = selection_bg_A
    data_C_oversampled['is_BG'] = selection_bg_C
    data_A_oversampled['is_det'] = ~selection_bg_A
    data_C_oversampled['is_det'] = ~selection_bg_C
    data_A_oversampled.write(os.path.join(esass_dir, 'all_AGN_evts_oversampled.fits'), overwrite=True)
    data_C_oversampled.write(os.path.join(esass_dir, 'all_GAS_evts_oversampled.fits'), overwrite=True)

    # Remove statistically the A and C is_BG from the oversampled data_B
    data_B_oversampled = remove_events_binned(data_B_oversampled, data_A_oversampled[data_A_oversampled['is_BG']])
    data_B_oversampled = remove_events_binned(data_B_oversampled, data_C_oversampled[data_C_oversampled['is_BG']])
    # Now we have removed from the oversampled background list contributions
    # from AGN and CLU (statistically speaking) that were undetected by eRASS1   RS: I don't think this is needed, we can use directly the estimated erassn limit
    # ==>  data_B_oversampled contains pure eRASS1 background rate, with exposure  #RS: so now this is oversampled BKG than contains expected erassN bkg rate with Texp=1.5*TexpN
    # time that is 1.5 longer than the eRASS:n exposure time.

    # We want to match the level of the eRASS:n background in the 0.2-2.3 keV
    # with the BgMap3 map.
    # It will take CLU and AGN up to the eRASS:n detection limit and complete
    # with events from the clean oversampled background list

    # verify the background levels between the simulated and real data
    # select events that will belong top detected sources
    # deduce events that will belong to the background

    #A_DetLimit_eRASSn = 3
    #C_DetLimit_eRASSn = 7
    #now get AGN and CLU evts from the right exptime window and add flag to the ones expected to fall back as BKG
    src_id_A, n_evt_A = np.unique(data_A['SRC_ID'].T[0], return_counts=True)
    not_detected_src_id_A = src_id_A[n_evt_A < A_DetLimit_eRASSn]
    selection_bg_A = np.isin(data_A['SRC_ID'].T[0], not_detected_src_id_A)
    BG_A = data_A[selection_bg_A]
    src_id_C, n_evt_C = np.unique(data_C['SRC_ID'].T[0], return_counts=True)
    not_detected_src_id_C = src_id_C[n_evt_C < C_DetLimit_eRASSn]
    selection_bg_C = np.isin(data_C['SRC_ID'].T[0], not_detected_src_id_C)
    BG_C = data_C[selection_bg_C]
    data_A['is_BG'] = selection_bg_A
    data_C['is_BG'] = selection_bg_C
    data_A['is_det'] = ~selection_bg_A
    data_C['is_det'] = ~selection_bg_C
    data_A.write(os.path.join(esass_dir, 'all_AGN_evts.fits'), overwrite=True)
    data_C.write(os.path.join(esass_dir, 'all_GAS_evts.fits'), overwrite=True)

    #then match to 0.2-2.3 bgmap
    emin_bgmap = 0.2
    emax_bgmap = 2.3
    wcs = WCS(expmap[0].header)
    shape_xy = bg3mapD[0].data.shape
    pix_mat_X, pix_mat_Y = np.meshgrid(np.arange(shape_xy[0]), np.arange(shape_xy[1]))
    pix_mat_X_ravel, pix_mat_Y_ravel = pix_mat_X.ravel(), pix_mat_Y.ravel()
    rade_mat = wcs.pixel_to_world(pix_mat_X, pix_mat_Y)
    is_in_unique_area = (rade_mat.ra.deg > sky_tile['RA_MIN']) & (rade_mat.ra.deg <= sky_tile['RA_MAX']) & (
                rade_mat.dec.deg > sky_tile['DE_MIN']) & (rade_mat.dec.deg <= sky_tile['DE_MAX'])
    bg_unique_area = bg3mapD[0].data[pix_mat_X[is_in_unique_area], pix_mat_Y[is_in_unique_area]]
    print('\nTile {0} - BG in the unique area:\n min {1}\n max {2}\n median {3}\n mean {4}\n std {5}'.format(str_field, bg_unique_area.min(), bg_unique_area.max(), np.median(bg_unique_area), np.mean(bg_unique_area), np.std(bg_unique_area)))
    BG_CT_val_target = np.median(bg_unique_area)  # This value is the eRASS:n background (including CLU and AGN below eRASS:n threshold)

    # locate background event on the BG map :
    x_pix_B, y_pix_B = wcs.wcs_world2pix(data_B_oversampled['RA'], data_B_oversampled['DEC'], 0)
    x_pix_B = np.round(x_pix_B).astype(int)
    y_pix_B = np.round(y_pix_B).astype(int)
    # all possible pixels in the unique area
    all_x_pix = np.arange(x_pix_B.min(), x_pix_B.max(), 1)
    all_y_pix = np.arange(y_pix_B.min(), y_pix_B.max(), 1)
    # all_X,all_Y = np.meshgrid(all_x_pix, all_y_pix)
    # all_pixels = np.transpose([all_X.ravel(), all_Y.ravel()])
    # N_pixels = all_pixels.shape[0]
    N_pixels = len(all_x_pix) * len(all_y_pix)
    #take BG candidates from AGN and CLU in the 0.2-2.3 range
    BGcandidate_in_A_eRange = data_A['is_BG'] & (data_A['SIGNAL'] > emin_bgmap) & (data_A['SIGNAL'] < emax_bgmap)
    BGcandidate_in_C_eRange = data_C['is_BG'] & (data_C['SIGNAL'] > emin_bgmap) & (data_C['SIGNAL'] < emax_bgmap)
    BGcandidate_in_B_eRange = (data_B_oversampled['SIGNAL'] > emin_bgmap) & (data_B_oversampled['SIGNAL'] < emax_bgmap)

    N_BGcandidate_in_A_eRange = np.count_nonzero(BGcandidate_in_A_eRange)
    N_BGcandidate_in_C_eRange = np.count_nonzero(BGcandidate_in_C_eRange)
    N_BGcandidate_in_B_eRange = np.count_nonzero(BGcandidate_in_B_eRange)
    N_BG_prediction = N_BGcandidate_in_A_eRange + N_BGcandidate_in_B_eRange + N_BGcandidate_in_C_eRange
    print('\nTile {0} - N_BGcandidate_in_A_eRange {1}, N_BGcandidate_in_B_eRange {2}, N_BGcandidate_in_C_eRange {3}'.format( str_field, N_BGcandidate_in_A_eRange, N_BGcandidate_in_B_eRange, N_BGcandidate_in_C_eRange))
    print('\nTile {0} - N_BG_prediction {1}'.format(str_field, N_BG_prediction))
    print('\nTile {0} - BG CT per pixel predicted by the model {1}'.format(str_field, N_BG_prediction / N_pixels))
    print('\nTile {0} - Median observed BG CT per pixel in the data {1}'.format(str_field, BG_CT_val_target))
    print('\nTile {0} - We predict {1} We need {2}'.format(str_field, N_BG_prediction, int(N_pixels * BG_CT_val_target)))

    if N_BG_prediction < int(N_pixels * BG_CT_val_target):
        print('\nTile {0} - Need to add more BG events. This should never be the case !'.format(str_field))
    else:
        print('\nTile {0} - Need to remove some BG events from the oversampled background list'.format(str_field))
        # Calculate the number of extra BG events in (emin, emax)
        Nextra_eRange = int(N_pixels * BG_CT_val_target) - N_BG_prediction
        # Assume all these extra events come from the BG file.
        # Apply correction by assuming a constant removal factor throughout energy distribution
        ratio_full_over_eRange = len(data_B_oversampled) / N_BGcandidate_in_B_eRange
        Nextra_fullRange = int(Nextra_eRange * ratio_full_over_eRange)
        # Remove Nextra_fullRange events from the oversampled clean BG list
        data_B = data_B_oversampled[:-Nextra_fullRange]

    # From here it continues as in original JC's code.......
    # ##
    # ##
    # ##
    print('\nTile {0} - Number of:\n AGNs {1}\n Clusters {2}\n Background {3}\n Total simulated data {4}\n Number of observed events {5}'.format(str_field, len(data_A), len(data_C), len(data_B), len(data_A) + len(data_C) + len(data_B), N_ev_OBS))
    if len(data_A) + len(data_C) + len(data_B) == N_ev_OBS:
        print('\nTile {0} - Exactly the number of events needed, perfect match!'.format(str_field))
        fi_up = ['RA', 'DEC', 'RAWX', 'RAWY', 'PHA']
        for fn in fi_up:
            hdul['EVENTS'].data[fn][ids_to_replace] = np.hstack((data_C[fn], data_A[fn], data_B[fn]))
        fn = 'SIGNAL'
        hdul['EVENTS'].data['PI'][ids_to_replace] = 1000. * np.hstack((data_C[fn], data_A[fn], data_B[fn]))
        hdul.writeto(path_2_event_file, overwrite=True)
        print('\nTile {0} - File written to:\n {1}'.format(str_field, path_2_event_file))

    elif len(data_A) + len(data_C) + len(data_B) < N_ev_OBS:
        N_available = len(data_A) + len(data_C) + len(data_B)
        N_too_many = N_ev_OBS - N_available
        print('\nTile {0} - Less simulated events than observed events (difference is {1}), need to move some true events out of the unique area.'.format(str_field, N_too_many))
        ids_to_replace = np.arange(len(to_replace))[to_replace]
        np.random.shuffle(ids_to_replace)

        fi_up = ['RA', 'DEC', 'RAWX', 'RAWY', 'PHA']
        for fn in fi_up:
            hdul['EVENTS'].data[fn][ids_to_replace[:N_available]] = np.hstack((data_C[fn], data_A[fn], data_B[fn]))
            hdul['EVENTS'].data[fn][ids_to_replace[N_available:]] = hdul['EVENTS'].data[fn].min() * np.ones(N_too_many)
        fn = 'SIGNAL'
        hdul['EVENTS'].data['PI'][ids_to_replace[:N_available]] = 1000. * np.hstack(
            (data_C[fn], data_A[fn], data_B[fn]))
        hdul['EVENTS'].data['PI'][ids_to_replace[N_available:]] = hdul['EVENTS'].data['PI'].min() * np.ones(N_too_many)
        hdul.writeto(path_2_event_file, overwrite=True)
        print('\nTile {0} - File written to:\n {1}'.format(str_field, path_2_event_file))

    else:
        N_additional = len(data_A) + len(data_C) + len(data_B) - N_ev_OBS
        print('\nTile {0} - More simulated events than observed events (difference is {1}), need to take some true events outside of the unique area.'.format(str_field, N_additional))
        print('\nTile {0} - Nsim {1} Nobs {2} Nadditional {3}'.format(str_field, len(data_A) + len(data_C) + len(data_B), N_ev_OBS, N_additional))
        extra_ids = np.arange(len(to_replace))[np.isin(np.arange(len(to_replace)), ids_to_replace, invert=True)]
        np.random.shuffle(extra_ids)
        ids_to_replace2 = np.hstack((np.arange(len(to_replace))[to_replace], extra_ids[:N_additional]))
        enough_events = len(ids_to_replace2)-(len(data_A) + len(data_C) + len(data_B))>0

        fi_up = ['RA', 'DEC', 'RAWX', 'RAWY', 'PHA']
        if enough_events:
            for fn in fi_up:
                hdul['EVENTS'].data[fn][ids_to_replace2] = np.hstack((data_C[fn], data_A[fn], data_B[fn]))
            fn = 'SIGNAL'
            hdul['EVENTS'].data['PI'][ids_to_replace2] = 1000. * np.hstack((data_C[fn], data_A[fn], data_B[fn]))

        else:
            print('\nTile {0} - Not enough events outside the unique area. Need to create a new bin table to accomodate all simulated events.'.format(str_field))
            N_sim = len(data_A) + len(data_C) + len(data_B)
            DIFF = N_sim - ids_to_replace2.size
            print('\nTile {0} - Not enough rows to overwrite: need {1}, have {2}. Growing EVENTS by {3} rows.'.format(str_field, N_sim,ids_to_replace2.size, DIFF))
            #grow EVENTS table by DIFF rows by copying last rows as template
            old = hdul['EVENTS'].data
            n_old = len(old)
            n_new = n_old + DIFF

            new = np.empty(n_new, dtype=old.dtype)
            new[:n_old] = old

            # pad using copies of last rows (repeat if DIFF > n_old)
            take = min(DIFF, n_old)
            reps = (DIFF + take - 1) // take
            pad = np.tile(old[-take:], reps)[:DIFF]
            new[n_old:] = pad

            # replace EVENTS HDU
            ev_hdu = fits.BinTableHDU(data=new, header=hdul['EVENTS'].header, name='EVENTS')
            hdul[hdul.index_of('EVENTS')] = ev_hdu

            # Now we can overwrite: use all previous ids plus the appended rows
            appended_ids = np.arange(n_old, n_new, dtype=int)
            ids_to_replace2 = np.hstack((ids_to_replace2, appended_ids))

            assert ids_to_replace2.size >= N_sim

            for fn in fi_up:
                hdul['EVENTS'].data[fn][ids_to_replace2[:N_sim]] = np.hstack((data_C[fn], data_A[fn], data_B[fn]))
            fn = 'SIGNAL'
            hdul['EVENTS'].data['PI'][ids_to_replace2] = 1000. * np.hstack((data_C[fn], data_A[fn], data_B[fn]))

        hdul.writeto(path_2_event_file, overwrite=True)
        print('\nTile {0} - File written to:\n {1}'.format(str_field, path_2_event_file))

    data_A.write(path_2_simeventAGN_file, overwrite=True)
    data_C.write(path_2_simeventCLU_file, overwrite=True)
    data_B.write(path_2_simeventBKG_file, overwrite=True)
    print('\nTile {2} - {1} AGNs written to {0}'.format(path_2_simeventAGN_file, len(data_A), str_field))
    print('\nTile {2} - {1} CLUs written to {0}'.format(path_2_simeventCLU_file, len(data_C), str_field))
    print('\nTile {2} - {1} BKGs written to {0}'.format(path_2_simeventBKG_file, len(data_B), str_field))
    N_sim = len(data_A) + len(data_C) + len(data_B)
    print('\nTile {1} - AGNs fraction: {0}'.format(np.round(len(data_A) / N_sim, 4), str_field))
    print('\nTile {1} - BKGs fraction: {0}'.format(np.round(len(data_B) / N_sim, 4), str_field))
    print('\nTile {1} - CLUs fraction: {0}'.format(np.round(len(data_C) / N_sim, 4), str_field))
    t1 = time.time()
    good.append(str_field)
    print('\nTile {0} - It took {1} sec in total.'.format(str_field, t1-t0))

print('{0} tiles out of {1} processed, {2} fails, total {3}'.format(len(good), len(sky_map_hdu[(sky_map_hdu['OWNER'] == 2) | (sky_map_hdu['OWNER'] == 0)]), len(fails), len(good)+len(fails)))