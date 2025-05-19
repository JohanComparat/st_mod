import sys, os, glob
import numpy as np
from astropy.table import Table,vstack
from astropy.io import fits
from sklearn.neighbors import BallTree
from astropy.io import ascii
from scipy import stats

MinTotalCts = 2.
poiss_ppf = 0.8
MinAperCts = 1  # There has to be at least one to know the SimputID

basedir = '/home/idies/workspace/erosim/Uchuu/LCerass'

subdir = sys.argv[1]
tile = sys.argv[2]

print('Directory:', subdir)
print('Tile:', tile)

if not os.path.isdir(os.path.join(basedir, tile, subdir)):
    print(os.path.join(basedir, tile, subdir), 'does not exist!')
    sys.exit()


# Paths
p_2_esass = os.path.join(basedir, tile, subdir, 'eSASS', tile + '_024_Sc1Cat.fits')
p_2_evt_clu  = os.path.join(basedir, tile, subdir, 'simCLUevt_'+tile+'.fits')
p_2_evt_agn  = os.path.join(basedir, tile, subdir, 'simAGNevt_'+tile+'.fits')
p_2_evt_sta  = os.path.join(basedir, tile, subdir, 'simSTAevt_'+tile+'.fits')
p_2_evt_bkg  = os.path.join(basedir, tile, subdir, 'simBKGevt_'+tile+'.fits')

p_2_clu      = os.path.join(basedir, tile, 'Xgas_bHS0.8_simput_final.fits')
p_2_agn      = os.path.join(basedir, tile, 'AGN_list_sigma_0.8_fsat_8.0_simput.fits')
p_2_star     = os.path.join(basedir, tile, 'SIMPUT_STA.fits')

# Read simput source tables
Simput_IDs, Simput_RA, Simput_DEC, Simput_FLUX = [], [], [], []

if os.path.isfile(p_2_clu):
    clu = Table.read(p_2_clu, memmap=True)
    Simput_IDs.append(clu['SRC_ID'].astype('int32'))
    Simput_RA.append(clu['RA'])
    Simput_DEC.append(clu['DEC'])
    Simput_FLUX.append(clu['FLUX'])

if os.path.isfile(p_2_agn):
    agn = Table.read(p_2_agn, memmap=True)
    Simput_IDs.append(agn['SRC_ID'])
    Simput_RA.append(agn['RA'])
    Simput_DEC.append(agn['DEC'])
    Simput_FLUX.append(agn['FLUX'])

if os.path.isfile(p_2_star):
    star = Table.read(p_2_star, memmap=True)
    Simput_IDs.append(star['SRC_ID'])
    Simput_RA.append(star['RA'])
    Simput_DEC.append(star['DEC'])
    Simput_FLUX.append(star['FLUX'])

# Stack source properties
Simput_SrcID = np.hstack(Simput_IDs)
Simput_RA    = np.hstack(Simput_RA)
Simput_DEC   = np.hstack(Simput_DEC)
Simput_FLUX  = np.hstack(Simput_FLUX)

# Read event files
Evt_Src_list = []

if os.path.isfile(p_2_evt_clu):
    Evt_Clu = Table.read(p_2_evt_clu)
    Evt_Src_list.append(Evt_Clu)

if os.path.isfile(p_2_evt_agn):
    Evt_Agn = Table.read(p_2_evt_agn)
    Evt_Src_list.append(Evt_Agn)

if os.path.isfile(p_2_evt_sta):
    Evt_Sta = Table.read(p_2_evt_sta)
    Evt_Src_list.append(Evt_Sta)

# Combine and filter source events
if Evt_Src_list:
    Evt_Src = vstack(Evt_Src_list)
    sel_src = (Evt_Src['SIGNAL'] > 0.2) & (Evt_Src['SIGNAL'] < 2.3)
    Evt_Src = Evt_Src[sel_src]
else:
    Evt_Src = None

# Read and filter background
Evt_Bkg = None
if os.path.isfile(p_2_evt_bkg):
    Evt_Bkg = Table.read(p_2_evt_bkg)
    sel_bkg = (Evt_Bkg['ENERGY'] / 1000. > 0.2) & (Evt_Bkg['ENERGY'] / 1000. < 2.3)
    Evt_Bkg = Evt_Bkg[sel_bkg]

# Combine all events and assign SRC_ID = -1 to background
if Evt_Src is not None and Evt_Bkg is not None:
    event_energy = np.hstack((Evt_Src['SIGNAL'], Evt_Bkg['ENERGY'] / 1000.))
    event_src_id = np.hstack((Evt_Src['SRC_ID'].T[0], -1 * np.ones(len(Evt_Bkg))))
    event_ra     = np.hstack((Evt_Src['RA'], Evt_Bkg['RA']))
    event_dec    = np.hstack((Evt_Src['DEC'], Evt_Bkg['DEC']))
elif Evt_Src is not None:
    event_energy = Evt_Src['SIGNAL']
    event_src_id = Evt_Src['SRC_ID'].T[0]
    event_ra     = Evt_Src['RA']
    event_dec    = Evt_Src['DEC']
elif Evt_Bkg is not None:
    event_energy = Evt_Bkg['ENERGY'] / 1000.
    event_src_id = -1 * np.ones(len(Evt_Bkg))
    event_ra     = Evt_Bkg['RA']
    event_dec    = Evt_Bkg['DEC']
else:
    raise RuntimeError("No source or background event files found.")

#p_2_evt_sta = os.path.join(basedir, tile, 'simSTAevt_'+tile+'.fits')
#p_2_evt_bkg = os.path.join(basedir, tile, 'simBKGevt_'+tile+'.fits')
#p_2_evt_clu = os.path.join(basedir, tile, 'simCLUevt_'+tile+'.fits')
#p_2_evt_agn = os.path.join(basedir, tile, 'simAGNevt_'+tile+'.fits')

#p_2_clu = os.path.join(basedir, tile, 'Xgas_bHS0.8_simput_final.fits')
#p_2_agn = os.path.join(basedir, tile, 'AGN_list_sigma_0.8_fsat_8.0_simput.fits')
#p_2_star = os.path.join(basedir, tile, 'SIMPUT_STA.fits')

#clu = Table.read(p_2_clu, memmap=True)
#agn = Table.read(p_2_agn, memmap=True)
#if os.path.isfile(p_2_star):
#    star = Table.read(p_2_star, memmap=True)

#clu_id, agn_id, star_id = clu['SRC_ID'], agn['SRC_ID'], star['SRC_ID']
#clu_id, agn_id = clu['SRC_ID'], agn['SRC_ID']
#Simput_SrcID = np.hstack((clu_id, agn_id))#, star_id))
#Simput_RA = np.hstack((clu['RA'], agn['RA']))#, star['RA']))
#Simput_DEC = np.hstack((clu['DEC'], agn['DEC']))#, star['DEC']))
#Simput_FLUX = np.hstack((clu['FLUX'], agn['FLUX']))#, star['FLUX']))

# all events in 0.2 to 2.3
#Evt_Clu = Table.read(p_2_evt_clu)

#Evt_Agn = Table.read(p_2_evt_agn)

#Evt_Sta = Table.read(p_2_evt_sta)

#Evt_Src = vstack((Evt_Clu, Evt_Agn, Evt_Sta))
#sel_02_23 = (Evt_Src['SIGNAL'] > 0.2) & (Evt_Src['SIGNAL'] < 2.3)
#Evt_Src = Evt_Src[sel_02_23]

#Evt_Bkg = Table.read(p_2_evt_bkg)
#sel_02_23 = (Evt_Bkg['ENERGY']/1000. > 0.2) & (Evt_Bkg['ENERGY']/1000. < 2.3)
#Evt_Bkg = Evt_Bkg[sel_02_23]

# group photons from sources and bkg, select the 0.2-2.3 keV, assign srcid=-1 to bkg
event_energy = np.hstack((Evt_Src['SIGNAL'], Evt_Bkg['ENERGY']/1000.))
event_src_id = np.hstack((Evt_Src['SRC_ID'].T[0], -1 * np.ones_like(Evt_Bkg['ENERGY']/1000.)))
event_ra = np.hstack((Evt_Src['RA'], Evt_Bkg['RA']))
event_dec = np.hstack((Evt_Src['DEC'], Evt_Bkg['DEC']))

#open source catalog
hdu = fits.open(p_2_esass)#, memmap=True)
SrcCat = Table(hdu[1].data)

#To do the matching
AperRad = 20 * np.ones(len(SrcCat))
AperRad[np.where(SrcCat['EXT_LIKE'] > 0)] = 60
FittingRad = 60 #to study blending

# find counts associated to simulation source
EvtSrcID = event_src_id[event_src_id >= 0]
EvtRA_Src = Evt_Src['RA'] * np.pi / 180
EvtDEC_Src = Evt_Src['DEC'] * np.pi / 180
event_Npix_Src = Evt_Src['NPIXELS']

EvtCoord = np.transpose([EvtDEC_Src, EvtRA_Src])
EvtTree = BallTree(EvtCoord, metric='haversine')

EvtCoordBkg = np.transpose([Evt_Bkg['DEC'] * np.pi / 180, Evt_Bkg['RA'] * np.pi / 180])
EvtTreeBkg = BallTree(EvtCoordBkg, metric='haversine')


# Method1
u, inv, photoncts = np.unique(EvtSrcID, return_inverse=True, return_counts=True)
# store number of detected counts for each source
EvtSrcTotCts = {u[i]: photoncts[i] for i in range(len(u))}  # array of: SrcID:number of photons from that source

#for each simulated source id satisfying MinTotalCts, find its number of counts
SimputSrcTotCts = np.array([EvtSrcTotCts.get(i, 0) for i in Simput_SrcID])
# assert np.all(Simput['SRC_ID']>1) #-1: particle; 1:X-ray background

if not os.path.isfile(os.path.join(basedir, tile, subdir,
                                   'simput_srccts_' + tile + '_TEST_Mincts' + str(MinTotalCts) + '_ppf_' + str(
                                           poiss_ppf) + '.fits')):
    # u,inv,photoncts=np.unique(EvtSrcID[(EvtPI>=200)&(EvtPI<500)],return_inverse=True,return_counts=True)
    u, inv, photoncts = np.unique(EvtSrcID, return_inverse=True, return_counts=True)
    EvtSrcCtsB1 = {u[i]: photoncts[i] for i in range(len(u))}
    del u, inv, photoncts
    SimputSrcCtsB1 = np.array([EvtSrcCtsB1.get(i, 0) for i in Simput_SrcID])
    ok = SimputSrcTotCts > 0
    fits.BinTableHDU.from_columns([fits.Column(name='SRC_ID', format='J', array=Simput_SrcID[ok]),
                                   fits.Column(name='CtsTot', format='J', array=SimputSrcTotCts[ok]),
                                   fits.Column(name='CtsB1', format='J', array=SimputSrcCtsB1[ok])]).writeto(
        os.path.join(basedir, tile, subdir, 'simput_srccts_' + tile + '.fits'), overwrite=True)
    print('>>', os.path.join(basedir, tile, subdir, 'simput_srccts_' + tile + '.fits'))
    #os.system('chgrp erosim ' + os.path.join(basedir, tile, 'simput_srccts_' + tile + '.fits'))
    # exit()
else:
    SimputSrcCts = fits.getdata(os.path.join(basedir, tile, subdir, 'simput_srccts_' + tile + '.fits'))

# Now study blend: analyze photons around simulated sources
# I don't remember why in Sim12 I changed from ID_blend to ID_blend2 adding a Contamination_MaxFactor upper limit. Such an upper limit is totally unresonable if the blending is defined as an one-way contamination. It only makes sense if we need to consider contamination to each other. Now I define blending in a new way, considering only one-way contamination.

# Now create a file with IDSrc1 IDSrc2 and Cts1 Cts2 as the counts coming from each source inside 60 arcsec (FittingRad)
if not os.path.isfile(os.path.join(basedir, tile, subdir, 'ID_contam_' + tile + '.fits')):
    simOK = (SimputSrcTotCts > 0)
    Simput_SrcID = Simput_SrcID[simOK]
    SimputSrcID = Simput_SrcID
    Simput_RA, Simput_DEC = Simput_RA[simOK], Simput_DEC[simOK]
    Simput_FLUX = Simput_FLUX[simOK]
    SimputSrcTotCts = SimputSrcTotCts[simOK]
    SimputCoord = np.transpose([Simput_DEC * np.pi / 180, Simput_RA * np.pi / 180])

    Contamination_MinCts = 3
    EvtInd4Sim = EvtTree.query_radius(SimputCoord, r=FittingRad * np.pi / 180 / 3600)
    blend = []
    for n in range(len(SimputSrcID)):  # Each simput source
        if len(EvtInd4Sim[n]) > 0:
            #print('n:',n)
            #print('len(EvtInd4Sim[n]):', len(EvtInd4Sim[n]))
            u, inv, photoncts = np.unique(EvtSrcID[EvtInd4Sim[n]], return_inverse=True, return_counts=True)
            if SimputSrcID[n] not in u: continue
            if 0:
                eventcts = np.zeros_like(u)
                for i in range(len(EvtInd4Sim[n])): eventcts[inv[i]] += event_Npix_Src[EvtInd4Sim[n]][i]
            else:
                eventcts = photoncts
            if len(u) <= 1: continue
            thisone = u == SimputSrcID[n]
            assert np.sum(thisone) == 1
            ctsthisone = eventcts[thisone][0]
            eventcts = eventcts[~thisone]
            u = u[~thisone]
            second = np.argmax(eventcts)
            if eventcts[second] >= max(Contamination_MinCts, np.sqrt(ctsthisone)):
                blend.append([SimputSrcID[n], Simput_FLUX[n], u[second], ctsthisone, eventcts[second]])
            else:
                blend.append([SimputSrcID[n], Simput_FLUX[n], -1 * u[second], ctsthisone, eventcts[second]])
        # else: blend.append([SimputSrcID[n],-99,0,0])
    print('>>', os.path.join(basedir, tile, subdir, 'ID_contam_' + tile + '.csv'))
    np.savetxt(os.path.join(basedir, tile, subdir, 'ID_contam_' + tile + '.csv'), np.array(blend), fmt='%d %.4e %d %d %d',
               delimiter=' ', header='ID_1 FLUX ID_2 counts_1 counts_2', comments='#')
    # exit()
    a = ascii.read(os.path.join(basedir, tile, subdir, 'ID_contam_' + tile + '.csv'), format='fast_csv', delimiter=' ')
    if a.columns[0].name[0] == '#': a.columns[0].name = a.columns[0].name[1:]
    a.write(os.path.join(basedir, tile, subdir, 'ID_contam_' + tile + '.fits'), overwrite=True)
else:
    d = Table.read(os.path.join(basedir, tile, subdir, 'ID_contam_' + tile + '.fits'))
    ID_contam = d['ID_1'][d['ID_2'] > 0]

simOK = (SimputSrcTotCts >= MinTotalCts)
# Simput_FLUX = Simput_FLUX[simOK]
SimputSrcTotCts = SimputSrcTotCts[simOK]
SimputCoord = np.transpose([Simput_DEC[simOK] * np.pi / 180, Simput_RA[simOK] * np.pi / 180])
SimputSrcID = Simput_SrcID[simOK]
# SimputTree=BallTree(SimputCoord,metric='haversine')

SrcCoord = np.transpose([SrcCat['DEC'] * np.pi / 180, SrcCat['RA'] * np.pi / 180])
EvtInd4Cat = EvtTree.query_radius(SrcCoord, r=AperRad * np.pi / 180 / 3600)
EvtInd4Cat_bkg = EvtTreeBkg.query_radius(SrcCoord,
                                         r=AperRad * np.pi / 180 / 3600)  # array of arrays. Index of bkg photons around each source
SrcCat_N_Bkg = [len(arr) for arr in EvtInd4Cat_bkg]  # N of bkg photons around each detection

if not os.path.isfile(
        os.path.join(basedir, tile, subdir, f'Sc1_' + tile + '_IDMatch_Any_Tot' + str(MinTotalCts) + '.fits.cp')):
    SrcID = np.zeros(len(SrcCoord), dtype=np.int64) - 1
    SrcApeCts = np.zeros(len(SrcCoord), dtype=np.int32)  # Aperture counts
    SrcID2 = np.zeros(len(SrcCoord), dtype=np.int64) - 1
    SrcApeCts2 = np.zeros(len(SrcCoord), dtype=np.int32)
    ftmp = open('tmp.dat', 'w')
    for n in range(len(SrcCat)):  # Each detected source
        if len(EvtInd4Cat[n]) > 0:  # if there are photons around that entry:
            print(SrcCat['ID_SRC'][n], np.unique(EvtSrcID[EvtInd4Cat[n]]),
                  file=ftmp)  # esassid, id of source photons <AperRad from detection
            ok = np.array([EvtSrcTotCts.get(i, 0) for i in EvtSrcID[EvtInd4Cat[
                n]]]) >= MinTotalCts  # select photons from sources that have at least 3 photons observed
            u, inv, photoncts = np.unique(EvtSrcID[EvtInd4Cat[n][ok]], return_inverse=True, return_counts=True)
            uc = np.array([[u[i], photoncts[i]] for i in range(len(u)) if photoncts[
                i] >= MinAperCts])  # ID of simput source related to this entry and cts from it
            apebkg = SrcCat_N_Bkg[n]
            if len(uc) == 0: print('no', [[u[i], photoncts[i]] for i in range(len(u))],
                                   file=ftmp)  # no simput corresponding to detection
            if len(uc) == 1:  # one simput correspponding to detection
                if uc[0, 1] + apebkg > stats.poisson.ppf(0.97725, apebkg): SrcID[n], SrcApeCts[n] = uc[
                    0]  # if the simput is significant over the bkg, associate it to the esass detection with the simput id and its counts                
                # else: SrcID2[n],SrcApeCts2[n]=uc[0] #NO, don't keep it
                # print('1',SrcID[n],SrcApeCts[n],file=ftmp)
            elif len(uc) > 1:  # more than one simput corresponding to detection
                print('>1', uc, file=ftmp)
                uc = uc[np.argsort(uc[:, 1])]  # sort the simput corresponding to detection by number of counts
                if uc[-1, 1] + apebkg > stats.poisson.ppf(0.97725,
                                                          apebkg):  # select the one that contributes the largest number of events
                    # if uc[-1][1]==1: print(uc[-1][1],apebkg,stats.poisson.ppf(0.97725,apebkg))
                    SrcID[n], SrcApeCts[n] = uc[-1]
                    if uc[-2, 1] + apebkg > stats.poisson.ppf(0.97725, apebkg): SrcID2[n], SrcApeCts2[n] = uc[
                        -2]  # also save the one with second highest counts
        else:  # these are clearly bkg fluctuations
            print(SrcCat['ID_SRC'][n], np.unique(EvtSrcID[EvtInd4Cat[n]]), len(EvtInd4Cat[n]), n,
                  SrcCat['RA'][n], SrcCat['DEC'][n], SrcCoord[n], file=ftmp)

    ftmp.close()

    ok = SrcID >= 0  # esass detections with a simput counterpart
    SrcID[~ok] = -99
    SrcApeCts[~ok] = -99
    ok = SrcID2 >= 0  # esass detections with secondary simput counterpart
    SrcID2[~ok] = -99
    SrcApeCts2[~ok] = -99
    outdata = np.transpose(
        [SrcCat['ID_SRC'], SrcCat['RA'], SrcCat['DEC'], SrcCat['DET_LIKE_0'], SrcCat['EXT'], SrcCat['EXT_LIKE'],
         SrcCat['ML_CTS_0'], SrcCat['ML_RATE_0'], SrcCat['ML_FLUX_0'], SrcID, SrcApeCts, SrcID2, SrcApeCts2])
    print('>>', os.path.join(basedir, tile, subdir, f'Sc1_' + tile + '_IDMatch_Any_Tot' + str(MinTotalCts) + '.csv'))
    np.savetxt(os.path.join(basedir, tile, subdir, f'Sc1_' + tile + '_IDMatch_Any_Tot' + str(MinTotalCts) + '.csv'), outdata,
               fmt='%d %.5f %.5f %.3f %.3f %.3f %.4f %.7f %.7g %d %d %d %d',
               delimiter=' ',
               header='ID_cat RA DEC DET_LIKE_0 EXT EXT_LIKE ML_CTS_0 ML_RATE_0 ML_FLUX_0 ID_simput aperture_counts ID_simput_2 aperture_counts_2',
               comments='#')
    a = ascii.read(os.path.join(basedir, tile, subdir, f'Sc1_' + tile + '_IDMatch_Any_Tot' + str(MinTotalCts) + '.csv'),
                   format='fast_csv', delimiter=' ')
    if a.columns[0].name[0] == '#': a.columns[0].name = a.columns[0].name[1:]
    a.write(os.path.join(basedir, tile, subdir, f'Sc1_' + tile + '_IDMatch_Any_Tot' + str(MinTotalCts) + '.fits'),
            overwrite=True)

if not os.path.isfile(
        os.path.join(basedir, tile, subdir, f'Sc1_' + tile + '_IDMatch_Uniq_Tot' + str(MinTotalCts) + '.fits.cp')):
    selected = np.zeros(len(SrcCoord), dtype=bool)
    s = np.argsort(SrcApeCts)[::-1]  # sort by counts the simput associated to esass detection
    _uID, ind = np.unique(SrcID[s],
                          return_index=True)  # identify esass detection with ONE single simput counterpart, save the srcid and the first time they appear (ind)
    selected[s[ind]] = True  # save the ones that appear only once
    selected[SrcApeCts <= 0] = False  # save the esass detections that have no simput counterpart

    # find non-selected with detected counts
    for n in np.where(~selected & (SrcApeCts2 > 0))[
        0]:  # consider multiple detections pointing to the same simput and with a secondary match
        assert SrcApeCts[n] > 0
        assert SrcID[n] in SrcID[selected]
        # if np.sum((SrcID==SrcID[n])&selected)>1: print(np.where(SrcID==SrcID[n])) print(selected[SrcID==SrcID[n]]) print(SrcID[n],np.sum((SrcID==SrcID[n])&selected), SrcApeCts[(SrcID==SrcID[n])&selected])
        assert SrcApeCts[(SrcID == SrcID[n]) & selected] >= SrcApeCts[
            n]  # if the id_any is not unique and id_any2 has not been assigned already, take id_any2 as unique counterpart
        if SrcID2[n] not in SrcID[selected]:
            SrcID[n], SrcApeCts[n] = SrcID2[n], SrcApeCts2[n]
            SrcID2[n], SrcApeCts2[n] = 0, 0
            selected[n] = True

    def matchlines(list1, list2):
        list2 = dict((r, i) for i, r in enumerate(list2))
        return np.array([list2.get(x, -1) for x in list1])

    ind = matchlines(SrcID[selected], SimputSrcID)
    ok = (SrcID[selected] != SimputSrcID[ind])
    assert all(SrcID[selected] == SimputSrcID[ind])
    SrcTotCts = np.zeros(len(SrcCoord), dtype=np.int32)
    SrcTotCts[selected] = SimputSrcTotCts[ind]

    SrcID[~selected] = -99
    SrcApeCts[~selected] = -99
    outdata = np.transpose(
        [SrcCat['ID_SRC'], SrcCat['RA'], SrcCat['DEC'], SrcCat['DET_LIKE_0'], SrcCat['EXT'], SrcCat['EXT_LIKE'],
         SrcCat['ML_CTS_0'], SrcCat['ML_RATE_0'], SrcCat['ML_FLUX_0'], SrcID, SrcApeCts, SrcTotCts])
    print('>>', os.path.join(basedir, tile, subdir, f'Sc1_' + tile + '_IDMatch_Uniq_Tot' + str(MinTotalCts) + '.csv'))
    np.savetxt(os.path.join(basedir, tile, subdir, f'Sc1_' + tile + '_IDMatch_Uniq_Tot' + str(MinTotalCts) + '.csv'), outdata,
               fmt='%d %.5f %.5f %.3f %.3f %.3f %.4f %.7f %.7g %d %d %d',
               delimiter=' ',
               header='ID_cat RA DEC DET_LIKE_0 EXT EXT_LIKE ML_CTS_0 ML_RATE_0 ML_FLUX_0 ID_simput aperture_counts total_counts',
               comments='#')
    a = ascii.read(os.path.join(basedir, tile, subdir, f'Sc1_' + tile + '_IDMatch_Uniq_Tot' + str(MinTotalCts) + '.csv'),
                   format='fast_csv', delimiter=' ')
    if a.columns[0].name[0] == '#': a.columns[0].name = a.columns[0].name[1:]
    a.write(os.path.join(basedir, tile, subdir, f'Sc1_' + tile + '_IDMatch_Uniq_Tot' + str(MinTotalCts) + '.fits'),
            overwrite=True)

os.system('rm ' + basedir+'/' + tile + '/' + subdir + '/ID_*csv')  # remove the csv files
os.system('rm ' + basedir+'/' + tile + '/' + subdir + '/Sc1_*csv')
