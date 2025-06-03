import sys, os, glob
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from astropy.cosmology import FlatLambdaCDM
'''
Creates key figures to analyze input-output of eROSITA simulations
Takes the desired experiment and figures directory as input, for example:

python summary_plotting_RS.py GE_e4_merge_AGNseed001_SimBKG_CLUseed001 /home/idies/workspace/Storage/rseppi/persistent/figures

R. Seppi (20.05.2025)
'''

nl = lambda selection: len(selection.nonzero()[0])
deg_to_rad = np.pi/180.
cosmo = FlatLambdaCDM(H0=67.77, Om0=0.29)

basedir = '/home/idies/workspace/erosim/Uchuu/LCerass'
#basefig = '/home/idies/workspace/Storage/rseppi/persistent/figures'

subdir = sys.argv[1]
basefig = sys.argv[2]
print('Directory:', subdir)
print('Saving figures in:', basefig)
p_2_clu_matched = os.path.join(basedir, 'SummaryFiles', 'LC_eRASS_Clusters_Input_'+subdir+'.fits')
p_2_eSASS_matched = os.path.join(basedir, 'SummaryFiles', 'LC_eRASS_eSASS_Output_'+subdir+'.fits')


#Input -> Output
#For selection function
#Which real clusters are detected?
clu = Table.read(p_2_clu_matched, memmap=True)
clu = clu[clu['in_good_region']]

def make_plot(qty, edges, logscale=True):
    if qty=='redshift_S':
        sele = clu['CLUSTER_FX_soft_OBS_R500c_nHattenuated']>-12.5
        Nall = np.histogram(clu[sele][qty], bins=edges)[0]
        seldet = sele & (clu['DET_LIKE_0']>0)
        Npoint = np.histogram(clu[qty][seldet], bins=edges)[0]
        selext = sele & (clu['EXT_LIKE']>0)
        Next = np.histogram(clu[qty][selext], bins=edges)[0]
    else:
        Nall = np.histogram(clu[qty], bins=edges)[0]
        seldet = clu['DET_LIKE_0']>0
        Npoint = np.histogram(clu[qty][seldet], bins=edges)[0]
        selext = clu['EXT_LIKE']>0
        Next = np.histogram(clu[qty][selext], bins=edges)[0]

    bins = (edges[1:]+edges[:-1])/2.
    plt.figure(figsize=(10,7))
    plt.plot(bins, Nall/np.max(Nall), c='black', lw=3, ls='--', label='TOT:%d, normed'%Nall.sum())
    plt.plot(bins, Npoint/Nall, label='DET', lw=6)
    plt.plot(bins, Next/Nall, label='DET EXT', lw=6)

    plt.legend(fontsize=20)
    if logscale:
        plt.xscale('log')
    if qty=='redshift_S':
        plt.ylabel('Completeness Fx>-12.5', fontsize=24)
    else:
        plt.ylabel('Completeness', fontsize=24)
    plt.xlabel(qty, fontsize=24)
    plt.tick_params(labelsize=22, direction='in', which='both')
    plt.tight_layout()
    p2out = os.path.join(basefig, 'Pdet_'+subdir+'_'+qty+'.png')
    plt.savefig(p2out)
    print('Written', p2out)

Flux_edges = np.linspace(-15.5, -10.5, 25)
make_plot(qty='CLUSTER_FX_soft_OBS_R500c_nHattenuated', edges=Flux_edges, logscale=False)

Mass_edges = np.geomspace(1e13, 2e15, 25)
make_plot(qty='M500c', edges=Mass_edges)

z_edges = np.geomspace(0.02, 1.5, 20)
make_plot(qty='redshift_S', edges=z_edges)

#flux in different exposure bins
qty = 'CLUSTER_FX_soft_OBS_R500c_nHattenuated'
edges = np.linspace(-15.8, -10.2, 10)
bins = (edges[1:] + edges[:-1]) / 2.
TEXPs = [300, 400, 500]

plt.figure(figsize=(10, 7))
for jj in range(len(TEXPs)+1):
    if jj==0:
        sel = clu['ExpMap_eSASS']<TEXPs[0]
        label = 'TEXP<%d'%TEXPs[0]
    elif jj==len(TEXPs):
        sel = clu['ExpMap_eSASS']>=TEXPs[-1]
        label = 'TEXP>=%d'%TEXPs[-1]
    else:
        sel = (clu['ExpMap_eSASS']>=TEXPs[jj-1])&(clu['ExpMap_eSASS']<TEXPs[jj])
        label = '%d<=TEXP<%d'%(TEXPs[jj-1],TEXPs[jj])

    Nall = np.histogram(clu[qty][sel], bins=edges)[0]
    seldet = (sel) & (clu['DET_LIKE_0'] > 0)
    Npoint = np.histogram(clu[qty][seldet], bins=edges)[0]
    selext = (sel) & (clu['EXT_LIKE'] > 0)
    Next = np.histogram(clu[qty][selext], bins=edges)[0]

    plt.plot(bins, Npoint / Nall, lw=4, ls='--', c='C'+str(jj))
    plt.plot(bins, Next / Nall, label=label, lw=6, c='C'+str(jj))

custom_lines = [Line2D([0], [0], color='black', lw=7, ls='solid'),
                Line2D([0], [0], color='black', lw=7, ls='dashed')]

legend1 = plt.legend(custom_lines, ['EXT', 'All'], fontsize=22)
plt.legend(fontsize=20)
plt.gca().add_artist(legend1)
plt.ylabel('Completeness', fontsize=24)
plt.xlabel(qty, fontsize=24)
plt.tick_params(labelsize=22, direction='in', which='both')
plt.tight_layout()
p2out = os.path.join(basefig, 'Pdet_' + subdir + '_' + qty + '_TEXP.png')
plt.savefig(p2out)
print('Written', p2out)


#Output -> Input
#For contamination
#What's in the eSASS catalogue?
cat = Table.read(p_2_eSASS_matched, memmap=True)
cat = cat[cat['in_good_region']]

Ldet_cut = np.arange(5, 51)

ALL = np.zeros_like(Ldet_cut)
PNT = np.zeros_like(Ldet_cut)
PNT2 = np.zeros_like(Ldet_cut)
EXT = np.zeros_like(Ldet_cut)
EXT2 = np.zeros_like(Ldet_cut)
BKG = np.zeros_like(Ldet_cut)

for jj,Ldet in enumerate(Ldet_cut):
    sel = cat['DET_LIKE_0']>=Ldet
    ALL[jj] = nl(sel)
    PNT[jj] = nl( np.logical_and(sel, cat['class_PNT']) )
    PNT2[jj] = nl( np.logical_and(sel, cat['class_PNT2']) )
    EXT[jj] = nl( np.logical_and(sel, cat['class_EXT']) )
    EXT2[jj] = nl( np.logical_and(sel, cat['class_EXT2']) )
    BKG[jj] = nl( np.logical_and(sel, cat['class_BKG']) )

plt.figure(figsize=(10,7))

plt.plot(Ldet_cut, PNT/ALL, lw=6, label='PNT')
plt.plot(Ldet_cut, EXT/ALL,  lw=6, label='EXT')
plt.plot(Ldet_cut, PNT2/ALL, lw=6, label='PNT2')
plt.plot(Ldet_cut, EXT2/ALL,  lw=6, label='EXT2')
plt.plot(Ldet_cut, BKG/ALL, lw=6, label='BKG')

plt.ylabel('Fraction', fontsize=24)
plt.xlabel('> DET_LIKE', fontsize=24)
plt.yscale('log')
plt.legend(fontsize=20)
plt.tick_params(labelsize=22, direction='in', which='both')
plt.tight_layout()
p2out = os.path.join(basefig, 'eSASS_' + subdir + '_Ldet.png')
plt.savefig(p2out)
print('Written', p2out)



Ldet_cut = np.arange(5, 51)

ALL = np.zeros_like(Ldet_cut)
PNT = np.zeros_like(Ldet_cut)
PNT2 = np.zeros_like(Ldet_cut)
EXT = np.zeros_like(Ldet_cut)
EXT2 = np.zeros_like(Ldet_cut)
BKG = np.zeros_like(Ldet_cut)

for jj,Ldet in enumerate(Ldet_cut):
    sel = (cat['DET_LIKE_0']>=Ldet) & (cat['EXT_LIKE']>0)
    ALL[jj] = nl(sel)
    PNT[jj] = nl( np.logical_and(sel, cat['class_PNT']) )
    PNT2[jj] = nl( np.logical_and(sel, cat['class_PNT2']) )
    EXT[jj] = nl( np.logical_and(sel, cat['class_EXT']) )
    EXT2[jj] = nl( np.logical_and(sel, cat['class_EXT2']) )
    BKG[jj] = nl( np.logical_and(sel, cat['class_BKG']) )

plt.figure(figsize=(10,7))

plt.plot(Ldet_cut, PNT/ALL, lw=6, label='PNT')
plt.plot(Ldet_cut, EXT/ALL,  lw=6, label='EXT')
plt.plot(Ldet_cut, PNT2/ALL, lw=6, label='PNT2')
plt.plot(Ldet_cut, EXT2/ALL,  lw=6, label='EXT2')
plt.plot(Ldet_cut, BKG/ALL, lw=6, label='BKG')

plt.ylabel('Fraction EXTENDED', fontsize=24)
plt.xlabel('> DET_LIKE', fontsize=24)
plt.yscale('log')
plt.legend(fontsize=20)
plt.tick_params(labelsize=22, direction='in', which='both')
plt.tight_layout()
p2out = os.path.join(basefig, 'eSASS_' + subdir + '_Ldet_EXT.png')
plt.savefig(p2out)
print('Written', p2out)



Lext_cut = np.arange(5, 51)

ALL = np.zeros_like(Ldet_cut)
PNT = np.zeros_like(Ldet_cut)
PNT2 = np.zeros_like(Ldet_cut)
EXT = np.zeros_like(Ldet_cut)
EXT2 = np.zeros_like(Ldet_cut)
BKG = np.zeros_like(Ldet_cut)

for jj,Lext in enumerate(Lext_cut):
    sel = cat['EXT_LIKE']>=Lext
    ALL[jj] = nl(sel)
    PNT[jj] = nl( np.logical_and(sel, cat['class_PNT']) )
    PNT2[jj] = nl( np.logical_and(sel, cat['class_PNT2']) )
    EXT[jj] = nl( np.logical_and(sel, cat['class_EXT']) )
    EXT2[jj] = nl( np.logical_and(sel, cat['class_EXT2']) )
    BKG[jj] = nl( np.logical_and(sel, cat['class_BKG']) )

plt.figure(figsize=(10,7))

plt.plot(Lext_cut, PNT/ALL, lw=6, label='PNT')
plt.plot(Lext_cut, EXT/ALL,  lw=6, label='EXT')
plt.plot(Lext_cut, PNT2/ALL, lw=6, label='PNT2')
plt.plot(Lext_cut, EXT2/ALL,  lw=6, label='EXT2')
plt.plot(Lext_cut, BKG/ALL, lw=6, label='BKG')

plt.ylabel('Fraction', fontsize=24)
plt.xlabel('> EXT_LIKE', fontsize=24)
plt.yscale('log')
plt.legend(fontsize=20)
plt.tick_params(labelsize=22, direction='in', which='both')
plt.tight_layout()
p2out = os.path.join(basefig, 'eSASS_' + subdir + '_Lext.png')
plt.savefig(p2out)
print('Written', p2out)
