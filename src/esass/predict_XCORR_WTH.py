# mamba activate /home/idies/workspace/erosim/software/JC/env/clustering
import os, sys
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ---------------------------------------------------------------------------
# Figure style. Colour encodes the X-ray component, line style the galaxy
# sample, line width / zorder the emphasis. Palette: Okabe-Ito (CVD-safe).
# ---------------------------------------------------------------------------
C_AGN = '#E69F00'   # orange - AGN emission
C_GAS = '#984EA3'   # purple - hot gas (CGM/ICM)
C_BKG = '#009E73'   # green  - X-ray background
C_TOT = '#000000'   # black  - full model, measured on the joint catalogue
C_SUM = '#666666'   # grey   - sum of the individual components
C_CMP = '#0072B2'   # blue   - second model in the comparison corner plot

COMP_COLOR = {'AGN': C_AGN, 'GAS': C_GAS, 'BKG': C_BKG}
SAMPLE_LS  = {'all': 'solid', 'cen': 'dashed', 'sat': 'dotted'}
SAMPLE_LW  = {'all': 2.5,     'cen': 1.6,      'sat': 1.6}

def sty(comp, sample='all', **kw):
	"""matplotlib kwargs for one (component x galaxy-sample) curve."""
	out = dict(color=COMP_COLOR[comp], ls=SAMPLE_LS[sample],
			   lw=SAMPLE_LW[sample], zorder=3)
	out.update(kw)
	return out

# totals and data: neutral colours, separated by dash pattern
STY_GAB = dict(color=C_TOT, ls=(0, (6, 2)),           lw=2.2, zorder=6)
STY_GA  = dict(color=C_TOT, ls=(0, (4, 1.5, 1, 1.5)), lw=1.6, zorder=5)
STY_SUM = dict(color=C_SUM, ls='solid',               lw=2.0, zorder=2)
STY_OBS = dict(color='k', marker='*', ls='none', markersize=7, zorder=10)

L_OBS = 'Obs, C25'
L_AGN, L_AGN_C, L_AGN_S = r'Gal$\times$AGN', r'Cen$\times$AGN', r'Sat$\times$AGN'
L_GAS, L_GAS_C, L_GAS_S = r'Gal$\times$Gas', r'Cen$\times$Gas', r'Sat$\times$Gas'
L_BKG, L_BKG_C, L_BKG_S = r'Gal$\times$Bkg', r'Cen$\times$Bkg', r'Sat$\times$Bkg'
L_GAB, L_GA, L_SUM = r'Gal$\times$GAB', r'Gal$\times$GA', 'Sum'

LEGEND_GROUPS = [
	[L_OBS,   L_GAB,   L_GA,    L_SUM],   # data + totals
	[L_GAS,   L_GAS_C, L_GAS_S],		  # gas family
	[L_AGN,   L_AGN_C, L_AGN_S],          # AGN family
	[ L_BKG, L_BKG_C, L_BKG_S],   # gas family + background
]

def grouped_legend(title, groups=LEGEND_GROUPS, **kw):
	"""Legend with one physical family per column.

	Labels not present on the current axes are dropped, so the same grouping
	serves every figure. Blank handles pad short columns, because matplotlib
	fills legend columns top-to-bottom.
	"""
	handles, labels = plt.gca().get_legend_handles_labels()
	by_label = dict(zip(labels, handles))
	cols = [[g for g in grp if g in by_label] for grp in groups]
	cols = [c for c in cols if c]
	nrow = max(len(c) for c in cols)
	blank = Line2D([], [], ls='none')
	H, L = [], []
	for c in cols:
		H += [by_label[g] for g in c] + [blank] * (nrow - len(c))
		L += list(c)                  + ['']    * (nrow - len(c))
	opts = dict(loc='lower left', ncol=len(cols), fontsize=9, title=title,
				framealpha=0.9, handlelength=2.6, columnspacing=1.0,
				labelspacing=0.3)
	opts.update(kw)
	return plt.legend(H, L, **opts)

import glob
import numpy as np
from astropy.table import Table, Column, vstack, hstack
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from scipy.interpolate import interp1d
import emcee
import corner
cosmoUNIT = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
h = 0.6774
L_box = 1000.0 / h
cosmo = cosmoUNIT
sky_map_hdu = Table.read(os.path.join(os.environ['GIT_STMOD_DATA'], 'data/models/eROSITA', 'SKYMAPS.fits') )
is_eroDE = ( (sky_map_hdu['OWNER']==2)|(sky_map_hdu['OWNER']==0) ) & ( abs(sky_map_hdu['GLAT_CEN']) > 20 ) #& ( sky_map_hdu['DE_CEN'] <= 32 )
SRVMAP_exGAL_eroDE = sky_map_hdu['SRVMAP'][is_eroDE]
benchmark_data_dir = os.path.join( os.environ['GIT_STMOD_DATA'], 'data/benchmark/xcorr_comparat_2025' )
benchmark_dir = os.path.join( os.environ['GIT_STMOD_DATA'], 'data/benchmark/xcorr_comparat_2025' , 'XCORR_SHELLzlt04')
benchmark_dir = os.path.join( os.environ['GIT_STMOD_DATA'], 'data/benchmark/xcorr_comparat_2025' , 'XCORR')

dir_fig = os.path.join( os.environ['GIT_STMOD_DATA'], 'data/validation/validation_GasGal/XCORR_SHELL' )
dir_fig = os.path.join( os.environ['GIT_STMOD_DATA'], 'data/validation/validation_GasGal/XCORR' )
os.system('mkdir -p '+dir_fig)

agn_seed = '1' # sys.argv[1] # 1
clu_seed = '19' # sys.argv[2] # 1
LC_dir = 'LCerass'

# dir_2pcf = os.path.join(os.environ['UCHUU'], LC_dir, '??????',
#                          'GE_e4_merge_AGNseed' + agn_seed.zfill(3) + '_SimBKG_CLUseed' + clu_seed.zfill(3), 'XCORR')

XCORR = {}
XCORR['m10.0'] = Table.read( os.path.join(benchmark_data_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALxEVT_wtheta.fits' ) )
XCORR['m10.5'] = Table.read( os.path.join(benchmark_data_dir, 'LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_GALxEVT_wtheta.fits' ) )
XCORR['m11.0'] = Table.read( os.path.join(benchmark_data_dir, 'LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_GALxEVT_wtheta.fits' ) )

PRED = {}

PRED['m10.0_GxG'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALxEVT_wtheta_prediction_GxG.fits'))
PRED['m10.0_GxB'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALxEVT_wtheta_prediction_GxB.fits'))
PRED['m10.0_GxA'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALxEVT_wtheta_prediction_GxA.fits'))
PRED['m10.0_GxGA' ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALxEVT_wtheta_prediction_GxGA.fits'))
PRED['m10.0_GxGAB'] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALxEVT_wtheta_prediction_GxGAB.fits'))

PRED['m10.0C_GxG'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALCENxEVT_wtheta_prediction_GxG.fits'))
PRED['m10.0C_GxB'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALCENxEVT_wtheta_prediction_GxB.fits'))
PRED['m10.0C_GxA'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALCENxEVT_wtheta_prediction_GxA.fits'))
PRED['m10.0C_GxGA' ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALCENxEVT_wtheta_prediction_GxGA.fits'))
PRED['m10.0C_GxGAB'] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALCENxEVT_wtheta_prediction_GxGAB.fits'))

PRED['m10.0S_GxG'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALSATxEVT_wtheta_prediction_GxG.fits'))
PRED['m10.0S_GxB'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALSATxEVT_wtheta_prediction_GxB.fits'))
PRED['m10.0S_GxA'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALSATxEVT_wtheta_prediction_GxA.fits'))
PRED['m10.0S_GxGA' ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALSATxEVT_wtheta_prediction_GxGA.fits'))
PRED['m10.0S_GxGAB'] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238_GALSATxEVT_wtheta_prediction_GxGAB.fits'))



PRED['m10.5_GxG'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_GALxEVT_wtheta_prediction_GxG.fits'))
PRED['m10.5_GxB'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_GALxEVT_wtheta_prediction_GxB.fits'))
PRED['m10.5_GxA'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_GALxEVT_wtheta_prediction_GxA.fits'))
PRED['m10.5_GxGA' ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_GALxEVT_wtheta_prediction_GxGA.fits'))
PRED['m10.5_GxGAB'] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_GALxEVT_wtheta_prediction_GxGAB.fits'))

PRED['m10.5C_GxG'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_GALCENxEVT_wtheta_prediction_GxG.fits'))
PRED['m10.5C_GxB'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_GALCENxEVT_wtheta_prediction_GxB.fits'))
PRED['m10.5C_GxA'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_GALCENxEVT_wtheta_prediction_GxA.fits'))
PRED['m10.5C_GxGA' ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_GALCENxEVT_wtheta_prediction_GxGA.fits'))
PRED['m10.5C_GxGAB'] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_GALCENxEVT_wtheta_prediction_GxGAB.fits'))

PRED['m10.5S_GxG'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_GALSATxEVT_wtheta_prediction_GxG.fits'))
PRED['m10.5S_GxB'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_GALSATxEVT_wtheta_prediction_GxB.fits'))
PRED['m10.5S_GxA'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_GALSATxEVT_wtheta_prediction_GxA.fits'))
PRED['m10.5S_GxGA' ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_GALSATxEVT_wtheta_prediction_GxGA.fits'))
PRED['m10.5S_GxGAB'] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228_GALSATxEVT_wtheta_prediction_GxGAB.fits'))



PRED['m11.0_GxG'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_GALxEVT_wtheta_prediction_GxG.fits'))
PRED['m11.0_GxB'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_GALxEVT_wtheta_prediction_GxB.fits'))
PRED['m11.0_GxA'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_GALxEVT_wtheta_prediction_GxA.fits'))
PRED['m11.0_GxGA' ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_GALxEVT_wtheta_prediction_GxGA.fits'))
PRED['m11.0_GxGAB'] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_GALxEVT_wtheta_prediction_GxGAB.fits'))

PRED['m11.0C_GxG'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_GALCENxEVT_wtheta_prediction_GxG.fits'))
PRED['m11.0C_GxB'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_GALCENxEVT_wtheta_prediction_GxB.fits'))
PRED['m11.0C_GxA'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_GALCENxEVT_wtheta_prediction_GxA.fits'))
PRED['m11.0C_GxGA' ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_GALCENxEVT_wtheta_prediction_GxGA.fits'))
PRED['m11.0C_GxGAB'] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_GALCENxEVT_wtheta_prediction_GxGAB.fits'))

PRED['m11.0S_GxG'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_GALSATxEVT_wtheta_prediction_GxG.fits'))
PRED['m11.0S_GxB'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_GALSATxEVT_wtheta_prediction_GxB.fits'))
PRED['m11.0S_GxA'  ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_GALSATxEVT_wtheta_prediction_GxA.fits'))
PRED['m11.0S_GxGA' ] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_GALSATxEVT_wtheta_prediction_GxGA.fits'))
PRED['m11.0S_GxGAB'] = Table.read( os.path.join(benchmark_dir, 'LS10_VLIM_ANY_11.0_Mstar_12.0_0.05_z_0.35_N_1619838_GALSATxEVT_wtheta_prediction_GxGAB.fits'))



UNIT_WTH = dict(tag='', fx=1.0, fy=1.0, xlim=(0.0008, 3), ylim=(0.0001, 3),
														 xlabel=r'$\theta$ (degrees)', ylabel=r'$w(\theta)$')

_U = UNIT_WTH


def use_units(U):
	"""Select the unit system pl/plot_obs/band_fill/finish_wtheta_fig draw in."""
	global _U
	_U = U


def pl(x, y, **kw):
	"""plt.plot, with x and y converted into the current unit system."""
	return plt.plot(np.asarray(x) * _U['fx'], np.asarray(y) * _U['fy'], **kw)


def plot_obs(x_obs, y_obs, in_fit, ax=None, label=L_OBS, fx=None, fy=None):
	"""Observed points: filled where they enter the fit, hollow where they do not."""
	tgt = plt if ax is None else ax
	x_obs = np.asarray(x_obs) * (_U['fx'] if fx is None else fx)
	y_obs = np.asarray(y_obs) * (_U['fy'] if fy is None else fy)
	tgt.plot(x_obs[in_fit], y_obs[in_fit], label=label, **STY_OBS)
	tgt.plot(x_obs[~in_fit], y_obs[~in_fit], markerfacecolor='none', **STY_OBS)


def start_wtheta_fig(residual=True):
	"""Open the w(theta) panel, with a narrow residual panel below it if asked."""
	h = 6.5 if residual else 5.5
	fig = plt.figure(11, (6.5, h))
	fig.clf()
	fig.set_size_inches(6.5, h)
	if not residual:
		ax = fig.add_subplot(1, 1, 1)
		plt.sca(ax)
		return ax, None
	gs = fig.add_gridspec(2, 1, height_ratios=[3.2, 1], hspace=0.05)
	ax = fig.add_subplot(gs[0])
	ax_r = fig.add_subplot(gs[1], sharex=ax)
	plt.sca(ax)
	return ax, ax_r


def finish_wtheta_fig(p_2_2PCF_figure, str_title, ax, ax_r,
																					 x_obs, y_obs, in_fit, theta_model=None, w_model=None,
																					 w_lo=None, w_hi=None, extra_legend=None):
	"""Axes, legend and save shared by the w(theta) figures.

	When ax_r is given, the lower panel shows data/model - 1 against the summed
	model, restricted to the observed theta that the model grid actually covers.
	w_lo/w_hi, when given, shade the 68 per cent band of that same ratio. Pass
	ax_r=None (start_wtheta_fig(residual=False)) for a single-panel figure.
	"""
	x_obs = np.asarray(x_obs)
	y_obs = np.asarray(y_obs)
	if ax_r is not None:
		tm = np.asarray(theta_model)
		ok = (x_obs >= tm.min()) & (x_obs <= tm.max())
		def _ratio(w):
			return y_obs[ok] / interp1d(tm, np.asarray(w))(x_obs[ok]) - 1.0
		if w_lo is not None and w_hi is not None:
			ax_r.fill_between(x_obs[ok] * _U['fx'], _ratio(w_hi), _ratio(w_lo),
																				 color=C_SUM, alpha=0.25, lw=0, zorder=1)
		ax_r.axhline(0.0, color=C_SUM, lw=1.0, zorder=2)
		plot_obs(x_obs[ok], _ratio(w_model), np.asarray(in_fit)[ok], ax=ax_r,
									  label=None, fy=1.0)

	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlim(_U['xlim'])
	ax.set_ylim(_U['ylim'])
	ax.tick_params(labelsize=14, labelbottom=(ax_r is None))
	ax.set_ylabel(_U['ylabel'], fontsize=18)
	ax.set_title('Correlation Galaxies x 0.5-2 keV events')
	plt.sca(ax)
	leg = grouped_legend(str_title)
	if extra_legend is not None:
		# keep the component legend, then add the one naming the two models
		ax.add_artist(leg)
		ax.legend(extra_legend[0], extra_legend[1], loc='upper right',
										  fontsize=9, frameon=False)

	if ax_r is None:
		ax.set_xlabel(_U['xlabel'], fontsize=18)
	else:
		ax_r.set_xscale('log')
		ax_r.set_xlim(_U['xlim'])
		ax_r.tick_params(labelsize=12)
		ax_r.set_xlabel(_U['xlabel'], fontsize=18)
		ax_r.set_ylabel(r'$\frac{data}{model}-1$', fontsize=15)
		ax_r.grid(True, axis='y', alpha=0.25, lw=0.5)

	ax.get_figure().subplots_adjust(left=0.16, right=0.97, top=0.94,
																															  bottom=0.105 if ax_r is not None else 0.12,
																															  hspace=0.05)
	plt.savefig(p_2_2PCF_figure)#, transparent = True)
	plt.clf()
	print(p_2_2PCF_figure, 'written')



# ---------------------------------------------------------------------------
# Mixture parametrisation.
#
# A template is the cross-correlation of one galaxy sample g (centrals or
# satellites) with one X-ray event class X (AGN, hot gas, background). It enters
# any combination with the number of galaxy-photon pairs it stands for,
# (f_g N_g) x (f_X N_X): one weight rescaling the galaxy count of the sample,
# one rescaling the photon budget of the class. The model is therefore bilinear
# in the two weight vectors
#     u = (f_cen, f_sat)            galaxy samples
#     v = (f_AGN, f_Gas, f_Bkg)     event classes
# and reads w(theta) = (u . M(theta) . v) / (u . D . v), with
# M[i, j] = N_g N_X w_gX(theta) and D[i, j] = N_g N_X.
# ---------------------------------------------------------------------------
SAMPLES = ('cen', 'sat')
CLASSES = ('AGN', 'GAS', 'BKG')
TEMPLATE_TAG = {('cen', 'AGN'): 'C_GxA', ('sat', 'AGN'): 'S_GxA',
		   ('cen', 'GAS'): 'C_GxG', ('sat', 'GAS'): 'S_GxG',
		   ('cen', 'BKG'): 'C_GxB', ('sat', 'BKG'): 'S_GxB'}
L_TEMPLATE = {('cen', 'AGN'): L_AGN_C, ('sat', 'AGN'): L_AGN_S,
		 ('cen', 'GAS'): L_GAS_C, ('sat', 'GAS'): L_GAS_S,
		 ('cen', 'BKG'): L_BKG_C, ('sat', 'BKG'): L_BKG_S}

L_F_AGN, L_F_GAS, L_F_BKG = r'$f_{AGN}$', r'$f_{Gas}$', r'$f_{Bkg}$'
L_F_CEN, L_F_SAT = r'$f_{cen}$', r'$f_{sat}$'
L_F_CLASS = {'AGN': L_F_AGN, 'GAS': L_F_GAS, 'BKG': L_F_BKG}
L_F_SAMPLE = {'cen': L_F_CEN, 'sat': L_F_SAT}


def fit_labels(samples, classes):
	"""Corner-plot labels of one flavour, sample weights first."""
	return [L_F_SAMPLE[s] for s in samples] + [L_F_CLASS[c] for c in classes]


def bilinear_model(p, M, D, n_s):
	"""w(theta) = (u . M(theta) . v) / (u . D . v), from the packed vector (u, v)."""
	u = np.asarray(p[:n_s], dtype=float)
	v = np.asarray(p[n_s:], dtype=float)
	return np.einsum('i,ijk,j->k', u, M, v) / float(u @ D @ v)


# ---------------------------------------------------------------------------
# What the likelihood assumes about the uncertainty of the measurement.
#
# This file uses a flat 5 per cent per bin. predict_XCORR_WTH_err.py is the same
# analysis driven by the formal per-bin uncertainties tabulated with the
# measurement; the two are meant to be read against each other.
# ---------------------------------------------------------------------------
SIGMA_TAG   = '5pc'
SIGMA_LABEL = 'flat 5 per cent per bin'

def SIGMA_FUN(y_fit, err_fit):
	return 0.05 * np.asarray(y_fit, dtype=float)

# Filled by do_fit, one entry per (stellar mass sample, fit variant).
SUMMARY = {}


def split_rhat(chain):
	"""Split-Rhat (Gelman-Rubin) per parameter from a (n_step, n_walker, n_dim) chain.

	Each walker is cut in half and the halves are treated as independent chains,
	which catches a sampler that has not forgotten its starting point as well as
	one whose walkers disagree.
	"""
	n_step, n_walk, n_dim = chain.shape
	h = n_step // 2
	c = np.concatenate([chain[:h], chain[h:2 * h]], axis=1)   # (h, 2*n_walk, n_dim)
	m = c.shape[1]
	mean_j = c.mean(axis=0)                     # (m, n_dim)
	var_j = c.var(axis=0, ddof=1)               # (m, n_dim)
	W = var_j.mean(axis=0)
	B = h * mean_j.var(axis=0, ddof=1)
	var_hat = (h - 1) / h * W + B / h
	return np.sqrt(var_hat / np.where(W > 0, W, np.nan))


def fit_diagnostics(sampler, n_burn, y_obs, sigma, model_med, n_dim, labels):
	"""Goodness of fit at the posterior median, and sampler convergence.

	chi2 is evaluated with the SAME sigma the likelihood used, so chi2/ndof is
	directly the quantity the fit optimised; ndof = n_bins - n_sampled, the two
	derived weights costing no freedom of their own.
	"""
	resid = (y_obs - model_med) / sigma
	chi2 = float(np.sum(resid ** 2))
	ndof = int(len(y_obs) - n_dim)
	post = sampler.get_chain(discard=n_burn)                  # (n_step, n_walk, n_dim)
	tau = sampler.get_autocorr_time(quiet=True)
	n_eff = post.shape[0] * post.shape[1] / np.maximum(tau, 1.0)
	return dict(chi2=chi2, ndof=ndof, chi2_red=chi2 / ndof,
				acceptance=float(np.mean(sampler.acceptance_fraction)),
				tau=np.asarray(tau, dtype=float),
				n_eff=np.asarray(n_eff, dtype=float),
				rhat=np.asarray(split_rhat(post), dtype=float),
				n_bins=int(len(y_obs)), n_par=int(n_dim),
				labels=np.array(labels))


def run_mcmc_fit(M_fit, D, samples, classes, x_fit, y_fit, corner_path,
		n_walkers=32, n_steps=60000, n_burn=10000, seed=42, sigma=None):
	"""Sample the bilinear mixture weights with emcee.

	M_fit[i, j] holds N_g N_X w_gX already interpolated onto x_fit and D[i, j] the
	matching constant N_g N_X, for sample samples[i] and class classes[j];
	templates the flavour leaves out are zero in both. The model is invariant
	under u -> a u and, separately, under v -> b v, so only the RATIOS within each
	block are constrained and a flat box prior would leave the chain wandering
	along those two rays. Each block is therefore forced onto its own simplex:
	(n_s - 1) + (n_c - 1) directions are sampled and the two remaining weights are
	derived. Writes a corner plot plus an .npz holding the chain, and returns
	(posterior median of the full (u, v) vector, flat chain).
	"""
	M_fit = np.asarray(M_fit, dtype=float)
	D = np.asarray(D, dtype=float)
	n_s, n_c = len(samples), len(classes)
	labels = fit_labels(samples, classes)
	y_obs = np.asarray(y_fit, dtype=float)
	# Default: a flat 5 per cent per bin. The variant that uses the formal
	# uncertainties of the measurement passes them in instead.
	sigma = 0.05 * y_obs if sigma is None else np.asarray(sigma, dtype=float)
	rng = np.random.default_rng(seed)

	def unpack(theta):
		a, b = theta[:n_s - 1], theta[n_s - 1:]
		return np.concatenate([a, [1.0 - np.sum(a)], b, [1.0 - np.sum(b)]])

	def log_prob(theta):
		a, b = theta[:n_s - 1], theta[n_s - 1:]
		if np.any(theta < 0.0) or np.sum(a) > 1.0 or np.sum(b) > 1.0:
			return -np.inf
		m = bilinear_model(unpack(theta), M_fit, D, n_s)
		return -0.5 * np.sum(((y_obs - m) / sigma) ** 2)

	n_dim = (n_s - 1) + (n_c - 1)
	p0 = np.column_stack([rng.dirichlet(np.ones(n_s), size=n_walkers)[:, :n_s - 1],
						 rng.dirichlet(np.ones(n_c), size=n_walkers)[:, :n_c - 1]])
	sampler = emcee.EnsembleSampler(n_walkers, n_dim, log_prob)
	# emcee draws its proposals from its OWN rng, so seeding the walker start is
	# not enough to make a run reproducible -- seed the sampler as well
	sampler.random_state = np.random.RandomState(seed).get_state()
	sampler.run_mcmc(p0, n_steps, progress=False)
	chain = sampler.get_chain(discard=n_burn, flat=True)
	flat = np.array([unpack(t) for t in chain])

	pc = np.percentile(flat, [16, 50, 84], axis=0)
	dg = fit_diagnostics(sampler, n_burn, y_obs, sigma,
						 bilinear_model(pc[1], M_fit, D, n_s), n_dim, labels)
	print('   chi2 = %.1f / %d dof = %.2f      acceptance = %.3f'
		 % (dg['chi2'], dg['ndof'], dg['chi2_red'], dg['acceptance']))
	print('   tau  = %s' % np.array2string(dg['tau'], precision=0))
	print('   Neff = %s   Rhat = %s'
		 % (np.array2string(dg['n_eff'], precision=0),
			np.array2string(dg['rhat'], precision=3)))
	for i, lab in enumerate(labels):
		print('   %-16s = %.4f -%.4f +%.4f'
			 % (lab, pc[1, i], pc[1, i] - pc[0, i], pc[2, i] - pc[1, i]))

	fig = corner.corner(flat, labels=labels, quantiles=[0.16, 0.5, 0.84],
						show_titles=True, title_fmt='.3f')
	fig.savefig(corner_path)
	plt.close(fig)
	print(corner_path, 'written')
	np.savez(corner_path.replace('W_CORNER_', 'W_CHAIN_').replace('.png', '.npz'),
			chain=flat, percentiles=pc, sigma=sigma, autocorr=dg['tau'],
			**{k: v for k, v in dg.items() if k != 'labels'},
			labels=np.array(labels))
	return pc[1], flat, dg


# enclosed probability of the 0.5, 1, 1.5, 2, 2.5 and 3 sigma contours of a
# 2D Gaussian: 1 - exp(-n**2/2)
CONTOUR_LEVELS = tuple(1.0 - np.exp(-0.5 * np.arange(0.5, 3.01, 0.5) ** 2))


def corner_compare(chains, names, colors, labels, corner_path):
	"""Overlay several posteriors, sharing the same parameters, on one corner plot.

	A common range is imposed across the chains so the panels line up, each model
	keeps its own colour, and the legend names them. Medians are printed rather
	than written as panel titles, which would collide between models.
	"""
	stack = np.vstack(chains)
	rng = []
	for i in range(stack.shape[1]):
		lo, hi = float(stack[:, i].min()), float(stack[:, i].max())
		pad = 0.02 * (hi - lo) if hi > lo else 1e-3
		rng.append((lo - pad, hi + pad))
	fig = None
	for ch, col in zip(chains, colors):
		fig = corner.corner(ch, labels=labels, color=col, fig=fig, range=rng,
																				  levels=CONTOUR_LEVELS, quantiles=[0.16, 0.5, 0.84],
																				  show_titles=False, plot_datapoints=False,
																				  plot_density=False, fill_contours=False,
																				  hist_kwargs=dict(color=col, lw=1.6, density=True),
																				  contour_kwargs=dict(linewidths=1.4))
	fig.legend([Line2D([], [], color=c, lw=2.5) for c in colors], names,
											 loc='upper right', frameon=False, fontsize=13)
	fig.savefig(corner_path)
	plt.close(fig)
	print(corner_path, 'written')
	for nm, ch in zip(names, chains):
		pc = np.percentile(ch, [16, 50, 84], axis=0)
		print('   %-24s %s' % (nm, '  '.join(
			'%.4f -%.4f +%.4f' % (pc[1, i], pc[1, i] - pc[0, i], pc[2, i] - pc[1, i])
			for i in range(ch.shape[1]))))


def component_draws(chain, components_fun, n_draw=200, seed=7):
	"""Model components evaluated over n_draw random posterior samples."""
	rng = np.random.default_rng(seed)
	idx = rng.choice(len(chain), size=min(n_draw, len(chain)), replace=False)
	return np.array([np.asarray(components_fun(*chain[i])) for i in idx])


def band_fill(theta_mid, lo, hi, style):
	"""Shade a 68 per cent credible band beneath its component's line."""
	plt.fill_between(np.asarray(theta_mid) * _U['fx'],
																	 np.asarray(lo) * _U['fy'], np.asarray(hi) * _U['fy'],
																	 color=style['color'], alpha=0.25, lw=0, zorder=1)


# Published best-fit models of Comparat et al. 2025 (A&A 697, A173). The
# lsdr10_clustering data root is not exported here, so fall back to its usual
# location.
C25_DIR = os.path.join(os.environ.get('GIT_DR10W_DATA',
							 '/home/comparat/data/legacysurvey/lsdr10_clustering_data'),
							 'data', 'BGS-xcorr-res_XPS')

# Only samples whose stored curve reproduces Table 4 of the paper are listed.
# The 11.0 files on disk carry p0=1.630, p1=0.880 with chi2/ndof=105, and the
# companion parameter file says 1.610 / 2.210, none of which match the
# published 1.634 / 1.221 -- so that sample is deliberately absent and its
# comparison figure is skipped rather than drawn from an unpublished fit.
C25_BASENAME = {
	'10.0': 'LS10_VLIM_ANY_10.0_Mstar_12.0_0.05_z_0.18_N_2759238',
	'10.5': 'LS10_VLIM_ANY_10.5_Mstar_12.0_0.05_z_0.26_N_3263228',
	}

# C25 column -> the component of our model it pairs with
C25_PAIRS = [('1h_cen_xps', 'AGN', 'cen'),
												 ('1h_sat_xps', 'AGN', 'sat'),
												 ('1h_cen_evt', 'GAS', 'cen'),
												 ('1h_sat_evt', 'GAS', 'sat')]

STY_OURS = dict(lw=4.0, alpha=0.35, zorder=2)   # this work: thick, translucent
STY_C25 = dict(lw=1.8, alpha=1.0, zorder=4)     # Comparat+25: thin, crisp


def read_c25_model(m0_str):
	"""Comparat et al. 2025 best-fit model for one sample, or None if absent.

	The stored components are w(theta); the caller scales them to surface
	brightness with the same BG_val used for our own curves.
	"""
	if m0_str not in C25_BASENAME:
		return None
	p = os.path.join(C25_DIR, C25_BASENAME[m0_str] + '-DMm-best-model-halo.fits')
	if not os.path.isfile(p):
		print('   C25 model file missing:', p)
		return None
	return Table.read(p)


def pair_weight(m0_str, tag):
	"""Weight of one (galaxy sample x event class) template in the mixture.

	A template's share of the galaxy-photon pairs scales as the product of the
	two catalogue sizes, N_gal x N_evt (N2_data x N_data). The event catalogue is
	shared between the central and the satellite template of a given class, so it
	is N_gal -- N_cen or N_sat -- that sets their relative importance: satellites
	are 23 to 34 per cent of each stellar mass sample and enter with that much
	less weight.
	"""
	t = PRED['m' + m0_str + tag]
	return np.asarray(t['N_data'], dtype=float) * np.asarray(t['N2_data'], dtype=float)


def make_c25_comparison_fig(m0_str, str_title, N_TOT, x_data, y_data, in_fit, out_path):
	"""Our prediction against the Comparat et al. 2025 published best fit.

	Both are drawn in the active (surface brightness) unit system and share the
	component colours: ours thick and translucent, theirs thin and crisp, so each
	pair can be read off directly. C25 supplies x already in kpc, so only its
	y values are converted.
	"""
	c25 = read_c25_model(m0_str)
	if c25 is None:
		print('   no published C25 model for m' + m0_str + ', comparison skipped')
		return
	ax, ax_r = start_wtheta_fig(residual=False)
	plot_obs(x_data, y_data, in_fit)
	OURS = [('C_GxB', 'BKG', 'cen', L_BKG_C), ('S_GxB', 'BKG', 'sat', L_BKG_S),
			('C_GxA', 'AGN', 'cen', L_AGN_C), ('S_GxA', 'AGN', 'sat', L_AGN_S),
			('C_GxG', 'GAS', 'cen', L_GAS_C), ('S_GxG', 'GAS', 'sat', L_GAS_S)]
	for tag, comp, samp, lab in OURS:
		t = PRED['m' + m0_str + tag]
		pl(t['theta_mid'], t['wtheta'] * pair_weight(m0_str, tag) / N_TOT, label=lab,
		   **sty(comp, samp, **STY_OURS))
	# our total, to pair with the C25 al_gal_evt curve
	our_sum = sum(PRED['m' + m0_str + tg]['wtheta'] * pair_weight(m0_str, tg)
				   for tg, _c, _s, _l in OURS) / N_TOT
	pl(PRED['m' + m0_str + '_GxG']['theta_mid'], our_sum, label=L_SUM,
	   **dict(STY_SUM, **STY_OURS))
	x25 = np.asarray(c25['x_kpc'])
	for col, comp, samp in C25_PAIRS:
		plt.plot(x25, np.asarray(c25[col]) * _U['fy'], **sty(comp, samp, **STY_C25))
	plt.plot(x25, np.asarray(c25['al_gal_evt']) * _U['fy'], **dict(STY_SUM, **STY_C25))
	proxy = [Line2D([], [], color='k', lw=STY_OURS['lw'], alpha=STY_OURS['alpha']),
								 Line2D([], [], color='k', lw=STY_C25['lw'])]
	finish_wtheta_fig(out_path, str_title, ax, ax_r, x_data, y_data, in_fit,
																	  extra_legend=(proxy, ['this work',
																																							  r'Comparat+25 ($\alpha_{SR}$=%.3f, $\Delta_{M_{min}}$=%.3f)'
																																							  % (c25['p0'][0], c25['p1'][0])]))


def make_wtheta_figures(m0, z1, z_mean, BG_val_0):
	"""Write the three w(theta) figures for one stellar mass selected sample.

	m0, z1 : lower stellar mass bound and upper redshift bound of the sample.

	Produces, with m0_str = str(round(m0, 1)) :
	  W_PRED_GALxEVT_m*.png         components predicted from the mock, each
	                                weighted by its N_gal x N_evt share of all pairs
	  W_FIT_GALxEVT_m*.png          3 parameter mixture fit (f_A, f_B, f_G),
	                                centrals and satellites for every component
	  W_FITnoSatAGN_GALxEVT_m*.png  the same fit without the satellite AGN and
	                                satellite background terms
	"""
	m0_str = str(np.round(m0,1))
	z1_str = str(np.round(z1,2))
	str_title = r'$z<$'+z1_str+', '+m0_str+r'$<\log_{10}(M*[M_\odot])$'

	def W(tag):
		return pair_weight(m0_str, tag)

	WTH = XCORR['m'+m0_str]
	#plt.errorbar(WTH['theta'], WTH['wtheta'], yerr= WTH['wtheta']*0.05, label='Obs', color='black', marker='*', markersize=10)
	x_data = WTH['theta']
	y_data = WTH['wtheta']

	# The first 3 and last 4 observed points are excluded from the fit. They are
	# still drawn, with hollow symbols, so the model can be judged outside the
	# fitted range as well.
	FIT_SLICE = slice(3, -4)
	x_fit = WTH['theta'][FIT_SLICE]
	y_fit = WTH['wtheta'][FIT_SLICE]
	e_fit = WTH['wtheta_err'][FIT_SLICE]
	in_fit = np.zeros(len(x_data), dtype=bool)
	in_fit[FIT_SLICE] = True

	# Every figure is drawn twice: as w(theta) against theta, and converted to
	# surface brightness against proper separation. The conversion is a pure
	# rescaling of both axes, so the fits are untouched.
	BG_val = BG_val_0 * (1 + z_mean) ** 2
	cv_r = cosmo.kpc_proper_per_arcmin(z_mean).to(u.kpc / u.deg).value
	UNITS = [UNIT_WTH,
					  dict(tag='SX_', fx=cv_r, fy=BG_val, xlim=(9, 5000), ylim=(3e33, 3e37),
							  xlabel=r'$r_p$ (kpc)', ylabel=r'$S_x$ (erg/kpc$^2$/s)')]

	# N_TOT = W('_GxA')+W('_GxB')+W('_GxG')
	N_TOT = (W('C_GxA')+
			 W('S_GxA')+
			 W('C_GxB')+
			 W('S_GxB')+
			 W('C_GxG')+
			 W('S_GxG'))

	# -------------------------------------------------------------------
	# prediction : every component weighted by its share of the pairs
	# -------------------------------------------------------------------
	for U in UNITS:
		use_units(U)
		ax, ax_r = start_wtheta_fig(residual=False)
		plot_obs(x_data, y_data, in_fit)

		pl(PRED['m'+m0_str+'C_GxB']['theta_mid'], PRED['m'+m0_str+'C_GxB']['wtheta']*W('C_GxB')/N_TOT, label=L_BKG_C, **sty('BKG', 'cen'))
		pl(PRED['m'+m0_str+'S_GxB']['theta_mid'], PRED['m'+m0_str+'S_GxB']['wtheta']*W('S_GxB')/N_TOT, label=L_BKG_S, **sty('BKG', 'sat'))

		pl(PRED['m'+m0_str+'C_GxA']['theta_mid'], PRED['m'+m0_str+'C_GxA']['wtheta']*W('C_GxA')/N_TOT, label=L_AGN_C, **sty('AGN','cen'))
		pl(PRED['m'+m0_str+'S_GxA']['theta_mid'], PRED['m'+m0_str+'S_GxA']['wtheta']*W('S_GxA')/N_TOT, label=L_AGN_S, **sty('AGN','sat'))

		pl(PRED['m'+m0_str+'C_GxG']['theta_mid'], PRED['m'+m0_str+'C_GxG']['wtheta']*W('C_GxG')/N_TOT, label=L_GAS_C, **sty('GAS','cen'))
		pl(PRED['m'+m0_str+'S_GxG']['theta_mid'], PRED['m'+m0_str+'S_GxG']['wtheta']*W('S_GxG')/N_TOT, label=L_GAS_S, **sty('GAS','sat'))

		# pl(PRED['m'+m0_str+'_GxGAB']['theta_mid'], PRED['m'+m0_str+'_GxGAB']['wtheta']*W('_GxGAB')/N_TOT, label=L_GAB, **STY_GAB)
		# pl(PRED['m'+m0_str+'_GxGA']['theta_mid'], PRED['m'+m0_str+'_GxGA']['wtheta']*W('_GxGA')/N_TOT, label=L_GA, **STY_GA)

		all_3 = (PRED['m'+m0_str+'C_GxA']['wtheta']*W('C_GxA')+
				PRED['m'+m0_str+'S_GxA']['wtheta']*W('S_GxA')+
				PRED['m'+m0_str+'C_GxB']['wtheta']*W('C_GxB')+
				PRED['m'+m0_str+'S_GxB']['wtheta']*W('S_GxB')+
				PRED['m'+m0_str+'C_GxG']['wtheta']*W('C_GxG')+
				PRED['m'+m0_str+'S_GxG']['wtheta']*W('S_GxG') ) / (
					N_TOT
					)
		pl(PRED['m'+m0_str+'_GxG']['theta_mid'], (all_3), label=L_SUM, **STY_SUM)
		# pl(PRED['m'+m0_str+'_GxG']['theta_mid'], (all_3)/2, label='all/2', color='blue')
		finish_wtheta_fig(os.path.join(dir_fig, U['tag']+'W_PRED_GALxEVT_m'+m0_str+'.png'), str_title,
																		  ax, ax_r, x_data, y_data, in_fit)

	# the published Comparat et al. 2025 model is tabulated in surface
	# brightness units, so the comparison is only drawn there
	use_units(UNITS[1])
	make_c25_comparison_fig(m0_str, str_title, N_TOT, x_data, y_data, in_fit,
																							  os.path.join(dir_fig, 'SX_W_C25_GALxEVT_m'+m0_str+'.png'))

	# -------------------------------------------------------------------
	# mixture fits
	#
	# A weight is the normalised effective count of its population, so D carries
	# only the mask and M only the templates. M is interpolated onto the fit grid
	# once instead of rebuilding interp1d inside every likelihood call.
	# -------------------------------------------------------------------
	_tm = PRED['m'+m0_str+'_GxG']['theta_mid']
	def _ip(arr):
		return interp1d(_tm, np.asarray(arr))(x_fit)

	# the AGN template is optionally truncated to the PSF core
	null_psf = np.ones_like(PRED['m'+m0_str+'C_GxA']['wtheta'])
	null_psf[PRED['m'+m0_str+'_GxA']['theta_mid']>=0.02] = 0

	def _template(s, c, psf_core):
		"""The template w_gX(theta) as it enters the model, on the native grid.

		The catalogue sizes are NOT folded in here: a weight is itself the
		normalised effective count of its population, so the pair weight of a
		template is the bare product f_g f_X.
		"""
		tag = TEMPLATE_TAG[(s, c)]
		w = np.asarray(PRED['m'+m0_str+tag]['wtheta'], dtype=float)
		if c == 'BKG':
			# uncorrelated with the galaxies: a constant, taken as the median.
			# The outermost bins are NaN (no random pairs), hence nanmedian.
			w = np.full_like(w, np.nanmedian(w))
		if psf_core and c == 'AGN':
			w = w * null_psf
		return w

	def _blocks(samples, classes, drop=(), psf_core=False):
		"""(M, D, M_fit, keys) of one flavour. Dropped templates are zero in both."""
		M = np.array([[np.zeros_like(np.asarray(_tm, dtype=float)) if (s, c) in drop
					  else _template(s, c, psf_core)
					  for c in classes] for s in samples])
		D = np.array([[0.0 if (s, c) in drop else 1.0
					  for c in classes] for s in samples])
		M_fit = np.array([[_ip(m) for m in row] for row in M])
		keys = [(s, c) for s in samples for c in classes if (s, c) not in drop]
		return M, D, M_fit, keys

	def _components_fun(M, D, samples, classes, keys):
		"""Per-template contribution to the model, in the order of keys."""
		si = {s: i for i, s in enumerate(samples)}
		ci = {c: j for j, c in enumerate(classes)}
		def fun(*p):
			u = np.asarray(p[:len(samples)], dtype=float)
			v = np.asarray(p[len(samples):], dtype=float)
			den = float(u @ D @ v)
			return np.array([u[si[s]] * v[ci[c]] * M[si[s], ci[c]] / den
							for s, c in keys])
		return fun

	def native_weights(samples, classes):
		"""The weights the unmodified mock corresponds to: the count fractions.

		A weight is the normalised effective count of its population, so the direct
		prediction sits at f_g = N_g / sum(N_g) and f_X = N_X / sum(N_X), not at a
		common value.
		"""
		ng = np.array([float(PRED['m'+m0_str+TEMPLATE_TAG[(s, 'GAS')]]['N2_data'][0])
					   for s in samples])
		nx = np.array([float(PRED['m'+m0_str+TEMPLATE_TAG[('cen', c)]]['N_data'][0])
					   for c in classes])
		return np.concatenate([ng / ng.sum(), nx / nx.sum()])

	def do_fit(name, tag, samples, classes, drop=(), psf_core=False):
		"""Run one flavour, draw its two w(theta) figures, return its flat chain."""
		M, D, M_fit, keys = _blocks(samples, classes, drop, psf_core)
		print('fitting', m0_str, ':', name)
		print('   native  ' + '  '.join('%s=%.4f' % (l, v) for l, v in zip(
			fit_labels(samples, classes), native_weights(samples, classes))))
		p_med, chain, dg = run_mcmc_fit(M_fit, D, samples, classes, x_fit, y_fit,
								   os.path.join(dir_fig, 'W_CORNER_'+tag+'_GALxEVT_m'+m0_str+'.png'),
								   sigma=SIGMA_FUN(y_fit, e_fit))
		dg['native'] = native_weights(samples, classes)
		dg['median'] = np.asarray(p_med, dtype=float)
		SUMMARY[(m0_str, tag)] = dg
		comps = _components_fun(M, D, samples, classes, keys)
		# the components and the likelihood are separate rewrites of the same
		# model -- verify they agree before either is trusted
		d = np.max(np.abs(_ip(comps(*p_med).sum(axis=0))
						 - bilinear_model(p_med, M_fit, D, len(samples))))
		if not d < 1e-10:
			raise RuntimeError('components disagree with the model by %.3e' % d)
		out = comps(*p_med)
		draws = component_draws(chain, comps)
		lo, hi = np.percentile(draws, [16, 84], axis=0)
		slo, shi = np.percentile(draws.sum(axis=1), [16, 84], axis=0)
		for U in UNITS:
			use_units(U)
			ax, ax_r = start_wtheta_fig()
			plot_obs(x_data, y_data, in_fit)
			for k, (s, c) in enumerate(keys):
				band_fill(_tm, lo[k], hi[k], sty(c, s))
			band_fill(_tm, slo, shi, STY_SUM)
			for k, (s, c) in enumerate(keys):
				pl(_tm, out[k], label=L_TEMPLATE[(s, c)], **sty(c, s))
			pl(_tm, out.sum(axis=0), label=L_SUM, **STY_SUM)
			finish_wtheta_fig(os.path.join(dir_fig, U['tag']+'W_'+tag+'_GALxEVT_m'+m0_str+'.png'),
							  str_title, ax, ax_r, x_data, y_data, in_fit, _tm,
							  out.sum(axis=0), slo, shi)
		return chain

	# five weights: one per galaxy sample, one per event class
	chain_FULL = do_fit('centrals and satellites', 'FIT', SAMPLES, CLASSES)

	# no background, and the AGN kept only around centrals and inside the PSF core
	do_fit('no background', 'FITnoBKG', SAMPLES, ('AGN', 'GAS'),
		   drop=(('sat', 'AGN'),), psf_core=True)

	# no satellite AGN and no satellite background: only the gas keeps a
	# satellite term
	chain_nS = do_fit('no satellite AGN', 'FITnoSatAGN', SAMPLES, CLASSES,
					  drop=(('sat', 'AGN'), ('sat', 'BKG')))

	# the two five-parameter models share a parameter space, so their posteriors
	# can be read against each other directly
	print('comparing', m0_str, ': cen+sat vs no satellite AGN')
	corner_compare([chain_FULL, chain_nS],
				   ['centrals and satellites', 'no satellite AGN'],
				   [C_TOT, C_CMP],
				   fit_labels(SAMPLES, CLASSES),
				   os.path.join(dir_fig, 'W_CORNER_COMPARE_GALxEVT_m'+m0_str+'.png'))


def write_summary(path=None):
	"""Dump every fit of this run as one plain-text table.

	Three blocks per fit: goodness of fit, sampler convergence, and the weights
	with the effective catalogue sizes they imply. The effective counts are the
	shares multiplied by the total of their block, which is what makes a fitted
	weight readable as a number of galaxies or of photons.
	"""
	path = path or os.path.join(dir_fig, 'W_FIT_SUMMARY_' + SIGMA_TAG + '.txt')
	L = []
	L.append('cross-correlation mixture fits')
	L.append('likelihood sigma : ' + SIGMA_LABEL)
	L.append('')
	L.append('%-6s %-12s %8s %5s %9s %7s %8s %8s'
			 % ('sample', 'variant', 'chi2', 'dof', 'chi2/dof', 'accept', 'max tau', 'min Neff'))
	L.append('-' * 74)
	for (m0, tag), d in SUMMARY.items():
		L.append('%-6s %-12s %8.1f %5d %9.2f %7.3f %8.0f %8.0f'
				 % (m0, tag, d['chi2'], d['ndof'], d['chi2_red'],
					d['acceptance'], np.max(d['tau']), np.min(d['n_eff'])))
	L.append('')
	L.append('convergence, per sampled direction')
	L.append('%-6s %-12s %-16s %8s %10s %7s' % ('sample', 'variant', 'weight', 'tau', 'Neff', 'Rhat'))
	L.append('-' * 74)
	for (m0, tag), d in SUMMARY.items():
		# the sampler explores (n_s - 1) + (n_c - 1) directions; the remaining
		# two weights are derived and carry no diagnostic of their own
		for i in range(len(d['tau'])):
			L.append('%-6s %-12s %-16s %8.0f %10.0f %7.3f'
					 % (m0, tag, d['labels'][i], d['tau'][i], d['n_eff'][i], d['rhat'][i]))
	L.append('')
	L.append('weights (posterior median) and effective catalogue sizes')
	L.append('%-6s %-12s %-16s %10s %10s' % ('sample', 'variant', 'weight', 'native', 'fitted'))
	L.append('-' * 74)
	for (m0, tag), d in SUMMARY.items():
		for i, lab in enumerate(d['labels']):
			L.append('%-6s %-12s %-16s %10.4f %10.4f'
					 % (m0, tag, lab, d['native'][i], d['median'][i]))
	txt = '\n'.join(L) + '\n'
	with open(path, 'w') as f:
		f.write(txt)
	print(txt)
	print(path, 'written')


if __name__ == '__main__':
	# One sample can be selected from the command line, e.g. "10.5"; with no
	# argument all three are run.
	_want = sys.argv[1:] if len(sys.argv) > 1 else None
	for _m0, _z1, _zm, _bg in [(10.0, 0.18, 0.136, 6.859*10**36),
							   (10.5, 0.26, 0.191, 6.772*10**36),
							   (11.0, 0.35, 0.252, 6.693*10**36)]:
		if _want is None or str(_m0) in _want:
			make_wtheta_figures(_m0, _z1, _zm, _bg)
	write_summary()
