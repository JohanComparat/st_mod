
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as plt
import os, sys
import glob
import numpy as np

fig_dir = os.path.join(os.environ['GIT_STMOD'], 'data', 'benchmark', 'Zu_Mandelbaum_2016')
BM_dir = os.path.join(os.environ['GIT_STMOD'], 'data', 'benchmark', 'Zu_Mandelbaum_2016')


figure_name = os.path.join(fig_dir, 'ESDfig_M_106_110_BCRS.png')
plt.figure(0, (4.5,4))
x, y0, y2, y1 = np.loadtxt( os.path.join(BM_dir, 'esd_M106_110_shift_4_BC.txt'), unpack = True)
plt.plot(x, y0/1e4, color='darkblue', ls='dashed')
plt.fill_between(x, y1/1e4, y2/1e4, color='darkblue', alpha=0.5)
x, y2, y1 = np.loadtxt( os.path.join(BM_dir, 'esd_M106_110_shift_5_RS.txt'), unpack = True)
plt.plot(x, (y1/1e5+ y2/1e5)/2., color='darkred', ls='dashed')
plt.fill_between(x, y1/1e5, y2/1e5, color='darkred', alpha=0.5)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$r_p$ [$h^{-1}$Mpc]")
plt.ylabel(r"$\Delta\Sigma(r_p)$ [M$_\odot h $pc$^2$]")
plt.title('Milky Way analog galaxies')
plt.tight_layout()

plt.savefig( figure_name )
plt.clf()
print(figure_name, 'written')



figure_name = os.path.join(fig_dir, 'ESDfig_M_106_110_BCRS_nt.png')
plt.figure(0, (4.5,4))
x, y0, y2, y1 = np.loadtxt( os.path.join(BM_dir, 'esd_M106_110_shift_4_BC.txt'), unpack = True)
plt.plot(x, y0/1e4, color='darkblue', ls='dashed')
plt.fill_between(x, y1/1e4, y2/1e4, color='darkblue', alpha=0.5, label='Star-forming')
x, y2, y1 = np.loadtxt( os.path.join(BM_dir, 'esd_M106_110_shift_5_RS.txt'), unpack = True)
plt.plot(x, (y1/1e5+ y2/1e5)/2., color='darkred', ls='dashed')
plt.fill_between(x, y1/1e5, y2/1e5, color='darkred', alpha=0.5, label='Quiescent')
plt.ylim((0.1, 300))
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$r_p$ [$h^{-1}$Mpc]")
plt.ylabel(r"$\Delta\Sigma(r_p)$ [M$_\odot h $pc$^2$]")
plt.legend()
plt.title('Galaxy-galaxy lensing')
plt.tight_layout()

plt.savefig( figure_name )
plt.clf()
print(figure_name)
