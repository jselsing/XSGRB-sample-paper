#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

from astropy.io import fits
import pandas as pd
import matplotlib; matplotlib.use('TkAgg')
import matplotlib.pyplot as pl
import seaborn as sns; sns.set_style('ticks')
import numpy as np
from scipy import stats, signal, interpolate
import matplotlib as mpl

params = {
   'axes.labelsize': 14,
   'font.size': 14,
   'legend.fontsize': 14,
   'xtick.labelsize': 14,
   'ytick.labelsize': 14,
   'text.usetex': True,
   'figure.figsize': [8.27, 11.7]
   }
mpl.rcParams.update(params)

def main():
    """
    # Script to produce plot of GRB121024A
    """
    # Load data from file
    data = np.genfromtxt("../data/GRB121024A_OB1_stitched_spectrum.dat")
    wl = data[:, 0]
    flux = data[:, 1]
    error = data[:, 2]
    n_elem = len(wl)

    # Get linelist
    import lineid_plot
    fit_line_positions = np.genfromtxt("../data/grblines_Christensen2011.dat", dtype=None)

    linelist = [0]*len(fit_line_positions)
    linenames = [0]*len(fit_line_positions)
    for ii, nn in enumerate(fit_line_positions):
        linelist[ii] = float(nn[0])
        linenames[ii] = nn[1].decode("utf-8")
    linelist = np.array(linelist)
    linenames = np.array(linenames)
    # print(linenames)
    # Redshifts of systems
    z = [2.300]
    ls = ["dashed"]
    # Make plot
    n_axes = 7
    fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7) = pl.subplots(n_axes, 1)

    for ii, ax in enumerate([ax1, ax2, ax3, ax4, ax5, ax6, ax7]):

        cut = slice(int(np.round(n_elem*(ii/n_axes))), int(np.round(n_elem*(ii/n_axes + 1/n_axes))))
        ax.plot(wl[cut], flux[cut], color="black", lw=0.5, linestyle="steps-mid", rasterized=True, zorder=2)
        interp = interpolate.interp1d(wl[cut], signal.medfilt(flux[cut], 11))
        for oo, pp in enumerate(z):
            idx = [ii for ii, kk in enumerate(linelist*(1 + pp)) if (kk > wl[cut][0] and kk < wl[cut][-1])]
            # for aa, ss in zip(linelist[idx]*(1 + pp), linenames[idx]):
                # ax.axvline(aa, ymin=interp(aa)/2, linestyle=ls[oo], color="black", alpha=0.7, lw=0.5)
                # ax.text(aa, 1.5, ss)
            lineid_plot.plot_line_ids(wl[cut], flux[cut], linelist[idx]*(1 + pp), linenames[idx], arrow_tip=1.2, box_loc=1.6, ax = ax)


        ax.plot(wl[cut][::10], signal.medfilt(error[cut][::10], 11), color=sns.color_palette("muted")[0], alpha=0.9, zorder=1)
        ax.axhline(1, linestyle="dashed", color=sns.color_palette("muted")[2], alpha=0.9, zorder=10)
        ax.axhline(0, linestyle="dotted", color="black", alpha=0.7)
        ax.set_xlim(min(wl[cut]), max(wl[cut]))
        ax.set_ylim(-0.2, 2)
        ax.tick_params(axis='x', direction='in')
        # ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        # ax.axvspan(12600, 12800, color = "grey", alpha = 0.4)
        ax.axvspan(13500, 14500, color = "grey", alpha = 0.4)
        ax.axvspan(18000, 19500, color = "grey", alpha = 0.4)
    # Save figure for tex
    # pl.legend()
    ax4.set_ylabel(r"$F_{\lambda}\,\rm{(10^{-17}\,erg\,s^{-1}\,cm^{-2}\, \AA^{-1})}$")
    # ax7.set_xlabel(r"$\rm{Observed\,wavelength\,}$ (\AA)")
    ax7.set_xlabel(r"Observed Wavelength (\AA)")
    pl.tight_layout()
    pl.subplots_adjust(hspace=0.2)
    pl.savefig("../document/figures/GRB121024A.pdf", dpi="figure", rasterize=True)
    pl.show()

if __name__ == '__main__':
    main()