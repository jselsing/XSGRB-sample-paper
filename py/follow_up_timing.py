#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

from astropy.io import fits
import matplotlib.pyplot as pl
import seaborn; seaborn.set_style('ticks')
import numpy as np
import matplotlib
from matplotlib.ticker import FormatStrFormatter
params = {
   'axes.labelsize': 10,
   'font.size': 10,
   'legend.fontsize': 10,
   'xtick.labelsize': 10,
   'ytick.labelsize': 10,
   'text.usetex': True,
   'figure.figsize': [4.5, 4.5]
   }
matplotlib.rcParams.update(params)



def main():
    """
    Small script to plot delay time vs acqmag
    """


    # Read in burst list
    burst_table = np.array([[ii, kk, ll, pp, tt] for ii, kk, ll, pp, tt in np.genfromtxt("../data/burst_list.dat", dtype=None)])

    name, z, delay, mag = burst_table[:, 0], burst_table[:, 1].astype(float), burst_table[:, 2].astype(float), burst_table[:, 3].astype(float)

    # Plot
    fig, ax = pl.subplots()

    cmap = pl.get_cmap("plasma")

    sc = ax.scatter(delay, mag, c=z, cmap=cmap, s=25)

    # a, b, c = ax.errorbar(nH, OA_nH, xerr=nH_std, yerr=OA_nHe, fmt="none", marker=None, mew=0, alpha=0.5, lw=0.5)
    # e_color = clb.to_rgba(OA_z)
    # c[0].set_color(e_color)
    # c[1].set_color(e_color)

    # from matplotlib.patches import Ellipse
    # for jj in np.arange(1, 2):
    #     for ii, (kk, ll) in enumerate(zip(nH, OA_nH)):
    #         ax.add_artist(Ellipse((kk, ll), 2*jj*nH_std[ii], 2*jj*OA_nHe[ii], fill=False, linestyle='dashed', lw = 0.5, alpha = 1.0/(2.0*jj) , color=cmap(OA_z[ii]/max(OA_z))))


    # for ii, txt in enumerate(names):
    #     ax.annotate(txt, (nH[ii], OA_nH[ii]), size="x-small")


    # ax.annotate(r"Z/Z$_\odot$ = 1", xy=(20, 20.5), rotation=50)
    ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)
    # ax.spines['left'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.tick_params(axis='x', direction='out')
    ax.tick_params(axis='y', direction='out')
    # offset the spines
    for spine in ax.spines.values():
        spine.set_position(('outward', 3))
    # put the grid behind
    ax.set_axisbelow(True)

    # ax.set_xlim((19.5, 23))
    ax.set_xlim((0.05, 110))
    ax.set_xlabel(r"Follow-up delay [hours]")
    ax.set_ylabel(r"Acquisition camera magnitude [r-band]")

    pl.gca().invert_yaxis()
    ax.semilogx()
    pl.tight_layout()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    #create a colorbar axis
    cax, kw = matplotlib.colorbar.make_axes(ax, location='top')
    clb = pl.colorbar(sc, orientation="horizontal", cax=cax, ticks=[0, 1, 2, 3, 4, 5, 6])
    clb.ax.set_title("Redshift")
    # pl.show()
    pl.savefig("../document/figures/timing.pdf", dpi="figure")


if __name__ == '__main__':
    main()