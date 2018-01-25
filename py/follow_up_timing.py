#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

from astropy.io import fits
import pandas as pd
import matplotlib; matplotlib.use('TkAgg')

import matplotlib.pyplot as pl
import seaborn as sns; sns.set_style('ticks')
import numpy as np
import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter
params = {
   'axes.labelsize': 15,
   'font.size': 15,
   'legend.fontsize': 15,
   'xtick.labelsize': 15,
   'ytick.labelsize': 15,
   'text.usetex': True,
   'figure.figsize': [6, 6]
   }
mpl.rcParams.update(params)
from astropy.cosmology import Planck15 as cosmo
import astropy.units as u

def main():
    """
    Small script to plot delay time vs acqmag
    """


    # Read in burst list
    burst_table = pd.read_csv("../data/Burst list - observed.csv")
    name, z, delay, mag = burst_table["GRB"].values, burst_table["z"].values, burst_table["Follow-up delay"].values, burst_table["Acquisition mag"].values




    # for kk, ll in list(zip(name, z)):
    #     print(kk, cosmo.luminosity_distance(ll).to(u.cm))
    # exit()

    for ii, ll in enumerate(mag):
        # print(ll)
        if "R" in str(ll):
            print(ll)
            mag[ii] = float(ll.split("=")[1]) + 0.21


    idx_limits = ~(mag == ">24")
    magnondet = 24 * (~idx_limits).astype("int")

    delay_sort = np.argsort(delay)
    sorted_delay = delay[delay_sort]
    sorted_z = z[delay_sort]
    print(len(sorted_z))
    print(sorted_z)
    # exit()
    fractional_completeness = (1 - np.cumsum(np.isnan(sorted_z).astype("int"))/len(sorted_z))*100
    # print(len(sorted_z), np.sum(np.isnan(sorted_z).astype("int")))
    # print(fractional_completeness[-1])
    # exit()
    idx_hasz = ~np.isnan(z)

    # Plot
    fig, ax1 = pl.subplots()


    color = sns.color_palette()[0]
    color_rgb = mpl.colors.colorConverter.to_rgb(color)
    colors = [sns.set_hls_values(color_rgb, l=l) for l in np.linspace(1, 0, 12)]
    cmap = sns.blend_palette(colors, as_cmap=True)
    # cmap = pl.get_cmap("plasma")
    # With redshift and acq mag
    sc = ax1.scatter(delay[idx_hasz & idx_limits], mag[idx_hasz & idx_limits], c=z[idx_hasz & idx_limits], cmap=cmap, s=35)

    # Without redshift and with acq mag
    ax1.scatter(delay[~idx_hasz & idx_limits], mag[~idx_hasz & idx_limits], color=sns.color_palette()[2], s=45, marker="x")
    # With redshift, but without acq mag
    ax1.scatter(delay[idx_hasz & ~idx_limits], magnondet[idx_hasz & ~idx_limits], c=z[idx_hasz & ~idx_limits], cmap=cmap, s=175, marker=r'$\downarrow$')
    # No redshift and without acq mag
    ax1.scatter(delay[~idx_hasz & ~idx_limits], magnondet[~idx_hasz & ~idx_limits], color=sns.color_palette()[2], s=175, marker=r'$\downarrow$')
    ax2 = pl.twinx()
    ax2.plot(sorted_delay, fractional_completeness, color=sns.color_palette()[2])



    ax1.spines['top'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()
    ax1.tick_params(axis='x', direction='out')
    ax1.tick_params(axis='y', direction='out')
    ax2.tick_params(axis='y', direction='out')
    # offset the spines
    for spine in ax1.spines.values():
        spine.set_position(('outward', 3))
    for spine in ax2.spines.values():
        spine.set_position(('outward', 3))
    # put the grid behind
    ax1.set_axisbelow(True)
    ax2.set_axisbelow(True)

    # ax.set_xlim((19.5, 23))
    ax1.set_xlim((0.05, 110))
    ax2.set_ylim((70, 105))
    ax2.set_yticks([75, 80, 85, 90, 95, 100])
    ax1.set_xlabel(r"Follow-up delay (hours)")
    ax1.set_ylabel(r"Acquisition camera magnitude ($R, i$-band)")
    ax2.set_ylabel(r"Redshift completeness (\%)")

    ax1.invert_yaxis()
    ax1.semilogx()
    # pl.tight_layout()
    ax1.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    #create a colorbar axis

    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('top', size='5%', pad=0.05)
    cbar = fig.colorbar(sc, cax=cax, orientation='horizontal')

    # cbar.ax.tick_params(direction='out')
    cax.xaxis.set_ticks_position("top")
    cax.xaxis.set_label_position("top")
    # cbar.make_axes(location="top")
    cbar.set_label("Redshift")


    pl.tight_layout()
    pl.savefig("../document/figures/timing.pdf", dpi="figure")
    pl.show()

if __name__ == '__main__':
    main()