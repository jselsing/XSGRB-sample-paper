#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

from astropy.io import fits
import pandas as pd
import matplotlib.pyplot as pl
import seaborn as sns; sns.set_style('ticks')
import numpy as np
import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter
params = {
   'axes.labelsize': 12,
   'font.size': 12,
   'legend.fontsize': 12,
   'xtick.labelsize': 12,
   'ytick.labelsize': 12,
   'text.usetex': False,
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
    # burst_table = np.array([[ii, kk, ll, pp, tt] for ii, kk, ll, pp, tt in np.genfromtxt("../data/burst_list.dat", dtype=None)])
    burst_table = pd.read_csv("../data/Burst list - CSV_observed.csv")

    name, z, delay, mag = burst_table["GRB"].values, burst_table["z"].values, burst_table["Follow-up delay"].values, burst_table["Acquisition mag"].values
    



    for kk, ll in list(zip(name, z)):
        print(kk, cosmo.luminosity_distance(ll).to(u.cm))

    exit()

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
    sc = ax1.scatter(delay[idx_hasz & idx_limits], mag[idx_hasz & idx_limits], c=z[idx_hasz & idx_limits], cmap=cmap, s=25)
    # Without redshift and with acq mag
    ax1.scatter(delay[~idx_hasz & idx_limits], mag[~idx_hasz & idx_limits], color="black", s=25)
    # With redshift, but without acq mag
    ax1.scatter(delay[idx_hasz & ~idx_limits], magnondet[idx_hasz & ~idx_limits], c=z[idx_hasz & ~idx_limits], cmap=cmap, s=125, marker=u'$\u2193$')
    # No redshift and without acq mag
    ax1.scatter(delay[~idx_hasz & ~idx_limits], magnondet[~idx_hasz & ~idx_limits], color="black", s=125, marker=u'$\u2193$')

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
    ax1.set_xlabel(r"Follow-up delay [hours]")
    ax1.set_ylabel(r"Acquisition camera magnitude [r-band]")
    ax2.set_ylabel(r"Redshift completeness [%]")

    ax1.invert_yaxis()
    ax1.semilogx()
    pl.tight_layout()
    ax1.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    #create a colorbar axis
    cax, kw = mpl.colorbar.make_axes(ax2, location='top')
    clb = pl.colorbar(sc, orientation="horizontal", cax=cax, ticks=[0, 1, 2, 3, 4, 5, 6])
    cax, kw = mpl.colorbar.make_axes(ax1, location='top')
    clb = pl.colorbar(sc, orientation="horizontal", cax=cax, ticks=[0, 1, 2, 3, 4, 5, 6])
    clb.ax.set_title("Redshift")
    # pl.show()
    pl.savefig("../document/figures/timing.pdf", dpi="figure")
    pl.show()

if __name__ == '__main__':
    main()