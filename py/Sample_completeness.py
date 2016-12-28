#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

from astropy.io import fits
import pandas as pd
import matplotlib.pyplot as pl
import seaborn as sns; sns.set_style('ticks')
import numpy as np
from scipy import stats
import matplotlib
from matplotlib.ticker import FormatStrFormatter
params = {
   'axes.labelsize': 10,
   'font.size': 10,
   'legend.fontsize': 10,
   'xtick.labelsize': 10,
   'ytick.labelsize': 10,
   'text.usetex': False,
   'figure.figsize': [9, 4.5]
   }
matplotlib.rcParams.update(params)

def filt_nan(input_array):
    """
    Small functino to remove values in array, not convertible to float. Return the filtered array
    """
    holder = []
    for ii in input_array:
        try:
            holder.append(float(ii))
        except:
            pass
    return np.array(holder)


def main():
    """
    Small script to produce sample completeness plots.
    """


    # Read in burst list
    # burst_table = np.array([[ii, kk, ll, pp, tt] for ii, kk, ll, pp, tt in np.genfromtxt("../data/burst_list.dat", dtype=None)])
    burst_table = pd.read_csv("../data/Burst list - CSV_master.csv")

    name, z, observed, OA, BAT, XRT, HI = burst_table["GRB"].values, burst_table["Redshift"].values, burst_table["Observed"].values, burst_table["Afterglow"].values, burst_table["BAT Fluence (15-150 keV) [10^-7 erg/cm^2]"].values.astype("float"), burst_table["XRT 11 Hour Flux (0.3-10 keV) [10^-11 erg/cm^2/s]"].values.astype("float"), burst_table["XRT Column Density (NH) [10^21 cm^-2]"].values.astype("float")

    BAT_o = BAT[observed == "Yes"]
    XRT_o = XRT[observed == "Yes"]
    HI_o = HI[observed == "Yes"]

    swift_table = pd.read_table("../data/grb_table_1482495106.txt", delimiter="\t", dtype=None)
    name_s, BAT_s, XRT_s, HI_s = swift_table["GRB"].values, swift_table["BAT Fluence (15-150 keV) [10^-7 erg/cm^2]"].values, swift_table["XRT 11 Hour Flux (0.3-10 keV) [10^-11 erg/cm^2/s]"].values, swift_table["XRT Column Density (NH) [10^21 cm^-2]"].values
    print(len(name_s))
    BAT_s = filt_nan(BAT_s)
    XRT_s = filt_nan(XRT_s)
    HI_s = filt_nan(HI_s)


    BAT, BAT_o, BAT_s = BAT[(BAT > 0)], BAT_o[(BAT_o > 0)], BAT_s[(BAT_s > 0)]
    XRT, XRT_o, XRT_s = XRT[(XRT > 0)], XRT_o[(XRT_o > 0)], XRT_s[(XRT_s > 0)]
    HI, HI_o, HI_s = HI[(HI > 0)], HI_o[(HI_o > 0)], HI_s[(HI_s > 0)]

    # Plot
    fig, (ax1, ax2, ax3) = pl.subplots(ncols=3)



    # BAT fluence
    sns.distplot(np.log10(1e-7*BAT), ax=ax1, kde=True)
    sns.distplot(np.log10(1e-7*BAT_s), ax=ax1, kde=True)
    sns.distplot(np.log10(1e-7*BAT_o), ax=ax1, kde=True)
    print(len(BAT[~np.isnan(BAT)]))
    l, m, h = np.percentile(np.log10(1e-7*BAT), [16, 50, 84])
    print(m, m - l, h - m)
    print(len(BAT_o[~np.isnan(BAT_o)]))
    l, m, h = np.percentile(np.log10(1e-7*BAT_o), [16, 50, 84])
    print(m, m - l, h - m)
    print(len(BAT_s[~np.isnan(BAT_s)]))
    l, m, h = np.percentile(np.log10(1e-7*BAT_s), [16, 50, 84])
    print(m, m - l, h - m)

    print(stats.ks_2samp(BAT_s[~np.isnan(BAT_s)], BAT)[1])
    print(stats.ks_2samp(BAT, BAT_o)[1])
    print(stats.ks_2samp(BAT_s[~np.isnan(BAT_s)], BAT_o)[1])




    # XRT flux
    sns.distplot(np.log10(1e-11*XRT), ax=ax2, kde=True)
    sns.distplot(np.log10(1e-11*XRT_s), ax=ax2, kde=True)
    sns.distplot(np.log10(1e-11*XRT_o), ax=ax2, kde=True)



    print(stats.ks_2samp(XRT_s[~np.isnan(XRT_s)], XRT)[1])
    print(stats.ks_2samp(XRT, XRT_o)[1])
    print(stats.ks_2samp(XRT_s[~np.isnan(XRT_s)], XRT_o)[1])

    print(len(XRT[~np.isnan(XRT)]))
    l, m, h = np.percentile(np.log10(1e-11*XRT), [16, 50, 84])
    print(m, m - l, h - m)

    print(len(XRT_s[~np.isnan(XRT_s)]))
    l, m, h = np.percentile(np.log10(1e-11*XRT_s), [16, 50, 84])
    print(m, m - l, h - m)

    print(len(XRT_o[~np.isnan(XRT_o)]))
    l, m, h = np.percentile(np.log10(1e-11*XRT_o), [16, 50, 84])
    print(m, m - l, h - m)



    # HI column
    sns.distplot(np.log10(1e21*HI), ax=ax3, kde=True)
    sns.distplot(np.log10(1e21*HI_s), ax=ax3, kde=True)
    sns.distplot(np.log10(1e21*HI_o), ax=ax3, kde=True)



    print(stats.ks_2samp(HI_s[~np.isnan(HI_s)], HI)[1])
    print(stats.ks_2samp(HI, HI_o)[1])
    print(stats.ks_2samp(HI_s[~np.isnan(HI_s)], HI_o)[1])
    print(len(HI[~np.isnan(HI)]))
    l, m, h = np.percentile(np.log10(1e21*HI), [16, 50, 84])
    print(m, m - l, h - m)
    print(len(HI_s[~np.isnan(HI_s)]))
    l, m, h = np.percentile(np.log10(1e21*HI_s), [16, 50, 84])
    print(m, m - l, h - m)
    print(len(HI_o[~np.isnan(HI_o)]))
    l, m, h = np.percentile(np.log10(1e21*HI_o), [16, 50, 84])
    print(m, m - l, h - m)
    exit()

    for ax in (ax1, ax2, ax3):
        ax.spines['top'].set_visible(False)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        ax.tick_params(axis='x', direction='out')
        ax.tick_params(axis='y', direction='out')

        # # offset the spines
        for spine in ax.spines.values():
            spine.set_position(('outward', 3))

        # put the grid behind
        ax.set_axisbelow(True)



    ax1.set_xlabel(r"log(15-150 keV Fluence) [erg/cm$^2$]")
    ax2.set_xlabel(r"log(0.3-10 keV Flux) [erg/cm$^2$/s]")
    ax3.set_xlabel(r"log(N$_H$) [cm$^2$]")
    ax1.set_ylabel(r"N")

    ax2.set_xlim(-18, -9)
    ax3.set_xlim(18, 24)
    pl.tight_layout()
    ax1.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))


    # pl.show()
    pl.savefig("../document/figures/completeness_BAT.pdf", dpi="figure")
    pl.show()

if __name__ == '__main__':
    main()