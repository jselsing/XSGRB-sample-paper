#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

from astropy.io import fits
import pandas as pd
import matplotlib.pyplot as pl
import seaborn as sns; sns.set_style('ticks')
import numpy as np
from scipy import stats
import matplotlib as mpl

params = {
   'axes.labelsize': 10,
   'font.size': 10,
   'legend.fontsize': 10,
   'xtick.labelsize': 10,
   'ytick.labelsize': 10,
   'text.usetex': False,
   'figure.figsize': [7.281, 4.5]
   }


def iqr(a):
    """Calculate the IQR for an array of numbers."""
    a = np.asarray(a)
    q1 = stats.scoreatpercentile(a, 25)
    q3 = stats.scoreatpercentile(a, 75)
    return q3 - q1


def _freedman_diaconis_bins(a):
    """Calculate number of hist bins using Freedman-Diaconis rule."""
    # From http://stats.stackexchange.com/questions/798/
    a = np.asarray(a)
    h = 2 * iqr(a) / (len(a) ** (1 / 3))
    # fall back to sqrt(a) bins if iqr is 0
    if h == 0:
        return int(np.sqrt(a.size))
    else:
        return int(np.ceil((a.max() - a.min()) / h))


def filt_nan(input_array):
    """
    Small functino to remove values in array, not convertible to float. Return the filtered array
    """
    holder = []
    for ii in input_array:
        try:
            holder.append(float(ii))
        except:
            holder.append(1e99)
    return np.array(holder)


def main():
    """
    # Script to produce NH distribution plots for dark bursts.
    """
    # Read in burst list
    burst_table = pd.read_csv("../data/Burst list - CSV_master.csv")
    name, z, observed, OA, BAT, XRT, HI = burst_table["GRB"].values, burst_table["Redshift"].values, burst_table["Observed"].values, burst_table["Afterglow"].values, burst_table["BAT Fluence (15-150 keV) [10^-7 erg/cm^2]"].values.astype("float"), burst_table["XRT 11 Hour Flux (0.3-10 keV) [10^-11 erg/cm^2/s]"].values.astype("float"), burst_table["XRT Column Density (NH) [10^21 cm^-2]"].values.astype("float")
    # Bursts that are observed
    name_o = name[observed == "Yes"]
    BAT_o = BAT[observed == "Yes"]
    XRT_o = XRT[observed == "Yes"]
    HI_o = HI[observed == "Yes"]
    # Bursts that are not
    name_c = name[observed == "No"]
    BAT_c = BAT[observed == "No"]
    XRT_c = XRT[observed == "No"]
    HI_c = HI[observed == "No"]
    # Swift catalog
    swift_table = pd.read_table("../data/grb_table_1482495106.txt", delimiter="\t", dtype=None)
    name_s, BAT_s, XRT_s, HI_s = swift_table["GRB"].values, swift_table["BAT Fluence (15-150 keV) [10^-7 erg/cm^2]"].values, swift_table["XRT 11 Hour Flux (0.3-10 keV) [10^-11 erg/cm^2/s]"].values, swift_table["XRT Column Density (NH) [10^21 cm^-2]"].values

    # Greiner Table
    Greiner_table = pd.read_csv("../data/Burst list - Greiner Table.csv")
    name_g, inst, OT = Greiner_table["GRB"].values, Greiner_table["Instrument"].values, Greiner_table["OT"].values

    # Bursts with optical afterglows
    name_g_ot = name_g[(OT == "y") & (inst == "Swift")]
    idx_ot = [ii for ii, kk in enumerate(name_s) if kk in name_g_ot]
    BAT_ot = np.log10(1e-7*filt_nan(BAT_s[idx_ot]))
    XRT_ot = np.log10(1e-11*filt_nan(XRT_s[idx_ot]))
    HI_ot = np.log10(1e21*filt_nan(HI_s[idx_ot]))

    # Bursts without optical afterglows
    name_g_not = name_g[~(OT == "y") & (inst == "Swift")]
    idx_not = [ii for ii, kk in enumerate(name_s) if kk in name_g_not]
    BAT_not = np.log10(1e-7*filt_nan(BAT_s[idx_not]))
    XRT_not = np.log10(1e-11*filt_nan(XRT_s[idx_not]))
    HI_not = np.log10(1e21*filt_nan(HI_s[idx_not]))

    idx = [ii for ii, kk in enumerate(name_g) if kk in name]
    name_g_ot = name_g[idx][(OT[idx] == "y") & (inst[idx] == "Swift")]
    name_g_not = name_g[idx][~(OT[idx] == "y") & (inst[idx] == "Swift")]

    print(idx)

    print((len(name_g_not))/(len(name_g_ot) + len(name_g_not)))
    exit()
    print(len(BAT_ot))
    mask_ot = (BAT_ot > -8) & (BAT_ot < -3) & (HI_ot > 19) & (HI_ot < 24)
    g = sns.JointGrid(x=BAT_ot[mask_ot], y=HI_ot[mask_ot], xlim = (-8, -3), ylim = (19, 24), space=0)
    g = g.plot_marginals(sns.distplot, hist=True, kde=False, norm_hist=True)
    # Import directly from Seaborn source to control
    color = sns.color_palette()[0]
    color_rgb = mpl.colors.colorConverter.to_rgb(color)
    colors = [sns.set_hls_values(color_rgb, l=l) for l in np.linspace(1, 0, 12)]
    cmap = sns.blend_palette(colors, as_cmap=True)
    x_bins = _freedman_diaconis_bins(g.x)
    y_bins = _freedman_diaconis_bins(g.y)
    gridsize = int(np.mean([x_bins, y_bins]))
    g = g.plot_joint(pl.hexbin, gridsize=gridsize, cmap=cmap)


    mask_not = (BAT_not > -8) & (BAT_not < -3) & (HI_not > 19) & (HI_not < 24)
    g.x = BAT_not[mask_not]
    g.y = HI_not[mask_not]
    g = g.plot_marginals(sns.distplot, hist=False)
    g = g.plot_joint(sns.kdeplot)

    print(stats.ks_2samp(HI_not[mask_not], HI_ot[mask_ot]))
    print(len(HI_not[mask_not]))
    l, m, h = np.percentile(HI_not[mask_not], [16, 50, 84])
    print(m, m - l, h - m)
    print(len(HI_ot[mask_ot]))
    l, m, h = np.percentile(HI_ot[mask_ot], [16, 50, 84])
    print(m, m - l, h - m)

    print("XRT KS")
    print(stats.ks_2samp(XRT_not[mask_not], XRT_ot[mask_ot]))
    print("BAT KS")
    print(stats.ks_2samp(BAT_not[mask_not], BAT_ot[mask_ot]))

    g.set_axis_labels(r"log(15-150 keV Fluence) [erg/cm$^2$]", r"log(N$_H$) [cm$^2$]")
    pl.tight_layout()


    # Save figure for tex
    pl.legend()
    pl.savefig("../document/figures/NH_dark.pdf", dpi="figure")
    pl.show()

if __name__ == '__main__':
    main()