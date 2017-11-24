#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

from astropy.io import fits
import pandas as pd
import matplotlib; matplotlib.use('TkAgg')
import matplotlib.pyplot as pl
import seaborn as sns; sns.set_style('ticks')
import numpy as np
from scipy import stats
import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter

params = {
   'axes.labelsize': 10,
   'font.size': 10,
   'legend.fontsize': 10,
   'xtick.labelsize': 10,
   'ytick.labelsize': 10,
   'text.usetex': True,
   'figure.figsize': [7.281, 4.5]
   }
mpl.rcParams.update(params)

def filt_nan(input_array):
    """
    Small functino to remove values in array, not convertible to float. Return the filtered array
    """
    holder = []
    for ii in input_array:
        try:
            holder.append(float(ii))
        except:
            holder.append(np.nan)

    return np.array(holder)


def main():
    """
    Small script to produce sample completeness plots.
    """


    # Read in burst list
    # burst_table = np.array([[ii, kk, ll, pp, tt] for ii, kk, ll, pp, tt in np.genfromtxt("../data/burst_list.dat", dtype=None)])
    burst_table = pd.read_csv("../data/Burst list - master.csv")

    name, z, observed, OA, BAT, XRT, HI = burst_table["GRB"].values, burst_table["Redshift"].values, burst_table["Observed"].values, burst_table["Afterglow"].values, burst_table["BAT Fluence (15-150 keV) [10^-7 erg/cm^2]"].values.astype("float"), burst_table["XRT 11 Hour Flux (0.3-10 keV) [10^-11 erg/cm^2/s]"].values.astype("float"), burst_table["XRT Column Density (NH) [10^21 cm^-2]"].values.astype("float")

    BAT_o = BAT[observed == "Yes"]
    XRT_o = XRT[observed == "Yes"]
    HI_o = HI[observed == "Yes"]

    BAT_c = BAT#[observed == "No"]
    XRT_c = XRT#[observed == "No"]
    HI_c = HI#[observed == "No"]

    swift_table = pd.read_table("../data/grb_table_1511519199.txt", delimiter="\t", dtype=None)
    name_s, BAT_s, XRT_s, HI_s = swift_table["GRB"].values, swift_table["BAT Fluence (15-150 keV) [10^-7 erg/cm^2]"].values, swift_table["XRT 11 Hour Flux (0.3-10 keV) [10^-11 erg/cm^2/s]"].values, swift_table["XRT Column Density (NH) [10^21 cm^-2]"].values
    print(len(name_s))
    exit()
    # Exclude sample bursts
    idx = [ii for ii, kk in enumerate(name_s) if kk not in name]

    BAT_s = filt_nan(BAT_s[idx])
    XRT_s = filt_nan(XRT_s[idx])
    HI_s = filt_nan(HI_s[idx])


    BAT_c, BAT_o, BAT_s = BAT_c[(BAT_c > 0)], BAT_o[(BAT_o > 0)], BAT_s[(BAT_s > 0)]
    XRT_c, XRT_o, XRT_s = XRT_c[(XRT_c > 0)], XRT_o[(XRT_o > 0)], XRT_s[(XRT_s > 0)]
    HI_c, HI_o, HI_s = HI_c[(HI_c > 0)], HI_o[(HI_o > 0)], HI_s[(HI_s > 0)]

    # Plot
    fig, (ax1, ax2, ax3) = pl.subplots(ncols=3)



    # BAT fluence
    sns.distplot(np.log10(1e-7*BAT_s), ax=ax1, kde=False, norm_hist=True, hist_kws={"histtype": "step", "alpha": 1, "linewidth": 3, "linestyle": "dashed"}, label="Full Swift sample")
    sns.distplot(np.log10(1e-7*BAT_c), ax=ax1, kde=False, norm_hist=True, hist_kws={"histtype": "step", "alpha": 1, "linewidth": 3, "linestyle": "dashed"}, label="Complete sample")
    sns.distplot(np.log10(1e-7*BAT_o), ax=ax1, kde=False, norm_hist=True, hist_kws={"histtype": "step", "alpha": 1, "linewidth": 3, "linestyle": "dashed"}, label="Observed sample")
    print(len(BAT_c[~np.isnan(BAT_c)]))
    l, m, h = np.percentile(np.log10(1e-7*BAT_c), [16, 50, 84])
    print(m, m - l, h - m)
    print(len(BAT_o[~np.isnan(BAT_o)]))
    l, m, h = np.percentile(np.log10(1e-7*BAT_o), [16, 50, 84])
    print(m, m - l, h - m)
    print(len(BAT_s[~np.isnan(BAT_s)]))
    l, m, h = np.percentile(np.log10(1e-7*BAT_s), [16, 50, 84])
    print(m, m - l, h - m)

    print(stats.ks_2samp(BAT_s[~np.isnan(BAT_s)], BAT_c)[1])
    print(stats.ks_2samp(BAT_c, BAT_o)[1])
    print(stats.ks_2samp(BAT_s[~np.isnan(BAT_s)], BAT_o)[1])
    print()



    # XRT flux
    sns.distplot(np.log10(1e-11*XRT_s), ax=ax2, kde=False, norm_hist=True, hist_kws={"histtype": "step", "alpha": 1, "linewidth": 3, "linestyle": "dashed"}, label="Full Swift sample")
    sns.distplot(np.log10(1e-11*XRT_c), ax=ax2, kde=False, norm_hist=True, hist_kws={"histtype": "step", "alpha": 1, "linewidth": 3, "linestyle": "dashed"}, label="Complete sample")
    sns.distplot(np.log10(1e-11*XRT_o), ax=ax2, kde=False, norm_hist=True, hist_kws={"histtype": "step", "alpha": 1, "linewidth": 3, "linestyle": "dashed"}, label="Observed sample")



    print(stats.ks_2samp(XRT_s[~np.isnan(XRT_s)], XRT_c)[1])
    print(stats.ks_2samp(XRT_c, XRT_o)[1])
    print(stats.ks_2samp(XRT_s[~np.isnan(XRT_s)], XRT_o)[1])

    print(len(XRT_c[~np.isnan(XRT_c)]))
    l, m, h = np.percentile(np.log10(1e-11*XRT_c), [16, 50, 84])
    print(m, m - l, h - m)

    print(len(XRT_s[~np.isnan(XRT_s)]))
    l, m, h = np.percentile(np.log10(1e-11*XRT_s), [16, 50, 84])
    print(m, m - l, h - m)

    print(len(XRT_o[~np.isnan(XRT_o)]))
    l, m, h = np.percentile(np.log10(1e-11*XRT_o), [16, 50, 84])
    print(m, m - l, h - m)

    burst_table = pd.read_csv("../data/Burst list - master.csv")

    name, z, observed, OA, NHX, NHXh, NHXl = burst_table["GRB"].values, burst_table["Redshift"].values, burst_table["Observed"].values, burst_table["Afterglow"].values, burst_table["Ave_NH_PC"].values.astype("float"), burst_table["Ave_dNH_PC+"].values.astype("float"), burst_table["Ave_dNH_PC-"].values.astype("float")

    print(len(z), len(observed))
    z = filt_nan(z)
    print(len(z))

    NHX_o = NHX[(observed == "Yes") & (~np.isnan(z))]
    NHXh_o = NHXh[(observed == "Yes") & (~np.isnan(z))]
    NHXl_o = NHXl[(observed == "Yes") & (~np.isnan(z))]

    NHX_c = NHX[~np.isnan(z)]#[observed == "No"]
    NHXh_c = NHXh[~np.isnan(z)]#[observed == "No"]
    NHXl_c = NHXl[~np.isnan(z)]#[observed == "No"]

    swift_table = pd.read_table("../data/allSwiftGRBs_NH.txt", delimiter="\t", dtype=None)
    # print(swift_table)
    name_s, z_s, NHX_s, NHXh_s, NHXl_s = swift_table["#GRB"].values, swift_table["Redshift"].values, swift_table["Ave_NH_PC"].values, swift_table["Ave_dNH_PC+"].values, swift_table["Ave_dNH_PC-"].values
    # Exclude sample bursts
    idx = [ii for ii, kk in enumerate(name_s) if kk not in name and ~np.isnan(z_s[ii])]

    NHX_s = filt_nan(NHX_s[idx])
    NHXh_s = filt_nan(NHXh_s[idx])
    NHXl_s = filt_nan(NHXl_s[idx])


    NHX_c, NHX_o, NHX_s = NHX_c[(NHX_c > 0)], NHX_o[(NHX_o > 0)], NHX_s[(NHX_s > 0)]
    NHXh_c, NHXh_o, NHXh_s = NHXh_c[(NHXh_c > 0)], NHXh_o[(NHXh_o > 0)], NHXh_s[(NHXh_s > 0)]
    NHXl_c, NHXl_o, NHXl_s = NHXl_c[(NHXl_c > 0)], NHXl_o[(NHXl_o > 0)], NHXl_s[(NHXl_s > 0)]


    # Rescale
    NHX_s = np.log10(1e22*NHX_s)
    # NHX_s = NHX_s[NHX_s > 10]
    NHX_c = np.log10(1e22*NHX_c)
    # NHX_c = NHX_c[NHX_c > 10]
    NHX_o = np.log10(1e22*NHX_o)
    # NHX_o = NHX_o[NHX_o > 10]

    # HI column
    sns.distplot(NHX_s, ax=ax3, kde=False, norm_hist=True, hist_kws={"histtype": "step", "alpha": 1, "linewidth": 3, "linestyle": "dashed"}, label="Full Swift sample")
    sns.distplot(NHX_c, ax=ax3, kde=False, norm_hist=True, hist_kws={"histtype": "step", "alpha": 1, "linewidth": 3, "linestyle": "dashed"}, label="Complete sample")
    sns.distplot(NHX_o, ax=ax3, kde=False, norm_hist=True, hist_kws={"histtype": "step", "alpha": 1, "linewidth": 3, "linestyle": "dashed"}, label="Observed sample")



    print(stats.ks_2samp(NHX_s[~np.isnan(NHX_s)], NHX_c)[1])
    print(stats.ks_2samp(NHX_c, NHX_o)[1])
    print(stats.ks_2samp(NHX_s[~np.isnan(NHX_s)], NHX_o)[1])


    print(len(NHX_s[~np.isnan(NHX_s)]))
    l, m, h = np.percentile(NHX_s, [16, 50, 84])
    print(m, m - l, h - m)
    print(len(NHX_c[~np.isnan(NHX_c)]))
    l, m, h = np.percentile(NHX_c, [16, 50, 84])
    print(m, m - l, h - m)
    print(len(NHX_o[~np.isnan(NHX_o)]))
    l, m, h = np.percentile(NHX_o, [16, 50, 84])
    print(m, m - l, h - m)
    # exit()

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
        # ax.set_ylim(0, 1)



    ax1.set_xlabel(r"log(15-150 keV Fluence) [erg/cm$^2$]")
    ax2.set_xlabel(r"log(0.3-10 keV Flux) [erg/cm$^2$/s]")
    ax3.set_xlabel(r"log(N$_H$) [cm$^2$]")
    # ax1.set_ylabel(r"N")

    ax2.set_xlim(-18, -9)
    # ax1.set_xlim(-10, -3)
    # ax2.set_xlim(-20, -9)
    ax1.set_ylim(0, 1)
    ax2.set_ylim(0, 1)
    ax3.set_ylim(0, 1)
    ax3.set_xlim(19, 24)
    pl.tight_layout()
    ax1.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))


    # pl.show()
    # ax1.legend(loc=2)
    ax2.legend(loc=2)
    # ax3.legend(loc=2)
    pl.savefig("../document/figures/completeness_BAT.pdf", dpi="figure")
    pl.show()

if __name__ == '__main__':
    main()