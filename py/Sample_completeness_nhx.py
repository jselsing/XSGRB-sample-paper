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
   'text.usetex': False,
   #'figure.figsize': [4.5, 7.281]
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
            pass
    return np.array(holder)


def main():
    """
    Small script to produce sample completeness plots.
    """


    # Read in burst list
    # burst_table = np.array([[ii, kk, ll, pp, tt] for ii, kk, ll, pp, tt in np.genfromtxt("../data/burst_list.dat", dtype=None)])
    burst_table = pd.read_csv("../data/Burst list - master.csv")

    name, z, observed, OA, NHX, NHXh, NHXl = burst_table["GRB"].values, burst_table["Redshift"].values, burst_table["Observed"].values, burst_table["Afterglow"].values, burst_table["Ave_NH_PC"].values.astype("float"), burst_table["Ave_dNH_PC+"].values.astype("float"), burst_table["Ave_dNH_PC-"].values.astype("float")

    print(z)

    NHX_o = NHX[observed == "Yes"]
    NHXh_o = NHXh[observed == "Yes"]
    NHXl_o = NHXl[observed == "Yes"]

    NHX_c = NHX[observed == "No"]
    NHXh_c = NHXh[observed == "No"]
    NHXl_c = NHXl[observed == "No"]

    swift_table = pd.read_table("../data/allSwiftGRBs_NH.txt", delimiter="\t", dtype=None)
    # print(swift_table)
    name_s, NHX_s, NHXh_s, NHXl_s = swift_table["#GRB"].values, swift_table["Ave_NH_PC"].values, swift_table["Ave_dNH_PC+"].values, swift_table["Ave_dNH_PC-"].values
    # Exclude sample bursts
    idx = [ii for ii, kk in enumerate(name_s) if kk not in name]

    NHX_s = filt_nan(NHX_s[idx])
    NHXh_s = filt_nan(NHXh_s[idx])
    NHXl_s = filt_nan(NHXl_s[idx])


    NHX_c, NHX_o, NHX_s = NHX_c[(NHX_c > 0)], NHX_o[(NHX_o > 0)], NHX_s[(NHX_s > 0)]
    NHXh_c, NHXh_o, NHXh_s = NHXh_c[(NHXh_c > 0)], NHXh_o[(NHXh_o > 0)], NHXh_s[(NHXh_s > 0)]
    NHXl_c, NHXl_o, NHXl_s = NHXl_c[(NHXl_c > 0)], NHXl_o[(NHXl_o > 0)], NHXl_s[(NHXl_s > 0)]


    # Rescale
    NHX_s = np.log10(1e22*NHX_s)
    NHX_c = np.log10(1e22*NHX_c)
    NHX_o = np.log10(1e22*NHX_o)
    # Plot
    fig, (ax1) = pl.subplots(ncols=1)



    # BAT fluence
    sns.distplot(NHX_s, ax=ax1, kde=False, norm_hist=True, hist_kws={"histtype": "step", "alpha": 1, "linewidth": 3, "linestyle": "dashed"}, label="Full Swift sample")
    sns.distplot(NHX_c, ax=ax1, kde=False, norm_hist=True, hist_kws={"histtype": "step", "alpha": 1, "linewidth": 3, "linestyle": "dashed"}, label="Complete sample")
    sns.distplot(NHX_o, ax=ax1, kde=False, norm_hist=True, hist_kws={"histtype": "step", "alpha": 1, "linewidth": 3, "linestyle": "dashed"}, label="Observed sample")
    print(len(NHX_c[~np.isnan(NHX_c)]))
    l, m, h = np.percentile(NHX_c, [16, 50, 84])
    print(m, m - l, h - m)
    print(len(NHX_o[~np.isnan(NHX_o)]))
    l, m, h = np.percentile(NHX_o, [16, 50, 84])
    print(m, m - l, h - m)
    print(len(NHX_s[~np.isnan(NHX_s)]))
    l, m, h = np.percentile(NHX_s, [16, 50, 84])
    print(m, m - l, h - m)

    print(stats.ks_2samp(NHX_s[~np.isnan(NHX_s)], NHX_c)[1])
    print(stats.ks_2samp(NHX_c, NHX_o)[1])
    print(stats.ks_2samp(NHX_s[~np.isnan(NHX_s)], NHX_o)[1])



    # for ax in [ax1]:
    #     ax.spines['top'].set_visible(False)
    #     ax.get_xaxis().tick_bottom()
    #     ax.get_yaxis().tick_left()
    #     ax.tick_params(axis='x', direction='out')
    #     ax.tick_params(axis='y', direction='out')

    #     # # offset the spines
    #     for spine in ax.spines.values():
    #         spine.set_position(('outward', 3))

    #     # put the grid behind
    #     ax.set_axisbelow(True)
        # ax.set_ylim(0, 1)




    ax1.set_xlabel(r"log(N$_H$) [cm$^2$]")
    ax1.set_ylabel(r"N")
    ax1.set_xlim((19, 23))
    pl.legend(loc=2)
    pl.tight_layout()
    ax1.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))


    # pl.show()
    pl.savefig("../document/figures/completeness_NHX.pdf", dpi="figure")
    pl.show()

if __name__ == '__main__':
    main()