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
fbo = [3.2425, 2.8983, 1.2356, 3.9693, 1.38, 1.7102, 2.6147, 0.8278, 3.3467, 2.1995, 3.5328, 2.3, 4.0559, 3.9133, 1.5026, 2.1, 3.2213, 3.0749, 0.1257, 3.424, 1.92, 2.7108, 0.5428, 1.8836, 5.4636, 1.2622, 0.3463, 0.7578, 3.4344, 1.3145, 2.3521, 1.4965, 0.8397, 0.9705, 2.309, 2.0394, 3.6298, 2.4541, 2.1462, 5.2, 2.6918, 0.8227, 1.3308, 2.6419, 2.4274, 0.9382, 1.9492, 1.5119, 1.1014, 1.5457, 2.6892, 1.4171, 1.6403, 3.0368, 1.2322, 0.8454, 2.5914, 2.2045, 1.5042, 3.3604, 2.3739, 6.7, 0.6887, 1.6919, 0.7029, 3.6856, 3.2086, 1.5471, 2.9538, 1.0301, 2.433]

tough = [1.95, 1.4436, 2.8983, 2.7, 0.6528, 5.2, 0.606, 3.96855, 1.38, 2.5043, 2.61469, 1.434, 0.8278, 6.295, 3.3467, 2.5273, 2.1992, 2.4296, 1.059, 0.9364, 0.481, 3.5328, 0.03351, 5.11, 3.2213, 2.1357, 3.773, 3.0749, 0.125, 3.424, 1.92, 2.7108, 1.532, 0.5428, 1.9229, 1.8836, 0.937, 5.4636, 1.2622, 0.3463, 0.7578, 3.4344, 1.3145, 2.6208, 2.3521, 2.3384, 1.49594, 0.8397, 1.9591, 2.309, 2.0394, 3.6298, 2.4541, 4.406, 0.0889, 0.787, 2.17, 2.0858]

bat = [1.44, 2.9, 0.65, 0.61, 1.71, 2.2, 4.05, 3.91, 3.5, 0.13, 1.92, 1.88, 0.94, 5.47, 1.26, 0.35, 1.31, 2.09, 1.5, 1.35, 2.15, 0.82, 1.33, 0.94, 1.95, 1.1, 0.77, 1.4, 2.69, 1.64, 3.04, 2.59, 2.2, 0.69, 0.53, 2.51, 2.1, 2.26, 2.77, 1.55, 0.54, 3, 2.45, 1.24, 0.97, 1.71, 0.49, 1.06, 0.54, 2.106, 2.22, 1.613]

bat2 = [1.44, 2.9, 0.65, 0.61, 1.71, 2.2, 4.05, 3.91, 3.5, 0.13, 1.92, 1.88, 0.94, 5.47, 1.26, 0.35, 1.31, 2.09, 1.5, 4, 1.35, 2.15, 0.82, 1.33, 0.94, 1.95, 1.1, 0.77, 1.4, 2.69, 1.64, 3.04, 2.59, 2.2, 0.69, 0.53, 2.51, 2.1, 2.26, 2.77, 1.55, 4, 0.54, 3.5, 3, 2.45, 1.24, 0.97, 1.71, 0.49, 1.06, 0.54, 2.106, 2.22, 1.613]

sho = [1.95, 1.4436, 3.2425, 2.8983, 0.606, 3.9693, 1.7102, 2.6147, 1.434, 6.295, 4.9, 2.1995, 2.4296, 1.059, 3.5328, 0.785, 2.3393, 3.9122, 0.0331, 1.559, 1.5026, 5.11, 3.2213, 3.0749, 3.424, 2.7108, 1.532, 0.5428, 1.9229, 1.8836, 0.937, 5.467, 1.2622, 0.3463, 0.7578, 3.4344, 1.3145, 2.253, 2.088, 2.3521, 2.3384, 1.6295, 1.4959, 0.84, 2.0627, 1.9588, 0.82, 2.0865, 3.6298, 1.35, 2.1462, 2.452, 0.8227, 2.72, 2.0858, 2.6419, 2.4274, 2.0265, 0.9382, 1.78, 1.0301, 2.433, 1.1014, 0.767, 2.6892, 1.6403, 3.0368, 0.8454, 2.5914, 2.2045, 1.5042, 3.3604, 0.6887, 1.6919, 1.967, 3.8479, 0.9787, 2.58, 2.512, 2.0631, 2.26, 2.77, 3.375, 0.345, 1.608, 0.544, 4.109, 3.85, 1.266, 0.54, 2.452, 0.696, 1.24, 0.971, 2.752, 3.076, 0.49, 1.0633, 1.398, 0.542, 2.106, 1.44, 2.22, 2.09, 1.728]

def main():
    """
    Small script to produce sample completeness plots.
    """


    # Read in burst list
    # burst_table = np.array([[ii, kk, ll, pp, tt] for ii, kk, ll, pp, tt in np.genfromtxt("../data/burst_list.dat", dtype=None)])
    burst_table = pd.read_csv("../data/Burst list - master.csv")

    name, z, observed, OA, BAT, XRT, HI, pure = burst_table["GRB"].values, burst_table["Redshift"].values, burst_table["Observed"].values, burst_table["Afterglow"].values, burst_table["BAT Fluence (15-150 keV) [10^-7 erg/cm^2]"].values.astype("float"), burst_table["XRT 11 Hour Flux (0.3-10 keV) [10^-11 erg/cm^2/s]"].values.astype("float"), burst_table["XRT Column Density (NH) [10^21 cm^-2]"].values.astype("float"), burst_table["In pure sample"].values.astype("string")

    # print(pure == "Yes")
    l, m, h = np.nanpercentile(filt_nan(z[pure == "Yes"]), (14, 50, 86))
    mm = np.nanmean(filt_nan(z[pure == "Yes"]))
    print(mm, m, l - m, h - m)

    l, m, h = np.nanpercentile(fbo, (14, 50, 86))
    mm = np.nanmean(fbo)
    print(mm, m, l - m, h - m)

    l, m, h = np.nanpercentile(tough, (14, 50, 86))
    mm = np.nanmean(tough)
    print(mm, m, l - m, h - m)


    l, m, h = np.nanpercentile(bat2, (14, 50, 86))
    mm = np.nanmean(bat2)
    print(mm, m, l - m, h - m)

    l, m, h = np.nanpercentile(sho, (14, 50, 86))
    mm = np.nanmean(sho)
    print(mm, m, l - m, h - m)

    print(stats.ks_2samp(sho, filt_nan(z[pure == "Yes"])[~np.isnan(filt_nan(z[pure == "Yes"]))]))
    print(stats.ks_2samp(fbo, filt_nan(z[pure == "Yes"])[~np.isnan(filt_nan(z[pure == "Yes"]))]))
    print(stats.ks_2samp(tough, filt_nan(z[pure == "Yes"])[~np.isnan(filt_nan(z[pure == "Yes"]))]))
    print(stats.ks_2samp(bat2, filt_nan(z[pure == "Yes"])[~np.isnan(filt_nan(z[pure == "Yes"]))]))
    exit()
    BAT_o = BAT[observed == "Yes"]
    XRT_o = XRT[observed == "Yes"]
    HI_o = HI[observed == "Yes"]

    BAT_c = BAT#[observed == "No"]
    XRT_c = XRT#[observed == "No"]
    HI_c = HI#[observed == "No"]

    swift_table = pd.read_table("../data/grb_table_1511519199.txt", delimiter="\t", dtype=None)
    name_s, BAT_s, XRT_s, HI_s = swift_table["GRB"].values, swift_table["BAT Fluence (15-150 keV) [10^-7 erg/cm^2]"].values, swift_table["XRT 11 Hour Flux (0.3-10 keV) [10^-11 erg/cm^2/s]"].values, swift_table["XRT Column Density (NH) [10^21 cm^-2]"].values
    # print(len(name_s))
    # exit()
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
    sns.distplot(np.log10(1e-7*BAT_s), ax=ax1, kde=False, norm_hist=True, hist_kws={"histtype": "step", "alpha": 1, "linewidth": 2, "linestyle": "dashed", 'cumulative': True}, label=r"\textit{Swift}")
    sns.distplot(np.log10(1e-7*BAT_c), ax=ax1, kde=False, norm_hist=True, hist_kws={"histtype": "step", "alpha": 1, "linewidth": 2, "linestyle": "dashed", 'cumulative': True}, label="Statistical")
    sns.distplot(np.log10(1e-7*BAT_o), ax=ax1, kde=False, norm_hist=True, hist_kws={"histtype": "step", "alpha": 1, "linewidth": 2, "linestyle": "dashed", 'cumulative': True}, label="Observed")
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
    sns.distplot(np.log10(1e-11*XRT_s), ax=ax2, kde=False, norm_hist=True, hist_kws={"histtype": "step", "alpha": 1, "linewidth": 2, "linestyle": "dashed", 'cumulative': True}, label=r"\textit{Swift}")
    sns.distplot(np.log10(1e-11*XRT_c), ax=ax2, kde=False, norm_hist=True, hist_kws={"histtype": "step", "alpha": 1, "linewidth": 2, "linestyle": "dashed", 'cumulative': True}, label="Statistical")
    sns.distplot(np.log10(1e-11*XRT_o), ax=ax2, kde=False, norm_hist=True, hist_kws={"histtype": "step", "alpha": 1, "linewidth": 2, "linestyle": "dashed", 'cumulative': True}, label="Observed")



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
    sns.distplot(NHX_s, ax=ax3, kde=False, norm_hist=True, hist_kws={"histtype": "step", "alpha": 1, "linewidth": 2, "linestyle": "dashed", 'cumulative': True}, label=r"\textit{Swift}")
    sns.distplot(NHX_c, ax=ax3, kde=False, norm_hist=True, hist_kws={"histtype": "step", "alpha": 1, "linewidth": 2, "linestyle": "dashed", 'cumulative': True}, label="Statistical")
    sns.distplot(NHX_o, ax=ax3, kde=False, norm_hist=True, hist_kws={"histtype": "step", "alpha": 1, "linewidth": 2, "linestyle": "dashed", 'cumulative': True}, label="Observed")



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



    ax1.set_xlabel(r"$\log(15-150~\mathrm{keV}~\mathrm{Fluence}/\mathrm{erg}~\mathrm{cm}^{-2})$")
    ax2.set_xlabel(r"$\log(0.3-10~\mathrm{keV}~\mathrm{Flux}/\mathrm{erg}~\mathrm{cm}^{-2}~\mathrm{s}^{-1})$")
    ax3.set_xlabel(r"$\log(N_{\mathrm{HI, X}}/\mathrm{cm}^{-2})$")
    # ax1.set_ylabel(r"N")

    ax1.set_xlim(-8, -4)
    ax2.set_xlim(-18, -10.4)
    # ax1.set_xlim(-10, -3)
    # ax2.set_xlim(-20, -9)
    ax1.set_ylim(-0.1, 1.1)
    ax2.set_ylim(-0.1, 1.1)
    ax3.set_ylim(-0.1, 1.1)
    ax3.set_xlim(10, 22.9)
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