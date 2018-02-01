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

params = {
   'axes.labelsize': 10,
   'font.size': 10,
   'legend.fontsize': 10,
   'xtick.labelsize': 10,
   'ytick.labelsize': 10,
   'text.usetex': True,
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


def filt_nan(input_array, fill_value=np.nan):
    """
    Small functino to remove values in array, not convertible to float. Return the filtered array
    """
    holder = []
    for ii in input_array:
        try:
            holder.append(float(ii))
        except ValueError:
            try:
                holder.append(float(ii.split("=")[1]))
            except IndexError:
                try:
                    holder.append(float(ii))
                except:
                    holder.append(fill_value)

    return np.array(holder)


def main():
    """
    # Script to produce beta OX.
    """
    # Read in burst list
    burst_table = pd.read_csv("../data/Burst list - CSV_observed.csv")
    name, sample, dt, acqmag, X_Fnu  = burst_table["GRB"].values, burst_table["In Sample"].values, burst_table["Follow-up delay"].values, burst_table["Acquisition mag"].values, burst_table["XRT flux density"].values
    name = np.array([ii[3:] for ii in name])
    # print(name)
    acqmag = filt_nan(acqmag, fill_value=24)
    idx_limit = (acqmag == 24)


    X_Fnu = filt_nan(X_Fnu)
    # Optical flux in microjansky
    O_Fnu = 10**((23.9 - acqmag)/2.5)

    O_nu0 = 3e8 / 658e-9 # Convert R-band rest to Hertz
    E_Xrt = ((10. - 0.3) / 2.)*1e3 # Middle energy of XRT in eV
    X_nu0 = E_Xrt / 4.135667516e-15 # Convert eV to Hertz
    # Calculate betaOX
    betaOX = np.log10(O_Fnu/X_Fnu)/np.log10(X_nu0/O_nu0)
    # Mask bad values and reverse order
    mask = np.isnan(betaOX) | (sample == "No")
    betaOX = betaOX[~mask][::-1]
    name = name[~mask][::-1]
    idx_limit = idx_limit[~mask][::-1]

    # for ii, kk in enumerate(betaOX):
    #     print(name[ii], kk)
    # exit()
    fynbo_table = np.genfromtxt("../data/Fynbo2009betaOX.dat", dtype=None)
    f_name = fynbo_table[:, 0]
    fynbo_table[:, 1] = [ii.replace('<', '') for ii in fynbo_table[:, 1]]
    f_betaOX = filt_nan(fynbo_table[:, 1])
    f_name = f_name[~np.isnan(f_betaOX)]
    f_betaOX = f_betaOX[~np.isnan(f_betaOX)]


    # Get Swift HI values
    swift_table = pd.read_table("../data/allSwiftGRBs_NH.txt", delimiter="\t", dtype=None)
    name_s, z_s, HI_s, HI_sh, HI_sl = swift_table["#GRB"].values, swift_table["Redshift"].values, swift_table["Ave_NH_PC"].values, swift_table["Ave_dNH_PC+"].values, -1*swift_table["Ave_dNH_PC-"].values
    # print(HI_s)


    # # Read in Buchner catalog
    # grbcat = fits.open("../data/grbcat-annotated.fits")
    # # print(grbcat[1].data.field)
    # name_s = np.array(["".join(kk.split(" "))[3:] for kk in grbcat[1].data.field("Name")])
    # HI_s, HI_sl, HI_sh, nH_std, z_s = grbcat[1].data.field("NH_sphere_mean"), grbcat[1].data.field("NH_sphere_hi"), grbcat[1].data.field("NH_sphere_lo"), grbcat[1].data.field("NH_sphere_std"), grbcat[1].data.field("z")
 


    # Find values with betaOX
    idx = [ii for ii, kk in enumerate(name_s) if kk in name and ~np.isnan(z_s[ii])]
    idx_2 = [ii for ii, kk in enumerate(name) if kk in name_s[idx]]
    f_idx = [ii for ii, kk in enumerate(name_s) if kk in f_name and ~np.isnan(z_s[ii])]
    f_idx_2 = [ii for ii, kk in enumerate(f_name) if kk in name_s[f_idx]]


    # HI = filt_nan(HI_s[idx], fill_value=0)
    # HIh = filt_nan(HI_sh[idx], fill_value=0)
    # HIl = filt_nan(HI_sl[idx], fill_value=0)
    # f_HI = filt_nan(HI_s[f_idx], fill_value=0)
    # f_HI[np.isnan(f_HI)] = 0

    HI = np.log10(1e22*filt_nan(HI_s[idx], fill_value=0))
    HIh = np.log10(1e22*filt_nan(HI_s[idx] + HI_sh[idx], fill_value=0))
    HIl = np.log10(1e22*filt_nan(HI_s[idx] - HI_sl[idx], fill_value=0))
    f_HI = np.log10(1e22*filt_nan(HI_s[f_idx], fill_value=0))
    f_HI[np.isnan(f_HI)] = 0

    # Plot new values


    print(stats.pearsonr(HI, betaOX[idx_2]))
    print(stats.spearmanr(HI, betaOX[idx_2]))
    print(stats.pointbiserialr(HI, betaOX[idx_2]))
    print(stats.kendalltau(HI, betaOX[idx_2]))
    
    # exit()


    g = sns.JointGrid(x=HI, y=betaOX[idx_2], xlim = (19.5, 23), ylim = (0, 1.3), space=0)
    color = sns.color_palette()[0]
    g = g.plot_marginals(sns.distplot, hist=True, kde=False, norm_hist=True, color=color)
    idx_limit = idx_limit[idx_2]
    g.x = HI[~idx_limit]
    g.y = betaOX[idx_2][~idx_limit]

    np.savetxt("betaOX.dat", list(zip(name_s[idx][~idx_limit] ,HI[~idx_limit], betaOX[idx_2][~idx_limit])), fmt = '%s %s %s')
    # exit()

    g = g.plot_joint(pl.scatter, color=color, label="This work")

    g.x = HI[idx_limit]
    g.y = betaOX[idx_2][idx_limit]
    g = g.plot_joint(pl.scatter, color=color, marker=r'$\downarrow$', s=150, lw=1.0)



    for ii, kk in enumerate(HI):
        g.y = np.array([betaOX[idx_2][ii], betaOX[idx_2][ii]])

        g.x = np.array([HIl[ii], HIh[ii]])
        g = g.plot_joint(pl.plot, color=color, lw=0.75)


    # Plot values from Fynbo et al. 2009
    color = sns.color_palette()[2]
    color_rgb = mpl.colors.colorConverter.to_rgb(color)
    colors = [sns.set_hls_values(color_rgb, l=l) for l in np.linspace(1, 0, 12)]
    cmap = sns.blend_palette(colors, as_cmap=True)

    g.x = f_HI
    g.y = f_betaOX[f_idx_2]

    # print(f_betaOX[f_idx_2], f_HI)

    g = g.plot_marginals(sns.distplot, hist=False, color=color)
    g = g.plot_joint(sns.kdeplot, cmap=cmap, label="Fynbo et al. 2009")
    g.x = 1
    g.y = 1
    g = g.plot_joint(pl.plot, color=color, label="Fynbo et al. 2009")
    ax = pl.gca()
    ax.axhline(0.5, color="black", linestyle="dashed", alpha=0.5)
    ax.annotate(r"$\beta_{OX} = 0.5$", (19.6, 0.45))
    g.set_axis_labels(r"$\log(N_{\mathrm{HI, X}}/\mathrm{cm}^{-2})$", r"Darkness [$\beta_{OX}$]")
    pl.tight_layout()

    # print(stats.ks_2samp(betaOX[idx_2], f_betaOX[f_idx_2]))
    # print(len(betaOX[idx_2]))
    # l, m, h = np.percentile(betaOX[idx_2], [16, 50, 84])
    # print(m, m - l, h - m)
    # print(len(f_betaOX[f_idx_2]))
    # l, m, h = np.percentile(f_betaOX[f_idx_2], [16, 50, 84])
    # print(m, m - l, h - m)

    # print(stats.ks_2samp(HI, f_HI))
    # print(len(HI))



    print("Above 0.5")
    l, m, h = np.percentile(HI[betaOX[idx_2] >= 0.5], [16, 50, 84])
    print(HI[betaOX[idx_2] >= 0.5])
    print(m, m - l, h - m)
    print("Below 0.5")
    l, m, h = np.percentile(HI[betaOX[idx_2] < 0.5], [16, 50, 84])
    print(HI[betaOX[idx_2] >= 0.5])
    print(m, m - l, h - m)

    #Fraction of dark bursts
    print(len(HI[betaOX[idx_2] < 0.5])/(len(HI[betaOX[idx_2] < 0.5]) + len(HI[betaOX[idx_2] >= 0.5])))
    print(stats.ks_2samp(HI[betaOX[idx_2] < 0.5], HI[betaOX[idx_2] >= 0.5]))
    print(len(HI))
    # exit()
    # print(len(f_HI))
    # l, m, h = np.percentile(f_HI, [16, 50, 84])
    # print(m, m - l, h - m)


    # Save figure for tex
    pl.legend()
    pl.savefig("../document/figures/betaOX.pdf", dpi="figure")
    pl.show()

if __name__ == '__main__':
    main()