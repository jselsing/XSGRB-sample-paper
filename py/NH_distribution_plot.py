#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Plotting
import matplotlib; matplotlib.use('TkAgg')
import matplotlib.pyplot as pl
import seaborn as sns; sns.set_style('ticks')
import matplotlib as mpl

from astropy.io import fits
import pandas as pd
# import matplotlib.pyplot as pl
import healpy as hp
# import seaborn as sns; sns.set_style('ticks')
import numpy as np
from scipy import stats


params = {
   'axes.labelsize': 10,
   'font.size': 10,
   'legend.fontsize': 10,
   'xtick.labelsize': 10,
   'ytick.labelsize': 10,
   'text.usetex': True ,
   'figure.figsize': [7.281, 4.5]
   }

def main():
    """
    # Script to produce NH distribution plot.
    """

    # Read DLA NH distribution from Noterdaeme et al. 2012b
    Noter_NH, Noter_z = np.zeros(12081), np.zeros(12081)
    count = 0
    with open("../data/J_A+A_547_L1/table3.dat") as fp:
        for line in fp:
            Noter_NH[count] = float(line[82:87])
            Noter_z[count] = float(line[51:56])
            count += 1

    # Read in GRB NH
    burst_table = pd.read_csv("../data/Burst list - HI columns.csv")
    name, z, NH, NHe, pub = np.array([ii[3:] for ii in burst_table["GRB"].values]), burst_table["z"].values, burst_table["N_H"].values, burst_table["N_H_err"].values, burst_table["Sample"].values
    XSGRB = (pub == "Selsing2017")
    # print(np.sum((~XSGRB).astype("int")))
    # print(len(name))
    # exit()

    # sns.set_palette("muted")
    fig, ax1 = pl.subplots()
    g = sns.jointplot(x=Noter_z, y=Noter_NH, kind="hex", stat_func=None, xlim = (1, 7), ylim = (20, 23), marginal_kws = {"kde": False, "norm_hist": True}, label="Noterdaeme et al. 2012b", space=0, color = "#4C72B0")

    # Plot old distribution
    NH_old = NH[~XSGRB]
    z_old = z[~XSGRB]
    g.x = z_old[NH_old > 20]
    g.y = NH_old[NH_old > 20]
    # g = g.plot_joint(pl.scatter, color="#55A868")
    g = g.plot_marginals(sns.distplot, hist=False, color="#55A868", bins = 25, norm_hist=True, kde_kws={"shade": True, "bw": 0.2})
    g = g.plot_joint(sns.kdeplot, color="#55A868", label="Tanvir et al. 2017")
    # pl.legend()
    # Plot new distribution
    NH_new = NH[XSGRB]
    z_new = z[XSGRB]
    g.x = z_new[NH_new > 20]
    g.y = NH_new[NH_new > 20]
    g = g.plot_marginals(sns.distplot, kde=False, color="#C44E52", bins = 25, norm_hist=True)
    g = g.plot_joint(pl.scatter, color="#C44E52", label="This work")
    g.x = 1
    g.y = 1
    g = g.plot_joint(pl.plot, color="#55A868", label="Tanvir et al. 2017")
    # print(stats.ks_2samp(NH_new, NH_old))
    # print(stats.ks_2samp(z_new, z_old))
    # l, m, h = np.percentile(NH_new, [16, 50, 84])
    # print(m, m - l, h - m)
    # l, m, h = np.percentile(NH_old, [16, 50, 84])
    # print(m, m - l, h - m)
    # l, m, h = np.percentile(z_new, [16, 50, 84])
    # print(m, m - l, h - m)
    # l, m, h = np.percentile(z_old, [16, 50, 84])
    # print(m, m - l, h - m)

    g.set_axis_labels("Redshift", r"$\log(N_{\mathrm{HI}}/\mathrm{cm}^{-2})$")

    pl.tight_layout()




    # Save figure for tex
    ax = pl.gca()
    handle, label = ax.get_legend_handles_labels()
    # from matplotlib.lines import Line2D     
    l1 = mpl.lines.Line2D([1], [1], color = "#4C72B0")
    l2 = mpl.lines.Line2D([1], [1], color = "#55A868")
    l3 = mpl.lines.Line2D([1], [1], color = "#C44E52")
    # print(handle, label)
    # pl.legend([l1, l2, l3], label)

    leg = pl.legend([l1, l2, l3], label, loc=1, handlelength=1)
    # leg = ax1.legend()
    # set the linewidth of each legend object
    for legobj in leg.legendHandles:
        legobj.set_linewidth(5.0)
        # legobj.handlelength(2.0)
    pl.savefig("../document/figures/NH_dist.pdf")
    # pl.show()

if __name__ == '__main__':
    main()