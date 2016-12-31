#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

from astropy.io import fits
import pandas as pd
import matplotlib.pyplot as pl
import healpy as hp
import seaborn as sns; sns.set_style('ticks')
import numpy as np
from scipy import stats
import corner

params = {
   'axes.labelsize': 10,
   'font.size': 10,
   'legend.fontsize': 10,
   'xtick.labelsize': 10,
   'ytick.labelsize': 10,
   'text.usetex': False,
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
    name, z, NH, NHe = np.array([ii[3:] for ii in burst_table["GRB"].values]), burst_table["z"].values, burst_table["N_H"].values, burst_table["N_H_err"].values

    fig, ax = pl.subplots()
    # gs = pl.GridSpec(6, 6)

    # ax = fig.add_subplot(gs[1:, :-1])
    # ax_marg_x = fig.add_subplot(gs[0, :-1], sharex=ax)
    # ax_marg_y = fig.add_subplot(gs[1:, -1], sharey=ax)
    g = sns.JointGrid(x=Noter_z, y=Noter_NH, ax=ax)
    d = sns.JointGrid(x=z, y=NH, ax=ax)
    # corner.hist2d(Noter_z, Noter_NH, bins=100, ax=ax, smooth=2, color="#C44E52")
    # ax.scatter(Noter_z, Noter_NH)
    ax.scatter(z, NH, color="#4C72B0")
    # sns.distplot(NH, hist=True, kde=False, ax=ax_marg_y, hist_kws={"bins": 20, "color": "#4C72B0", "gridsize": 200, "alpha": 0.2},
                 # vertical=True, axlabel=False)



    # Save figure for tex
    pl.savefig("../document/figures/NH_dist.pdf", dpi="figure")
    pl.show()

if __name__ == '__main__':
    main()