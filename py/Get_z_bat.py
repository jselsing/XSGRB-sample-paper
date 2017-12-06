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
from astropy.cosmology import Planck15 as cosmo
import astropy.units as u

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
    Small script to get Swift redshift BAT .
    """
    # swift_table = pd.read_table("../data/allSwiftGRBs_NH.txt", delimiter="\t", dtype=None)
    # # print(swift_table)
    # name_s, z_s, NHX_s, NHXh_s, NHXl_s = swift_table["#GRB"].values, swift_table["Redshift"].values, swift_table["Ave_NH_PC"].values, swift_table["Ave_dNH_PC+"].values, swift_table["Ave_dNH_PC-"].values

    swift_table = pd.read_table("../data/grb_table_1511519199.txt", delimiter="\t", dtype=None)
    name, BAT, z = swift_table["GRB"].values, swift_table["BAT Fluence (15-150 keV) [10^-7 erg/cm^2]"].values, swift_table["Redshift"].values
    z_out = []
    idxs = []
    for ii, kk in enumerate(z):
      try:
        z_out.append(kk.split()[0])
        idxs.append(ii)
      except:
        pass

    name = name[idxs]
    BAT = filt_nan(BAT[idxs])
    z = filt_nan(z_out)

    L_bat = []
    for kk, ll, pp in list(zip(name, z, BAT)):
        dL = cosmo.luminosity_distance(ll).to(u.cm)
        BAT_flu = pp * 1e-7 * u.erg * u.cm**(-2)
        L_bat.append(BAT_flu * 4 * np.pi * dL**2 / (1 + ll)/u.erg)
    # L_bat = np.array(L_bat)


    np.savetxt("../data/zbatl.dat", list(zip(z, np.log10(L_bat))))





    L_bat = []
    z = np.arange(0, 15, 0.01)
    for ll in z:
        dL = cosmo.luminosity_distance(ll).to(u.cm)
        BAT_flu = 3 * 1e-8 * u.erg * u.cm**(-2)
        L_bat.append(BAT_flu * 4 * np.pi * dL**2 / (1 + ll)/u.erg)
    # L_bat = np.array(L_bat)


    np.savetxt("../data/zbatlswift.dat", list(zip(z, np.log10(L_bat))))





if __name__ == '__main__':
  main()