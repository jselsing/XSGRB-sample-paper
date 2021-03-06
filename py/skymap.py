#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

from astropy.io import fits
import pandas as pd
import matplotlib; matplotlib.use('TkAgg')
import matplotlib.pyplot as pl
import healpy as hp
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
   'text.usetex': True,
   'figure.figsize': [9, 4.5]
   }
matplotlib.rcParams.update(params)
from astropy.coordinates import SkyCoord
from astropy.coordinates import Galactic, FK5
from astropy.coordinates import Angle, Latitude, Longitude
import astropy.units as u

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
    # Script to produce sky plots.
    """
    ebv_map = hp.read_map("../data/lambda_sfd_ebv.fits")
    A_Vmap = ebv_map*3.1
    A_Vmap_cut = A_Vmap.copy()
    A_Vmap_cut[A_Vmap >= 0.5] = 0
    ax = hp.mollview(A_Vmap_cut, max=0.5, title="", cbar=False)

    # , cbar=False


    fig = pl.gcf()
    ax = pl.gca()
    image = ax.get_images()[0]
    cmap = fig.colorbar(image, ax=ax, location="bottom", pad = 0.02, label=r"$A_V$", shrink=0.5)
    hp.graticule()

    # A_Vmap_cut = A_Vmap.copy()
    # A_Vmap_cut[A_Vmap <= 0.5] = np.nan
    # ax = hp.mollview(A_Vmap_cut, title="", cbar=False, norm="hist", cmap="viridis", ax=ax)

    # deep=["#4C72B0", "#55A868", "#C44E52",
    #       "#8172B2", "#CCB974", "#64B5CD"],

    # Overplot swift positions
    # Read in master sample
    burst_table = pd.read_csv("../data/Burst list - master.csv")
    name, z, observed, OA, = burst_table["GRB"].values, burst_table["Redshift"].values, burst_table["Observed"].values, burst_table["Afterglow"].values


    hasz = []
    for ii, kk in enumerate(z):
        try:
            if not np.isnan(float(kk)):
                hasz.append(name[ii])
        except:
            pass

    # Read in observed bursts
    burst_table = pd.read_csv("../data/Burst list - observed.csv")
    name_o = np.array([ii[3:] for ii in burst_table["GRB"].values])


    for ii, kk in enumerate(burst_table["z"].values):
        try:
            if not np.isnan(float(kk)):
                hasz.append(name_o[ii])
        except:
            pass

    # Read in Swift table
    swift_table = pd.read_table("../data/grb_table_1517061198.txt", delimiter="\t", dtype=None)
    name_s, RA, DEC, redshift, bat, xrt, nh = swift_table["GRB"].values, pd.to_numeric(swift_table["BAT RA (J2000)"], errors='coerce').values, pd.to_numeric(swift_table["BAT Dec (J2000)"], errors='coerce').values, swift_table["Redshift"].values, swift_table["BAT Fluence (15-150 keV) [10^-7 erg/cm^2]"].values, swift_table["XRT 11 Hour Flux (0.3-10 keV) [10^-11 erg/cm^2/s]"].values, swift_table["XRT Column Density (NH) [10^21 cm^-2]"].values

    # nammmm = ["170111A", "170113A", "170127A", "170127B", "170202A", "170208A", "170208B"]

    # idx = [ii for ii, kk in enumerate(name_s) if kk in nammmm]

    # for ii, kk in enumerate(nammmm):
    #     print(bat[idx][ii])

    # for ii, kk in enumerate(nammmm):
    #     print(xrt[idx][ii])

    # for ii, kk in enumerate(nammmm):
    #     print(nh[idx][ii])

    # exit()
    count = 0
    for ii, kk in enumerate(name_s):
        hej = ''.join(filter(lambda x: x.isdigit(), kk))
        # print(hej, ii)

        try:
            if float(hej) < 170331:
                count += 1
        except:
            pass
    print(count)

    l, b = np.ones_like(RA), np.ones_like(DEC)
    for ii, (kk, ll) in enumerate(zip(RA, DEC)):
        c = SkyCoord(ra=kk, dec=ll, frame='fk5', unit='deg')
        l[ii] = c.galactic.l.degree
        b[ii] = c.galactic.b.degree

    # Find bursts with redshifts
    holder = [0]*len(redshift)
    for ii, kk in enumerate(redshift):
        try:
            holder[ii] = ''.join(c for c in kk if c.isdigit())
            holder[ii] = 1
        except:
            holder[ii] = np.nan

    for ii, kk in enumerate(name_s):
        if kk in hasz:
            holder[ii] = 1

    # Create mask of bad values and redshifts and apply
    mask_z = (~(np.isnan(l) | np.isnan(b))) & ~np.isnan(holder)
    name_s_z = name_s[mask_z]
    l_z = l[mask_z]
    b_z = b[mask_z]


    # print(l_z, b_z)
    # Scatterplot with bursts position
    print(len(l_z))
    hp.visufunc.projscatter(l_z, b_z, facecolors='none', edgecolor="black", marker = "*", s = 50, lonlat=True)

    # Get indices of sample bursts and observed bursts
    idx = [ii for ii, kk in enumerate(name_s_z) if kk in name]
    idx_o = [ii for ii, kk in enumerate(name_s_z) if kk in name_o]
    idx_b = [ii for ii, kk in enumerate(name_s_z) if kk in name_o and kk in name]
    # Scatter
    print(len(l_z[idx]))
    print(len(l_z[idx_o]))
    print(len(l_z[idx_b]))

    hp.visufunc.projscatter(l_z[idx], b_z[idx], edgecolor="none", facecolors="#4C72B0", s = 100, marker = "*", lonlat=True)
    hp.visufunc.projscatter(l_z[idx_o], b_z[idx_o], facecolor="none", edgecolor="#C44E52", s = 110, lw=1, marker = "*", lonlat=True)

    # Create mask of bad values and redshifts and apply
    mask_noz = (~(np.isnan(l) | np.isnan(b))) & np.isnan(holder)
    name_s_noz = name_s[mask_noz]
    l_noz = l[mask_noz]
    b_noz = b[mask_noz]

    # Scatterplot with bursts position
    print(len(l_noz))
    hp.visufunc.projscatter(l_noz, b_noz, facecolors='none', edgecolor="black", s = 30, lonlat=True)

    # Get indices of sample bursts and observed bursts
    idx = [ii for ii, kk in enumerate(name_s_noz) if kk in name]
    idx_o = [ii for ii, kk in enumerate(name_s_noz) if kk in name_o]
    idx_b = [ii for ii, kk in enumerate(name_s_noz) if kk in name_o and kk in name]
    # print(name_s_noz[idx_o])
    # Scatter
    print(len(l_noz[idx]))
    print(len(l_noz[idx_o]))
    print(len(l_noz[idx_b]))
    hp.visufunc.projscatter(l_noz[idx], b_noz[idx], edgecolor="none", facecolors="#4C72B0", s = 40, lonlat=True)
    hp.visufunc.projscatter(l_noz[idx_o], b_noz[idx_o], facecolor="none", edgecolor="#C44E52", s = 50, lw=1, lonlat=True)
    # exit()




    # Save figure for tex
    pl.savefig("../document/figures/skymap.pdf", dpi="figure")
    pl.show()

if __name__ == '__main__':
    main()