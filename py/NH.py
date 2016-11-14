#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

from astropy.io import fits
import matplotlib.pyplot as pl
import seaborn; seaborn.set_style('ticks')
import numpy as np
import matplotlib
params = {
   'axes.labelsize': 8,
   'text.fontsize': 8,
   'legend.fontsize': 10,
   'xtick.labelsize': 10,
   'ytick.labelsize': 10,
   'text.usetex': False,
   'figure.figsize': [4.5, 4.5]
   }
matplotlib.rcParams.update(params)



def main():
    # Read in Buchner catalog
    grbcat = fits.open("../data/grbcat-annotated.fits")
    names = np.array(["".join(kk.split(" ")) for kk in grbcat[1].data.field("Name")])

    # Read in Optical afterglow (OA)
    OA = np.array([[ii, kk, ll, pp] for ii, kk, ll, pp in np.genfromtxt("../data/HIcoulmns.dat", dtype=None)])


    # print(grbcat[1].header)
    # print(grbcat[1].data)
    # exit()
    # Find indices of overlap
    idx = [ii for ii, kk in enumerate(names) if kk in OA[:, 0]]
    idx_OA = [ii for ii, kk in enumerate(OA[:, 0]) if kk in names[idx]]

    # Cut parameters in overlap
    OA_names, OA_z, OA_nH, OA_nHe = OA[idx_OA, 0], OA[idx_OA, 1].astype("float"), OA[idx_OA, 2].astype("float"), OA[idx_OA, 3].astype("float")

    names, t90, nH_gal, z = names[idx], grbcat[1].data.field("duration")[idx], grbcat[1].data.field("nhgal")[idx], grbcat[1].data.field("z")[idx]
    # nH, nH_l, hH_h, nH_std = grbcat[1].data.field("NH_sphere_mean")[idx], grbcat[1].data.field("NH_sphere_hi")[idx], grbcat[1].data.field("NH_sphere_lo")[idx], grbcat[1].data.field("NH_sphere_std")[idx]
    nH, nH_l, hH_h, nH_std = grbcat[1].data.field("NH_tbabs_mean")[idx], grbcat[1].data.field("NH_tbabs_hi")[idx], grbcat[1].data.field("NH_tbabs_lo")[idx], grbcat[1].data.field("NH_tbabs_std")[idx]
    Bx, Bx_std = grbcat[1].data.field("PhoIndex_sphere_mean")[idx], grbcat[1].data.field("PhoIndex_sphere_std")[idx]

    # Sorting lists to match bursts
    sort = np.argsort(names)
    sort_OA = np.argsort(OA_names)

    # for ii, kk in enumerate(OA_names[sort_OA]):
    #     print(kk, names[sort][ii])

    # Plot
    fig, ax = pl.subplots()

    x = np.arange(10, 30, 0.01)
    ax.plot(x, x)

    cmap = pl.get_cmap("plasma")

    sc = ax.scatter(nH[sort], OA_nH[sort_OA], c=OA_z[sort_OA], cmap=cmap, s=10)
    clb = pl.colorbar(sc)
    clb.set_label("Redshift")
    a, b, c = ax.errorbar(nH[sort], OA_nH[sort_OA], xerr=nH_std[sort], yerr=OA_nHe[sort_OA], fmt="none", marker=None, mew=0, alpha=0.5, lw=0.5)
    e_color = clb.to_rgba(OA_z[sort_OA])
    c[0].set_color(e_color)
    c[1].set_color(e_color)

    # from matplotlib.patches import Ellipse
    # for jj in np.arange(1, 2):
    #     for ii, (kk, ll) in enumerate(zip(nH[sort], OA_nH[sort_OA])):
    #         ax.add_artist(Ellipse((kk, ll), 2*jj*nH_std[sort][ii]/10, 2*jj*OA_nHe[sort_OA][ii], fill=False, linestyle='dashed', lw = 2, alpha = 1.0/(2.0*jj) , color=cmap(OA_z[sort_OA][ii]/max(OA_z[sort_OA]))))


    for ii, txt in enumerate(names[sort]):
        ax.annotate(txt, (nH[sort][ii], OA_nH[sort_OA][ii]), size="large")


    ax.annotate(r"Z/Z$_\odot$ = 1", xy=(20, 20.5), rotation=50)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # ax.spines['left'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.tick_params(axis='x', direction='out')
    ax.tick_params(axis='y', direction='out')
    # offset the spines
    for spine in ax.spines.values():
        spine.set_position(('outward', 3))
    # put the grid behind
    ax.set_axisbelow(True)

    ax.set_xlim((19.5, 23))
    ax.set_ylim((19, 23))
    ax.set_xlabel(r"X-ray derived column density N$_H$ [cm$^{-2}$]")
    ax.set_ylabel(r"Optically derived column density N$_H$ [cm$^{-2}$]")
    # pl.title("Optical vs X-ray derived column densities")

    pl.tight_layout()
    pl.show()
    pl.savefig("../document/figures/OAXR_NH.pdf", dpi="figure")


if __name__ == '__main__':
    main()