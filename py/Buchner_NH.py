#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

from astropy.io import fits
import matplotlib.pyplot as pl
import seaborn; seaborn.set_style('ticks')

def main():
    grbcat = fits.open("../data/grbcat-annotated.fits")
    print(grbcat[1].header)
    names = grbcat[1].data.field("Name")
    t90 = grbcat[1].data.field("duration")
    nH_gal = grbcat[1].data.field("nhgal")
    z = grbcat[1].data.field("z")
    nH = grbcat[1].data.field("NH_sphere_mean")
    nH_l, hH_h = grbcat[1].data.field("NH_sphere_hi"), grbcat[1].data.field("NH_sphere_lo")
    print(names, nH)
    # print(hH_h)
    # pl.errorbar(z, nH, yerr=[nH - nH_l, hH_h - nH], fmt=".k", capsize=0, elinewidth=0.5, ms=3)
    # pl.show()



if __name__ == '__main__':
    main()