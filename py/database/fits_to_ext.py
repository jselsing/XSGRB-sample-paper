#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

import numpy as np
import glob
from astropy.io import fits
import matplotlib.pyplot as pl


def main():
    """
    Small script to convert BinTableHDU to extensions to use with IRAF
    """

    files = glob.glob("../specdb/*")

    # Loop through files
    for ii in files:
        outname = ii[:3] + "fits" + ii[9:]
        bintab = fits.open(ii)
        header = bintab[0].header

        wl, flux, error, bpmap, cont, tell = bintab[1].data.field("WAVE"), bintab[1].data.field("FLUX"), bintab[1].data.field("ERR"), bintab[1].data.field("QUAL"), bintab[1].data.field("CONTINUUM"), bintab[1].data.field("TELL_CORR")

        minlam, maxlam, dlam = min(wl), max(wl), np.median(np.diff(wl))
        header["CUNIT1"] = "angstrom"
        header["CRVAL1"] = minlam
        header["CDELT"] = dlam
        header["CD1_1"] = dlam

        f = fits.PrimaryHDU(flux, header= header)
        e = fits.ImageHDU(error, name ="ERR", header= header)
        bp = fits.ImageHDU(bpmap, name ="QUAL", header= header)
        cont = fits.ImageHDU(cont, name ="CONTINUUM", header= header)
        tell = fits.ImageHDU(tell, name ="TELL_CORR", header= header)
        hdulist = fits.HDUList([f, e, bp, cont, tell])

        for kk in xrange(1, 5):
            hdulist[kk].header["XTENSION"]
            hdulist[kk].header["PCOUNT"]
            hdulist[kk].header["GCOUNT"]
            hdulist[kk].header["PCOUNT"]

        hdulist.writeto(outname, clobber=True)
        # exit()




    pass




if __name__ == '__main__':
    main()