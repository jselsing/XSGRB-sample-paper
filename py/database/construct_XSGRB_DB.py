#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb
import os, glob
import imp
import json

from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.table import Table, Column, vstack
from astropy.time import Time
from astropy.io import fits
from astropy import units as u

from linetools.spectra import io as lsio
from linetools import utils as ltu

from specdb.build.utils import chk_meta
from specdb.build import privatedb as pbuild
import specdb

def main():
    test = fits.open("/Users/jselsing/github/specdb/specdb/data/privateDB/XSGRBDB_ztbl.fits")
    # print(test[1].data.field)


    burst_list = np.genfromtxt("/Users/jselsing/github/XSGRB-sample-paper/data/burst_list.dat", dtype=None)
    names = [ii[0] for ii in burst_list][6:]
    redshifts = [float(ii[1]) for ii in burst_list][6:]
    z_source = [str("XSGRB")]*len(names)
    ra, dec, count = [0]*len(names), [0]*len(names), 0
    for kk, ll in zip(names, redshifts):
        bursts = glob.glob("../specdb/*"+kk+"*")

        try:
            hdr = fits.open(bursts[0])
            ra[count] = float(hdr[1].header["RA"])
            dec[count] = float(hdr[1].header["DEC"])
        except:
            pass
        count += 1

    prihdr = fits.Header()
    prihdu = fits.PrimaryHDU(header=prihdr)

    tbhdu = fits.BinTableHDU.from_columns(
        [fits.Column(name = str("RA"), format="D", array=ra),
         fits.Column(name = str("DEC"), format="D", array=dec),
         fits.Column(name = str("ZEM"), format="D", array=redshifts),
         fits.Column(name = str("ZEM_SOURCE"), format="A", array=z_source)])

    thdulist = fits.HDUList([prihdu, tbhdu])
    thdulist.writeto('/Users/jselsing/github/specdb/specdb/data/privateDB/XSGRBDB_ztbl.fits', clobber = True)

    # Read z table
    ztbl = Table.read(specdb.__path__[0]+'/data/privateDB/XSGRBDB_ztbl.fits')

    # Go
    tree = specdb.__path__[0]+'/data/privateDB'
    pbuild.mk_db('XSGRB', tree, '../XSGRB.hdf5', ztbl, fname=True)


if __name__ == '__main__':
    main()