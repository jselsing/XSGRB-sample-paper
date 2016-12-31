#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

from astropy.io import fits

import glob
import numpy as np
import matplotlib.pyplot as pl
import os

def load_array(data_array):

    # Check is reposnse has been written based on value in sixt column
    thomas = False

    if np.shape(data_array)[1] == 3:
        thomas = True


    if thomas:
        return data_array[:, 0], data_array[:, 1], data_array[:, 2], np.zeros_like(data_array[:, 0])

    has_resp = False
    if np.nanmean(data_array[:, 6]) < 1e-10:
        has_resp = True



    if has_resp:
        try:
            tellcorr = data_array[:, 8]
            tellcorr[np.isnan(tellcorr)] = 1
            slitcorr = np.median(data_array[:, 7])
            return data_array[:, 1], data_array[:, 2]*slitcorr, data_array[:, 3]*slitcorr, data_array[:, 4], tellcorr
        except:
            slitcorr = np.median(data_array[:, 7])
            return data_array[:, 1], data_array[:, 2]*slitcorr, data_array[:, 3]*slitcorr, data_array[:, 4]
    elif not has_resp:
        try:
            tellcorr = data_array[:, 7]
            tellcorr[np.isnan(tellcorr)] = 1
            slitcorr = np.median(data_array[:, 6])
            return data_array[:, 1], data_array[:, 2]*slitcorr, data_array[:, 3]*slitcorr, data_array[:, 4], tellcorr
        except:
            slitcorr = np.median(data_array[:, 6])
            return data_array[:, 1], data_array[:, 2]*slitcorr, data_array[:, 3]*slitcorr, data_array[:, 4]



def main():
    """
    Script to combine extractions, telluric correction and continua in a single fits-file for ingestion into specdb
    """

    extract_dir = "../data/"
    dirs = [ii for ii in glob.glob(extract_dir+"*") if "GRB" in ii.split("/")[-1]]
    burst_names = [ii.split("/")[-1] for ii in dirs]

    for kk, ll in zip(burst_names, dirs):
        # if not "GRB160117B" in ll:
        #     continue
        print(kk)
        fitsfiles = glob.glob(ll+"/*")
        OBs = list(set([ii.split("/")[-1][:3] for ii in fitsfiles]))
        arms = list(set([ii.split("/")[-1][3:6] for ii in fitsfiles]))

        for tt in OBs:
            for yy in arms:
                outname = "../specdb/"+kk+"_"+tt+yy+".fits"
                if not os.path.exists(outname):
                    # Finding extractions

                    data_file = [ii for ii in glob.glob(extract_dir+kk+"/*1D*") if tt in ii and yy in ii and "norm" not in ii]


                    if len(data_file) == 0:
                        print("Empty dataset ... Skipping")
                        continue

                    # print(data_file)
                    try:
                        wl, flux, error, bp, tellcorr = load_array(np.genfromtxt(data_file[0]))
                    except:
                        print("No telluric correction found ...")
                        wl, flux, error, bp = load_array(np.genfromtxt(data_file[0]))
                        tellcorr = np.ones_like(wl)

                    # Finding normalisations
                    norm_path = [ii for ii in glob.glob(extract_dir+kk+"/*norm.npy*") if tt in ii and yy in ii][0]

                    norm = np.load(norm_path)
                    cont, cont_err = norm[5, :], norm[4, :]

                    # Open correct fits based on arm and OB to get headers
                    if len([ii for ii in glob.glob(extract_dir+kk+"/*2D*") if tt in ii and yy in ii and "norm" not in ii]) == 0:
                        print("No 2D-file for", tt, yy)
                        continue

                    file = fits.open([ii for ii in glob.glob(extract_dir+kk+"/*2D*") if tt in ii and yy in ii and "norm" not in ii][0])


                    prihdu = fits.PrimaryHDU(header=file[0].header)

                    wl = wl/(1 + file[0].header["HIERARCH ESO QC VRAD BARYCOR"]/3e5)
                    # print(1 + (file[0].header["HIERARCH ESO QC VRAD BARYCOR"]/3e5))
                    # exit()
                    nelem = len(wl)
                    tbhdu = fits.BinTableHDU.from_columns(
                        [fits.Column(name = "WAVE", format="D", unit="angstrom", array=wl),
                         fits.Column(name = "FLUX", format="D", unit="erg cm**(-2) s**(-1) angstrom**(-1)", array=flux),
                         fits.Column(name = "ERR", format="D", unit="erg cm**(-2) s**(-1) angstrom**(-1)", array=error),
                         fits.Column(name = "QUAL", format="J", unit="", array=bp),
                         fits.Column(name = "CONTINUUM", format="D", unit="erg cm**(-2) s**(-1) angstrom**(-1)", array=cont),
                         fits.Column(name = "CONTINUUM_ERR", format="D", unit="erg cm**(-2) s**(-1) angstrom**(-1)", array=cont_err),
                         fits.Column(name = "TELL_CORR", format="D", unit="", array=tellcorr)])

                    tbhdu.header["OBJECT"] = file[0].header["OBJECT"]
                    tbhdu.header["TITLE"] = kk
                    tbhdu.header["RA"] = file[0].header["RA"]
                    tbhdu.header["DEC"] = file[0].header["DEC"]
                    tbhdu.header["TELAPSE"] = 4*file[0].header["EXPTIME"]
                    tbhdu.header["MJD-OBS"] = file[0].header["MJD-OBS"]
                    tbhdu.header["TDMIN1"] = str(min(wl))
                    tbhdu.header["TDMAX1"] = str(max(wl))
                    tbhdu.header["SPEC_BW"] = str(max(wl) - min(wl))
                    tbhdu.header["NELEM"] = str(nelem)
                    thdulist = fits.HDUList([prihdu, tbhdu])

                    outname = "../specdb/"+kk+"_"+tt+yy+".fits"
                    thdulist.writeto(outname, clobber = True)
                    # exit()

if __name__ == '__main__':
    main()