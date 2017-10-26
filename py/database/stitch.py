#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

import sys
# sys.path.append("/Users/jselsing/github/XSGRB_reduction_scrips/py/")
import stitch_arms, util

from glob import glob
import numpy as np
import matplotlib.pyplot as pl
from scipy import optimize
from astropy.io import fits
from scipy.interpolate import interp1d
import os

def pow(x, y, z):
    return y * x ** z


def main():
    """
    Script to loop though bursts and stitch arms
    """
    data_dir = "../specdb/"
    bursts = [ii.split("/")[-1] for ii in glob("../data/*")]

    for kk in bursts:
        # if not kk == "GRB161023A":
        #     continue

        files = glob(data_dir+kk+"*")

        UVB_path = [ll for ll in files if "UVB" in ll]
        VIS_path = [ll for ll in files if "VIS" in ll]
        NIR_path = [ll for ll in files if "NIR" in ll]


        for tt, ll in enumerate(UVB_path):
            # if "GRB100418A" in ll:
            outname = (UVB_path[tt].split("/")[-1])[:-8]
            if not os.path.exists("../stitched/"+outname+"_stitched_spectrum.dat"):
                print("Stitching files:")
                print(UVB_path[tt])
                print(VIS_path[tt])
                print(NIR_path[tt])
                # outname = (UVB_path[tt].split("/")[-1])[:-8]
                print(outname)

                # Load synthetic transmission
                trans = np.genfromtxt("static/trans.dat")
                wl_trans = 10*trans[:, 0] # In micron
                trans = trans[:, 1]
                f = interp1d(wl_trans, trans, bounds_error=False, fill_value=np.nan)

                # Load files
                UVB = fits.open(UVB_path[tt])
                UVB_wl, UVB_flux, UVB_error, UVB_bp, UVB_cont, UVB_tell = UVB[1].data.field("WAVE"), UVB[1].data.field("FLUX"),  UVB[1].data.field("ERR"), UVB[1].data.field("QUAL"), UVB[1].data.field("CONTINUUM"), UVB[1].data.field("TELL_CORR")

                VIS = fits.open(VIS_path[tt])
                VIS_wl, VIS_flux, VIS_error, VIS_bp, VIS_cont, VIS_tell = VIS[1].data.field("WAVE"), VIS[1].data.field("FLUX"),  VIS[1].data.field("ERR"), VIS[1].data.field("QUAL"), VIS[1].data.field("CONTINUUM"), VIS[1].data.field("TELL_CORR")
                if np.mean(VIS_tell) == 1:
                    # Convolve to observed grid
                    VIS_tell = 1/f(VIS_wl)

                NIR = fits.open(NIR_path[tt])
                NIR_wl, NIR_flux, NIR_error, NIR_bp, NIR_cont, NIR_tell = NIR[1].data.field("WAVE"), NIR[1].data.field("FLUX"),  NIR[1].data.field("ERR"), NIR[1].data.field("QUAL"), NIR[1].data.field("CONTINUUM"), NIR[1].data.field("TELL_CORR")
                if np.mean(NIR_tell) == 1:
                    # Convolve to observed grid
                    NIR_tell = 1/f(NIR_wl)

                # Construct lists
                UVB_mask = (UVB_wl > 3200) & (UVB_wl < 5600)
                waves = [UVB_wl[UVB_mask], VIS_wl, NIR_wl]
                flux = [(UVB_flux/UVB_cont)[UVB_mask], VIS_tell*VIS_flux/VIS_cont, NIR_tell*NIR_flux/NIR_cont]
                error = [(UVB_error/UVB_cont)[UVB_mask], VIS_tell*VIS_error/VIS_cont, NIR_tell*NIR_error/NIR_cont]
                bpmap = [UVB_bp[UVB_mask], VIS_bp, NIR_bp]


                # Stitch!
                wl, flux, error, bpmap = stitch_arms.stitch_XSH_spectra(waves, flux, error, bpmap, scale = False)




                np.savetxt("../stitched/"+outname+"_stitched_spectrum.dat", zip(wl, flux, error, bpmap), fmt = ['%10.6e', '%10.6e', '%10.6e', '%10.6f'], header=" wl flux error bpmap")

                hbin = 10
                wl_bin, flux_bin, error_bin, bp_bin = util.bin_spectrum(wl, flux, error, bpmap.astype("bool"), hbin, weight=True)
                np.savetxt("../stitched/"+outname+"_stitched_spectrum_bin"+str(hbin)+".dat", zip(wl_bin, flux_bin, error_bin, bp_bin), fmt = ['%10.6e', '%10.6e', '%10.6e', '%10.6f'], header=" wl flux error bpmap")
                pl.errorbar(wl_bin[::1], flux_bin[::1], yerr=error_bin[::1], fmt=".k", capsize=0, elinewidth=0.5, ms=3, alpha=0.3, rasterized=True)
                pl.plot(wl_bin, flux_bin, linestyle="steps-mid", lw=1.0, alpha=0.7, rasterized=True)
                pl.plot(wl_bin, error_bin, linestyle="steps-mid", lw=1.0, alpha=0.5, color = "grey", rasterized=True)
                pl.axhline(0, linestyle="dashed", color = "black", lw = 0.4)
                pl.axhline(1, linestyle="dashed", color = "red", lw = 0.4)

                # mask = (wl > 6000) & (wl < 10000)
                # epsilon, diff, msk = np.median(error[~np.isnan(flux)][mask]), 2*np.median(error[~np.isnan(flux)][mask]), np.ones_like(error[~np.isnan(flux)][mask]).astype("bool")
                # ii = 0
                # while diff > epsilon and ii < 5:
                #     popt, pcov = optimize.curve_fit(pow, wl[~np.isnan(flux)][mask][msk], flux[~np.isnan(flux)][mask][msk], sigma=error[~np.isnan(flux)][mask][msk], p0 = [5e-11, -1.6], maxfev=5000)
                #     msk = abs(pow(wl[~np.isnan(flux)][mask], *popt) - flux[~np.isnan(flux)][mask]) < (10 - ii)*error[~np.isnan(flux)][mask]
                #     diff = np.median(abs(pow(wl[~np.isnan(flux)][mask], *popt) - flux[~np.isnan(flux)][mask]))
                #     ii += 1
                # # popt, pcov = optimize.curve_fit(pow, wl[~np.isnan(flux)][mask], flux[~np.isnan(flux)][mask], sigma=error[~np.isnan(flux)][mask], maxfev=5000, p0 = [5e-11, -1.6])
                # print(popt, np.sqrt(np.diag(pcov)))
                # pow_fit = popt[0] * wl_bin ** (popt[1])
                # pl.plot(wl_bin, pow_fit)
                # pl.plot(wl_bin, flux_bin/(pow_fit))
                # pl.errorbar(wl_bin[::1], flux_bin[::1]/(pow_fit), yerr=error_bin[::1]/(pow_fit), fmt=".k", capsize=0, elinewidth=0.5, ms=3, alpha=0.3)

                # scale = np.median(flux_bin[~np.isnan(flux_bin)])
                pl.xlim(3200, 20000)
                pl.ylim(-0.1, 1.5)
                pl.xlabel(r"Wavelength / [$\mathrm{\AA}$]")
                pl.ylabel(r'Flux density [erg s$^{-1}$ cm$^{-1}$ $\AA^{-1}$]')
                pl.savefig("../stitched/"+outname+"_stitched_spectrum_bin"+str(hbin)+".pdf")
                pl.clf()
                # pl.show()

                # exit()


if __name__ == '__main__':
    main()