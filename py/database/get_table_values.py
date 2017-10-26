# -*- coding: utf-8 -*-
# Adding ppxf path

import numpy as np
import glob
from astropy.io import fits

def main():
    all_files = glob.glob("../specdb/*")
    outlist = [0]*len(all_files)
    double_check = []
    for ii, kk in enumerate(all_files):
        # if not "GRB120119A" in kk:
        #     continue
        # print(kk)
        burst_name = kk.split("/")[-1].split("_")[0]
        OB = kk.split("/")[-1].split("_")[1][:3]

        if burst_name+OB in double_check:
            continue
        double_check.append(burst_name+OB)
        fitsfile = fits.open(kk)
        # print(fitsfile[0].header)


        airmass_start = np.around(fitsfile[0].header["HIERARCH ESO TEL AIRM START"], 1)
        airmass_end = np.around(fitsfile[0].header["HIERARCH ESO TEL AIRM END"], 1)
        airmass = str(airmass_start) + "--" + str(airmass_end)

        seeing_start = fitsfile[0].header["HIERARCH ESO TEL AMBI FWHM START"]
        if seeing_start < 0:
            seeing_start = np.nan
        seeing_end = fitsfile[0].header["HIERARCH ESO TEL AMBI FWHM END"]
        if seeing_end < 0:
            seeing_end = np.nan
        seeing =  str(np.around(np.nanmean([seeing_start, seeing_end]), 1))
        # print(seeing)

        uvb_slitwidth = fitsfile[0].header["HIERARCH ESO INS OPTI3 NAME"].replace("x11", "")
        vis_slitwidth = fitsfile[0].header["HIERARCH ESO INS OPTI4 NAME"].replace("x11", "")
        nir_slitwidth = fitsfile[0].header["HIERARCH ESO INS OPTI5 NAME"].replace("x11", "")
        slitwidth = str(uvb_slitwidth) + "/" + str(vis_slitwidth) + "/" + str(nir_slitwidth)

        date_obs = fitsfile[0].header["DATE-OBS"]



        outlist[ii] = [burst_name, date_obs, seeing, airmass, slitwidth]

    outlist = [x for x in outlist if x != 0]
    for ii in outlist:
        print(str(ii))
    # print(np.array(outlist), sep=" :-) ")


if __name__ == '__main__':
    main()