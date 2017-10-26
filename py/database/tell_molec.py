# -*- coding: utf-8 -*-
# Adding ppxf path
import sys
sys.path.append('/Users/jselsing/github/XSspec/')
import molec
import glob
import numpy as np
from astropy.io import fits


def main():
    allfiles = glob.glob("/Users/jselsing/work/XSGRB/specdb/*.fits")
    outpath =  "../molec_tell/"


    for ii in allfiles:

        if "UVB" in ii:
            continue

        fitsfile = fits.open(ii)
        sn = np.median(fitsfile[1].data.field("FLUX"))/np.median(fitsfile[1].data.field("ERR"))
        if sn > 10:
            print(ii)
            print(sn)
            try:
                obj_name = ii.split("/")[-1].split("_")[0]
                m_fit = molec.molecFit(ii, obj_name, outpath)
                m_fit.setParams()
                m_fit.runMolec()
                m_fit.updateSpec()
            except:
                pass



if __name__ == '__main__':
    main()