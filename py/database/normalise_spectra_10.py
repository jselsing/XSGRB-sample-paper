#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

import sys

import stitch_arms, util

from glob import glob
import numpy as np
import matplotlib.pyplot as pl
from scipy import optimize
import os
import xsh_norm.automatic
from scipy.signal import medfilt

from numpy.polynomial import chebyshev
from scipy.interpolate import interp1d
from scipy.interpolate import splrep, splev


def main():
    """
    Script to loop though bursts and normalise arms
    """
    data_dir = "../data/"
    bursts = glob(data_dir+"*")

    # Loop through burst names
    for ii, kk in enumerate(bursts):
        # List of extractions in burst dir
        extractions = [ll for ll in glob(kk+"/*") if ".dat" in ll or ".spec" in ll]
        # print(extractions)

        # Loop through extractions and normalize
        for ll in extractions:

            if "GRB10" in ll:
                print("Normalizing file:")
                print(ll)
                burst_name = ll.split("/")[2]
                if "spec" in ll:
                    outname = (ll.split("/")[-1]).replace("skysub", "")[:-5]
                else:
                    outname = (ll.split("/")[-1]).replace("skysub", "")[:-4]
                print(kk+"/"+outname+"_norm.npy")

                if not os.path.exists(kk+"/"+outname+"_norm.npy"):
                    wl, flux, error, bpmap = stitch_arms.load_array(np.genfromtxt(ll))
                    mask = ~bpmap.astype("bool")
                    normalise = xsh_norm.automatic.xsh_norm(wl[mask], flux[mask], error[mask], bpmap, wl, flux, error, kk+"/"+outname)

                    # if burst_name == "GRB090809A":
                    #     normalise.leg_order = 4
                    #     normalise.endpoint_order = 4
                    # elif burst_name == "GRB091018A":
                    #     normalise.leg_order = 4
                    #     normalise.lover_mask = -5e-17
                    #     normalise.exclude_width = 2
                    #     normalise.sigma_mask = 10

                    # normalise.leg_order = 5 #float(raw_input('Filtering Chebyshev Order (default = 3) = '))
                    # normalise.endpoint_order = 5#float(raw_input('Order of polynomial used for placing endpoints (default = 3) = '))
                    normalise.exclude_width = 2#float(raw_input('Exclude_width (default = 5) = '))
                    normalise.sigma_mask = 3#float(raw_input('Sigma_mask (default = 5) = '))
                    normalise.lover_mask = -1e-16#float(raw_input('Lover bound mask (default = -1e-17) = '))
                    # normalise.tolerance = 0.5#float(raw_input('Filtering tolerance (default = 0.25) = '))
                    # normalise.leg_order = #float(raw_input('Filtering Chebyshev Order (default = 3) = '))
                    normalise.spacing = 400#float(raw_input('Spacing in pixels between spline points (default = 150) = '))
                    normalise.division = 5#float(raw_input('Maximum number of allowed points (default = 300) = '))
                    normalise.endpoint_order = 2#float(raw_input('Order of polynomial used for placing endpoints (default = 3) = '))
                    normalise.endpoint = "t"#str(raw_input('Insert endpoint before interpolation(y/n)? '))
                    try:
                        normalise.run()
                    except:
                        pass


if __name__ == '__main__':
    main()


