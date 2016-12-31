#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function


import numpy as np
from astropy.io import fits
import matplotlib.pyplot as pl
import seaborn as sns; sns.set_style('ticks')
import glob
import itertools
from scipy.signal import medfilt
from scipy import interpolate
from astropy.convolution import Gaussian1DKernel, Gaussian2DKernel, convolve


def mask_intervening_systems(wl, bpmap, z, mask_width = 5.):
    """Small script to update bad-pixel map to mask out lines in intervening systems.

    Args:
        wl (np.array): Wavelength array corresponding to bad pixel map
        bpmap (np.array): Bap pixel map of input spectrum
        z (float): Redshift of intervening system

    Returns:
        np.array: Updated bad pixel map, where pixels in the vicinity of lines in intervening systems has been replaced by the value 1216
    """
    # Read in the GRB linelist from Christensen et al. 2011
    linelist = np.genfromtxt('static/grblines_strong.dat', dtype=None)
    # Pick out all wavlengths and move to the redshift of the intervening system
    waves = [kk*(1. + z) for ii, (kk, jj) in enumerate(linelist)]

    # List comprehension for each line where pixels within 5 AA on each side of the line are replaced by ones.
    # 50 km/s mask intervening systems are
    c = 299792.458
    mask_list = [((wl < kk + mask_width*kk / (c))  & (wl > kk - mask_width*kk / (c))).astype("int") for kk in waves]

    # Sum all the individual line masks
    mask = np.sum(np.array(mask_list), axis = 0)

    # Replace all masked values by 1216
    bpmap[mask != 0] = 1216

    # Return updated bap pixel map
    return bpmap



def find_nearest(array, value):
        idx = (np.abs(array-value)).argmin()
        return idx


def main():
    """
    Script to recalibrate wavelengths by cross-correlation with preliminary composite
    """

    # Get stitched spectra
    objects = glob.glob("../stitched/*spectrum.dat")
    # Burst names
    names_inp = list(set([ii.split("/")[-1] for ii in glob.glob("../data/*")]))

    gold_names = ["GRB090926", "GRB091018", "GRB100316B", "GRB100418A", "GRB100814A", "GRB111008A", "GRB111209A", "GRB111211A", "GRB111228A", "GRB120119A", "GRB120327A", "GRB120712A", "GRB120815A", "GRB120909A", "GRB121024A", "GRB130408A", "GRB130418A", "GRB130427A", "GRB130606A0", "GRB131030A", "GRB131231A", "GRB140213A", "GRB140506A", "GRB141028A", "GRB141109A", "GRB150403A", "GRB150514A", "GRB150821A", "GRB150915A", "GRB151021A", "GRB151027B", "GRB151031A", "GRB160117B", "GRB160203A", "GRB160410A", "GRB160625B", "GRB161023A"]

    # gold_names = ["GRB161023A"]
    # Select gold bursts in stitched spectra
    targets = [ii for ii in names_inp if ii in gold_names]

    # Number of targets
    nobj = len(targets)
    # Read in burst list for redshifts and intervening systems
    burst_list = np.genfromtxt("/Users/jselsing/github/XSGRB-sample-paper/data/burst_list.dat", dtype=None)
    names_list = [ii[0] for ii in burst_list]
    redshifts_list = [float(ii[1]) for ii in burst_list]
    intervening_list = [ii[4].split(",") for ii in burst_list]
    # Remove empty elements
    for ii, kk in enumerate(intervening_list):
        if kk[0] == "nan":
            intervening_list[ii] = []
        else:
            intervening_list[ii] = [float(ll) for ll in kk]

    # Read in composite
    composite = np.genfromtxt("../gold_comp.dat")
    comp_wl, comp_f = composite[:, 0], composite[:, 1]

    linelist = np.genfromtxt("/Users/jselsing/github/H2sim/atom.dat")
    print(linelist)
    exit()

    # Find redshifts and intervening lists for selected sample
    z, names, intervening  = [0]*len(targets), [0]*len(targets), [0]*len(targets)
    for ii, kk in enumerate(targets):
        idx = list(itertools.chain.from_iterable(np.where(kk == np.array(names_list))))
        if len(idx) != 1:
            idx = [idx[0]]
        names[ii] =  names_list[idx[0]]
        z[ii] = redshifts_list[idx[0]]
        intervening[ii] = intervening_list[idx[0]]

    # Read in files and place in list
    wl, f, e, bpmap = [0]*nobj, [0]*nobj, [0]*nobj, [0]*nobj
    for ii, kk in enumerate(targets):
        files = glob.glob("../stitched/"+kk+"*spectrum.dat")
        print(kk)
        try:
            dat = np.genfromtxt(files[0])
            wl[ii] = dat[:, 0]
            tell_mask = ((wl[ii] > 6860) & (wl[ii] < 6960)) | ((wl[ii] > 7590) & (wl[ii] < 7700)) | ((wl[ii] > 13000) & (wl[ii] < 15000)) | ((wl[ii] > 18000) & (wl[ii] < 20000))
            f[ii] = dat[:, 1]
            e[ii] = dat[:, 2]
            bpmap[ii] = dat[:, 3]
            bpmap[ii][tell_mask] = 15
        except:
            print("Ooops")
            wl[ii] = np.zeros_like(wl[ii - 1])
            f[ii] = np.zeros_like(wl[ii - 1])
            e[ii] = np.zeros_like(wl[ii - 1])
            bpmap[ii] = np.ones_like(wl[ii - 1])

    # Mask lines from intervening systems
    for ii, kk in enumerate(names):
        for jj, pp in enumerate(intervening[ii]):
            if pp:
                bpmap[ii] = mask_intervening_systems(wl[ii], bpmap[ii], pp, mask_width = 100.)

    # Interpolate composite and move to burst redshift
    for ii, kk in enumerate(names):

        # wave = comp_wl * (1 + z[ii])
        f_spline = interpolate.interp1d(comp_wl, comp_f, bounds_error=False, fill_value=0)
        # composite = f_spline(wl[ii])

        # Cross correlate with redshifted spectrum and find velocity offset
        offsets = np.arange(-0.01, 0.01, 0.00001)
        correlation = [0]*len(offsets)
        m2 = (np.isnan(f[ii]) | bpmap[ii].astype("bool")) | ((wl[ii] > 13000) | (wl[ii]/z[ii] < 1000))
        for ll, pp in enumerate(offsets):
            corr_comp = f_spline(wl[ii][~m2]/(1 + z[ii] + pp))
            # m1 = np.isnan(corr_comp)
            correlation[ll] = np.median(np.correlate(corr_comp, f[ii][~m2]))



        # Index with maximum correlation
        max_idx = find_nearest(correlation, max(correlation))
        redshift_new = z[ii] + offsets[max_idx]

        print(redshift_new, kk)
        pl.plot(comp_wl*(1 + redshift_new), comp_f)
        pl.plot(wl[ii][~m2], f[ii][~m2])
        pl.ylim(-0.1, 1.5)
        pl.show()





    pass



if __name__ == '__main__':
    main()