#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

import sys
sys.path.append("/Users/jselsing/github/XSGRB_reduction_scrips/py/")
import stitch_arms, util

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as pl
import seaborn as sns; sns.set_style('ticks')
import glob
import itertools
from scipy.signal import medfilt
from scipy import interpolate
from NonnegMFPy import nmf

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

# Execution of scipt
if __name__ == '__main__':

    # Get stitched spectra
    objects = glob.glob("../stitched/*spectrum.dat")
    # Burst names
    names = list(set([ii.split("/")[-1] for ii in glob.glob("../data/*")]))
    # Define gold_sample of high S/N bursts - commented for speed
    # gold_sample = []
    # for ii, kk in enumerate(objects):
    #     obj = np.genfromtxt(kk)
    #     mask = ~obj[:, 3].astype(bool)
    #     SN = np.median(obj[:, 1][mask]/obj[:, 2][mask])
    #     if SN > 5:
    #         print(kk)
    #         gold_sample.append(kk.split("/")[-1][:10])

    # gold_sample = ['GRB090926', 'GRB091018', 'GRB100316B', 'GRB100316D', 'GRB100418A', 'GRB100814A','GRB110715A', 'GRB111005A', 'GRB111209A', 'GRB111211A', 'GRB111228A', 'GRB120119A', 'GRB120327A', 'GRB120815A', 'GRB120909A', 'GRB121024A', 'GRB130408A', 'GRB130418A', 'GRB130427A', 'GRB131030A', 'GRB131231A', 'GRB140213A', 'GRB141028A', 'GRB141109A', 'GRB150403A', 'GRB150514A', 'GRB150821A', 'GRB151021A', 'GRB151027B', 'GRB160203A', 'GRB160410A', 'GRB160625B', 'GRB161023A']

    # Manually define gold sample
    gold_sample = ["GRB090926", "GRB091018", "GRB100316B", "GRB100418A", "GRB100814A", "GRB111008A", "GRB111209A", "GRB111211A", "GRB111228A", "GRB120119A", "GRB120327A", "GRB120712A", "GRB120815A", "GRB120909A", "GRB121024A", "GRB130408A", "GRB130418A", "GRB130427A", "GRB130606A0", "GRB131030A", "GRB131231A", "GRB140213A", "GRB140506A", "GRB141028A", "GRB141109A", "GRB150403A", "GRB150514A", "GRB150821A", "GRB150915A", "GRB151021A", "GRB151027B", "GRB151031A", "GRB160117B", "GRB160203A", "GRB160410A", "GRB160625B", "GRB161023A"]

    # Define high hydrogen column to look for molecules.
    # gold_sample = ["GRB121027A", "GRB130612A", "GRB111008A", "GRB160203A", "GRB151021A", "GRB141109A", "GRB130427B", "GRB120815A", "GRB120716A", "GRB120327A"]

    # Define with detection of molecules.
    # gold_sample = ["GRB140506A", "GRB120119A", "GRB120327A", "GRB120815A", "GRB121024A"]




    # Select gold bursts in stitched spectra
    targets = [ii for ii in gold_sample if ii in names]
    # targets = [ii for ii in names]
    # Number of targets
    nobj = len(targets)
    print("Bulding composite out of: "+str(nobj)+" files")
    # Read in burst list for redshifts and intervening systems
    burst_list = np.genfromtxt("/Users/jselsing/github/XSGRB-sample-paper/data/burst_list.dat", dtype=None)
    names = [ii[0] for ii in burst_list]
    redshifts = [float(ii[1]) for ii in burst_list]
    intervening = [ii[4].split(",") for ii in burst_list]
    # Remove empty elements
    for ii, kk in enumerate(intervening):
        if kk[0] == "nan":
            intervening[ii] = []
        else:
            intervening[ii] = [float(ll) for ll in kk]

    # Find redshifts and intervening lists for selected sample
    gold_z, gold_intervening, gold_names  = [0]*len(targets), [0]*len(targets), [0]*len(targets)
    for ii, kk in enumerate(targets):
        idx = list(itertools.chain.from_iterable(np.where(kk == np.array(names))))
        if len(idx) != 1:
            idx = [idx[0]]
        gold_names[ii] =  names[idx[0]]
        gold_z[ii] = redshifts[idx[0]]
        gold_intervening[ii] = intervening[idx[0]]

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
    for ii, kk in enumerate(gold_names):
        for jj, pp in enumerate(gold_intervening[ii]):
            if pp:
                bpmap[ii] = mask_intervening_systems(wl[ii], bpmap[ii], pp, mask_width = 100.)

    # Iterate through wavelength array to find wavlength-range
    wl_start, wl_end = [], []
    for ii, kk in enumerate(gold_names):
        wave = wl[ii] / (1 + gold_z[ii])
        wl_start.append(min(wave))
        wl_end.append(max(wave))
    wl_new = np.arange(min(wl_start), max(wl_end), 0.1)

    # Interpolate spectra to common sampling at rest
    for ii, kk in enumerate(gold_names):
        wave = wl[ii] / (1 + gold_z[ii])
        f_spline = interpolate.interp1d(wave, f[ii], bounds_error=False, fill_value=0)
        e_spline = interpolate.interp1d(wave, e[ii], bounds_error=False, fill_value=1e2)
        bpmap_spline = interpolate.interp1d(wave, bpmap[ii], bounds_error=False, fill_value=666)

        f[ii] = f_spline(wl_new)
        e[ii] = e_spline(wl_new)
        bpmap[ii] = bpmap_spline(wl_new)
        bpmap[ii][f[ii] > 3] = 10
        bpmap[ii][f[ii] < -0.2] = 10

    f, e, bpmap = np.array(f), np.array(e), np.array(bpmap)


    # Mask regions in individual objects
    for ii, kk in enumerate(gold_names):
        if 'GRB091018A' in kk:
            bpmap[ii][(wl_new > 15000)] = 1e2
        elif 'GRB100316B' in kk:
            bpmap[ii][(wl_new > 9500)] = 1e2
        elif 'GRB100418A' in kk:
            bpmap[ii][(wl_new > 15000)] = 1e2
        elif 'GRB100814A' in kk:
            bpmap[ii][(wl_new > 9300)] = 1e2
        elif 'GRB111008A' in kk:
            bpmap[ii][(wl_new > 3500)] = 1e2
        elif 'GRB111211A' in kk:
            bpmap[ii][(wl_new > 16000)] = 1e2
        elif 'GRB111228A' in kk:
            bpmap[ii][(wl_new > 10500)] = 1e2
        elif 'GRB120327A' in kk:
            bpmap[ii][(wl_new > 6000)] = 1e2
        # elif 'GRB121024A' in kk:
            # bpmap[ii][(wl_new > 5250) & (wl_new < 5500)] = 1e2
        # elif 'GRB130408A' in kk:
        #     bpmap[ii][(wl_new > 400) & (wl_new < 1178)] = 1e2
        # elif 'GRB130418A' in kk:
        #     bpmap[ii][(wl_new > 1350) & (wl_new < 1400)] = 1e2
        # elif 'GRB131030A' in kk:
        #     bpmap[ii][(wl_new > 1300) & (wl_new < 1400)] = 1e2
        # elif 'GRB131231A' in kk:
        #     bpmap[ii][(wl_new > 1800) & (wl_new < 2000)] = 1e2
        # elif 'GRB140213A' in kk:
        #     bpmap[ii][(wl_new > 1300) & (wl_new < 1500)] = 1e2
        elif 'GRB151021A' in kk:
            bpmap[ii][(wl_new > 6800)] = 1e2
        elif 'GRB1602203A' in kk:
            bpmap[ii][(wl_new > 5000)] = 1e2

    # for ii, kk in enumerate(f):
    #     mask = ~bpmap[ii].astype("bool")
    #     pl.plot(wl_new[mask], kk[mask])
        # pl.plot(wl_new, bpmap[ii])
        # pl.title(targets[ii])
    # pl.ylim(-0.1, 1.3)
    # pl.xlim(2700, 2900)
    # pl.show()

    g = nmf.NMF(f, V=1/e**2., M=~bpmap.astype("bool"), n_components=5)
    chi2, time_used = g.SolveNMF() # Yes, it's that simple!
    print(g.W)
    exit()


    ma_arr = np.ma.array(f, mask=bpmap.astype("bool"))
    median = np.ma.median(ma_arr, axis=0)
    n = np.sum((~bpmap.astype("bool")).astype("int"), axis=0)
    pl.plot(wl_new, medfilt(n, 5))
    pl.xlabel(r'Rest wavelength [$\mathrm{\AA}$]', fontsize = 25)
    pl.ylabel(r'Normalised flux density', fontsize = 25)
    pl.savefig("../figures/n_bursts.pdf")
    pl.clf()
    ma_e_arr = np.ma.array(e, mask=bpmap.astype("bool"))
    median_err = np.sqrt((np.pi/2) * np.ma.sum(ma_e_arr**2, axis=0)/n**2)

    pl.plot(wl_new, median)
    pl.plot(wl_new, median_err)
    pl.show()

    std = medfilt(np.ma.std(ma_arr, axis=0), 11)
    wmean, wmeanerr, _ = util.avg(f, e, mask=bpmap.astype("bool"), axis=0, weight=True)


    # Saving to .dat file
    dt = [("wl", np.float64), ("median", np.float64), ("median_err", np.float64), ("wmean", np.float64), ("wmeanerr", np.float64)]
    data = np.array(zip(wl_new, median, median_err, wmean, wmeanerr), dtype=dt)
    np.savetxt("../gold_comp.dat", data, header="wavelength   median median_error   weighted_mean weighted_mean_err", fmt = ['%1.5e', '%1.5e', '%1.5e', '%1.5e', '%1.5e'], delimiter="\t")


    dib_spec_txt = np.genfromtxt("static/dibs_spec_for_dovi.txt")
    dibspec_wl, dibspec = dib_spec_txt[:,0], dib_spec_txt[:,1]


    molec = np.genfromtxt("/Users/jselsing/Work/spectral_libraries/H2_Abs/Synspec/Synspec_H2_20_z2.35.txt")
    wl = molec[:, 0]/3.3538
    mol = molec[:, 1]

    import matplotlib as mpl
    label_size = 20
    mpl.rcParams['xtick.labelsize'] = label_size
    mpl.rcParams['ytick.labelsize'] = label_size
    fig, ax = pl.subplots(figsize=(16, 6))
    skip = 5
    # ax.errorbar(wl_new[::skip], wmean[::skip], yerr=wmeanerr[::skip], fmt=".k", capsize=0, elinewidth=0.5, ms=3)
    ax.plot(wl_new, medfilt(med, 1), lw = 1, linestyle="steps-mid", alpha=0.7, color = "black", rasterized=True)
    ax.plot(dibspec_wl, dibspec)
    ax.plot(wl, mol)
    # ax.plot(wl_new, wmean/wmeanerr, lw = 0.5, linestyle="steps-mid", alpha=0.5)
    # ax.fill_between(wl_new, wmean - wmeanstd, wmean + wmeanstd, alpha = 0.5, lw = 0.0)
    ax.fill_between(wl_new, np.zeros_like(wl_new), std, alpha = 0.7, lw = 0.0, rasterized=True)
    ax.set_xlabel(r'Rest wavelength [$\mathrm{\AA}$]', fontsize = 25)
    ax.set_ylabel(r'Normalised flux density', fontsize = 25)
    # pl.tight_layout()
    ax.set_ylim(0, 1.5)
    ax.set_xlim(4490, 5390)
    # pl.tight_layout()
    fig.savefig('../figures/XSHcomposite.pdf')
    # pl.show()
    # exit()
    # mask = wl_new > 1300 & wl_new < 1700
    ax.set_xlim(1000, 5000)
    ax.set_ylim(-0.1, 1.3)

    # #Overplot lines
    import lineid_plot
    fit_line_positions = np.genfromtxt('static/grblines.dat', dtype=None)

    # pl.tight_layout()
    linelist = []
    linenames = []
    for n in fit_line_positions:
        linelist.append(float(n[0]))
        linenames.append(str(n[1]))
    # lineid_plot.plot_line_ids(wl_new, wmean, linelist, linenames, ax = ax)
    pl.gcf().subplots_adjust(bottom=0.15)
    # fig.savefig('../figures/XSHcompositezoom.pdf')
    # pl.show()
    ax.set_xlim(1020, 1040)
    lineid_plot.plot_line_ids(wl_new, wmean, linelist, linenames, ax = ax)
    # pl.show()
    fig.savefig('../figures/XSHcompositezoomzoom.pdf')
    # pl.show()
