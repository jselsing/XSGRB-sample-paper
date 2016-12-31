#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

import glob
from astropy.io import fits
from scipy.special import wofz, erf
from scipy.optimize import curve_fit
import seaborn; seaborn.set_style('ticks')
import numpy as np
import matplotlib.pyplot as pl
from scipy import interpolate


def voigt_base(x, amp=1, cen=0, sigma=1, gamma=0):
    """1 dimensional voigt function.
    see http://en.wikipedia.org/wiki/Voigt_profile
    """
    z = (x-cen + 1j*gamma)/ (sigma*np.sqrt(2.0))
    return amp * wofz(z).real / (sigma*np.sqrt(2*np.pi))


def multi_voigt(x, *params):
    # Multiple voigt-profiles for telluric resolution estimate
    sigma = params[0]
    gamma = params[1]
    c = params[2]
    a = params[3]



    multivoigt = 0
    for ii in range(4, int(4 + len(params[4:])/2)):
        if sigma < 0 or gamma < 0 or params[ii] > 0:
            a = 1e10
        multivoigt += voigt_base(x, params[ii], params[int(len(params[4:])/2) + ii], sigma, gamma)
    return multivoigt + c + a * x

# ------------------------------------------------------------------------------
# Copied these two function from ppxf_util.py by Michelle Cappellari
def vac_to_air(lam_vac):
    """
    Convert vacuum to air wavelengths using
    equation (1) of Ciddor 1996, Applied Optics 35, 1566
        http://dx.doi.org/10.1364/AO.35.001566

    :param lam_vac - Wavelength in Angstroms
    :return: lam_air - Wavelength in Angstroms

    """
    sigma2 = (1e4/lam_vac)**2
    fact = 1 + 5.792105e-2/(238.0185 - sigma2) + 1.67917e-3/(57.362 - sigma2)

    return lam_vac/fact

# ------------------------------------------------------------------------------

def air_to_vac(lam_air):
    """
    Convert air to vacuum wavelengths using
    equation (1) of Ciddor 1996, Applied Optics 35, 1566
        http://dx.doi.org/10.1364/AO.35.001566
    :param lam_air - Wavelength in Angstroms
    :return: lam_vac - Wavelength in Angstroms

    """
    sigma2 = (1e4/lam_air)**2
    fact = 1 + 5.792105e-2/(238.0185 - sigma2) + 1.67917e-3/(57.362 - sigma2)

    return lam_air*fact

#------------------------------------------------------------------------------



def main():
    # Get list of burst names
    bursts = list(set([ii.split("/")[-1].split("_")[0] for ii in glob.glob("../specdb/*")]))
    # Loop through bursts
    for ii in bursts:
        # if not "GRB161023A" in ii:
        #     continue
        files = glob.glob("../specdb/"+ii+"*")
        for ll in files:

            file = fits.open(ll)
            wl, flux = file[1].data.field("WAVE") / (1 + file[0].header['HIERARCH ESO QC VRAD BARYCOR']/3e5), file[1].data.field("FLUX")
            arm = ll[-8:-5]
            # pl.plot(wl, flux)
            # pl.show()
            if arm == "UVB":
                continue

            elif arm == "VIS":
                mask = (wl > 8996) & (wl < 9043)
                cens = [9003, 9006, 9009, 9014, 9019, 9025, 9028, 9034]
                amps = len(cens) * [min(flux[mask]) - max(flux[mask])]
            elif arm == "NIR":
                mask = (wl > 17580) & (wl < 17644)
                cens = [1.76100617e+04, 1.76245674e+04, 1.76304265e+04]
                amps = len(cens) * [-1e-17]

            min_idx, max_idx = min(*np.where(mask)), max(*np.where(mask))

            file2d = [kk for kk in glob.glob("../data/"+ii+"/*") if arm in kk and "2D" in kk][0]
            tell_file2D = fits.open(file2d)
            print(ll, file2d)
            v_len = np.shape(tell_file2D[0].data)[0]

            profile = np.median(tell_file2D[0].data[int(v_len/4):int(-v_len/4), min_idx:max_idx], axis= 1)
            xarr = np.arange(len(profile))



            p0 = [max(profile), len(xarr)/2, 5, 0]
            try:
                popt, pcuv = curve_fit(voigt_base, xarr, profile, p0=p0)
                fwhm_g, fwhm_l = 2.35 * popt[2], 2*popt[3]
                fwhm_g_var, fwhm_l_var = 2.35 * pcuv[2, 2], 2*pcuv[3, 3]
                fwhm = 0.5346 * fwhm_l + np.sqrt(0.2166 * fwhm_l**2 + fwhm_g**2)
                dfdl = 0.5346 - 0.5 * ((0.2166 * fwhm_l**2 + fwhm_g**2) ** (-3/2)) *(2 * 0.2166 * fwhm_l)
                dfdg = - 0.5 * ((0.2166 * fwhm_l**2 + fwhm_g**2) ** (-3/2)) *(2 * fwhm_g)
                fwhm_err = np.sqrt((dfdl**2) * (fwhm_g_var) + (dfdg**2) * (fwhm_l_var))
                seeing_fwhm = fwhm*tell_file2D[0].header["CD2_2"]
                seeing_fwhm_err = fwhm_err*tell_file2D[0].header["CD2_2"]
                p0 =  [0.2, 0.01, max(flux[mask]), 0] + amps + cens
                popt, pcuv = curve_fit(multi_voigt, wl[mask], flux[mask], p0=p0)
                midwl = np.median(wl[mask])
                R = midwl / (popt[0]*2.35)
                Rerr =   R - midwl /((popt[0] + np.sqrt(np.diag(pcuv)[0]))*2.35)
                x = np.arange(min(wl[mask]), max(wl[mask]), 0.01)
                pl.plot(x, multi_voigt(x, *popt), label="R = "+str(int(np.around(R, decimals = -2))) + " +- " + str(int(np.around(Rerr, decimals = -2))))
                pl.plot(wl[mask], flux[mask], label="Seeing FWHM = "+str(np.around(seeing_fwhm, decimals = 2)) + " +- " + str(np.around(seeing_fwhm_err, decimals = 2)))
                pl.xlabel(r"Wavelength / [$\mathrm{\AA}$]")
                pl.ylabel(r'Flux density [erg s$^{-1}$ cm$^{-1}$ $\AA^{-1}$]')
                pl.ylim((min(flux[mask]), max(flux[mask])*1.10))
                pl.legend()
                pl.savefig("../figures/res/"+ii+arm+"_resolution.pdf")
                # pl.show()
                pl.clf()
            except:
                continue




if __name__ == '__main__':
    main()