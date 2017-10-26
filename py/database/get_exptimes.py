# -*- coding: utf-8 -*-
# Adding ppxf path

import numpy as np
import glob
from astropy.io import fits


def main():

    input_file = np.genfromtxt("exptime.dat", dtype=None)


    names = np.array([ii[0].split("/")[0] for ii in input_file])
    burst_list = np.sort(list(set(names)))
    exptime = np.array([ii[1] for ii in input_file])
    arm = np.array([ii[2] for ii in input_file])
    date_obs = np.array([ii[3] for ii in input_file])

    name_check = []


    for ii in burst_list:
        idx = np.where(ii == names)

        uvb_tim, vis_tim, nir_tim = [], [], []
        for ll, kk in enumerate(names[idx]):
            if arm[idx][ll] == "UVB":
                uvb_tim.append(exptime[idx][ll])
            elif arm[idx][ll] == "VIS":
                vis_tim.append(exptime[idx][ll])
            elif arm[idx][ll] == "NIR":
                nir_tim.append(exptime[idx][ll])
        tim = str(np.sum(uvb_tim)/1000) + "/" + str(np.sum(vis_tim)/1000) + "/" + str(np.sum(nir_tim)/1000)
        # if ii == "GRB140622A":
        print(ii, tim, date_obs[idx][0])
        # exit()

    # for ii, (kk, ll, pp) in enumerate(zip(names, exptime, arm)):

        # if kk.split("/")[0] in name_check:
        # name_check.append(kk.split("/")[0])


        # exptime = ll
        # print(kk)
        # continue


if __name__ == '__main__':
    main()