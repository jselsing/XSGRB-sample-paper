#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

import numpy as np
import pandas as pd


def main():
    complete_sample = pd.read_csv("Burst list - CSV_master.csv")
    # idx_observed =
    # idx_afterglow = (complete_sample.Afterglow.values == "Yes")
    # rd = complete_sample.Redshift.values[idx_afterglow]
    # print(np.sum((complete_sample.Afterglow.values == "Yes").astype("int")), len(rd[~np.isnan(rd)]))
    # exit()
    swift_table = pd.read_table("grb_table_1482495106.txt", delimiter="\t", dtype=None)
    # print(swift_table)
    idx = [kk for kk, ii in enumerate(swift_table["GRB"]) if ii in complete_sample["GRB"].values]
    # for ii in swift_table["XRT 11 Hour Flux (0.3-10 keV) [10^-11 erg/cm^2/s]"][idx][::-1]:
    #     print(ii)
    # for ii in swift_table["BAT Fluence (15-150 keV) [10^-7 erg/cm^2]"][idx][::-1]:
    #     print(ii)
    # for ii in swift_table["XRT Column Density (NH) [10^21 cm^-2]"][idx][::-1]:
    #     print(ii)
    for ii, kk in zip(swift_table["Time [UT]"][idx][::-1], swift_table["GRB"][idx][::-1]):
        print("20"+str(kk[:2])+":"+str(kk[2:4])+":"+str(kk[4:6])+"T"+str(ii))
    # print(swift_table["BAT Fluence (15-150 keV) [10^-7 erg/cm^2]"][idx])
    # burst_name, UT, Fluence, XRT, NH = table[0, :], table[1, :], table[7, :], table[18, :], table[22, :]



# "XRT 11 Hour Flux (0.3-10 keV) [10-11 erg/cm2/s]"


if __name__ == '__main__':
    main()