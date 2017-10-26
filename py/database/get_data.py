#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

import os
import glob
from shutil import copyfile
import time

def main():
    """
    Small script to crawl through ERDA can collect burst extractions
    """
    dirs = [ii for ii in glob.glob("/Volumes/io.erda.dk/XSGRB/*") if "GRB" in ii.split("/")[-1]]
    burst_names = [ii.split("/")[-1] for ii in dirs]
    data_dir = "../data/"

    # Loop through data directories and copy files to local directories
    for ii, kk in zip(burst_names, dirs):
        if not "GRB120119A" in kk:
            continue
        # Create directory to contain data
        if not os.path.exists(data_dir+ii):
            os.makedirs(data_dir+ii)
        # Files to be moved to local dir
        print(kk)
        spec_filelist = [ll for ll in glob.glob(kk+"/*") if "ext.dat" in ll]
        d2_filelist = [ll for ll in glob.glob(kk+"/*") if "skysub.fits" in ll]

        OBs = list(set([ll.split("/")[-1][3:6] for ll in d2_filelist]))

        lst = [ll for ll in glob.glob(kk+"/*") if "combine" in ll]
        if len(lst) != 0:
            OBs.append("com")

        arms = list(set([ll.split("/")[-1][:3] for ll in d2_filelist]))
        copy_filelist, filenames, count = 2*[0]*len(OBs)*len(arms), 2*[0]*len(OBs)*len(arms), 0

        for dd in OBs:
            for ss in arms:
                try:
                    copy_filelist[count] = [ll for ll in glob.glob(kk+"/*") if ss in ll and dd in ll and ("skysub.fits" in ll or "combined.fits" in ll)][0]
                    filenames[count] = dd+ss+"_2D.fits"
                    count += 1
                except:
                    print("No 2D image for:" +ss +" in "+ kk)
                    copy_filelist[count] = 0
                    filenames[count] = 0
                    count += 1

        for dd in OBs:
            for ss in arms:
                copy_filelist[count] = [ll for ll in glob.glob(kk+"/*") if ss in ll and dd in ll and "ext.dat" in ll]
                if len(copy_filelist[count]) != 1:
                    for ww in copy_filelist[count]:
                        if "optext" in ww:
                            copy_filelist[count] = ww
                elif len(copy_filelist[count]) == 1:
                    for ww in copy_filelist[count]:
                        copy_filelist[count] = ww
                filenames[count] = dd+ss+"_1D.dat"
                count += 1



        # Move files
        for tt, pp in zip(copy_filelist, filenames):
            print(tt, pp)
            if tt != 0:
                if not os.path.exists(data_dir+ii+"/"+pp):
                    try:
                        print("Moving file: "+tt)
                        print("To file: "+data_dir+ii+"/"+pp)
                        copyfile(tt, data_dir+ii+"/"+pp)
                    except:
                        print("Skipping file:")
                        print(tt)


if __name__ == '__main__':
    main()