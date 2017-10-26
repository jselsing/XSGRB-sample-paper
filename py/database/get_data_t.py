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

        # Create directory to contain data
        if not os.path.exists(data_dir+ii):
            os.makedirs(data_dir+ii)
        # Files to be moved to local dir
        print(kk)
        spec_filelist = [ll for ll in glob.glob(kk+"/*/*") if ii[:-3] in ll.split("/")[-1] and ll[-4:] == "spec"]
        d2_filelist = [ll for ll in glob.glob(kk+"/*/*") if ii[:-3] in ll.split("/")[-1] and ll[-4:] == "fits"]
        # print(spec_filelist)
        OBs = ["OB1"]

        arms = list(set([ll.split("/")[-1].split(".")[0].upper()[-3:] for ll in d2_filelist]))

        copy_filelist, filenames, count = 2*[0]*len(OBs)*len(arms), 2*[0]*len(OBs)*len(arms), 0

        for dd in OBs:
            for ss in arms:
                copy_filelist[count] = [ll for ll in glob.glob(kk+"/*/*") if ss.lower() in ll and ll[-4:] == "fits"][0]
                filenames[count] = dd+ss+"_2D.fits"
                count += 1

        for dd in OBs:
            for ss in arms:
                try:
                    copy_filelist[count] = [ll for ll in glob.glob(kk+"/*/*") if ss.lower() in ll and ll[-4:] == "spec"][0]
                    filenames[count] = dd+ss+"_1D.spec"
                    count += 1
                except:
                    print("No extraction for:" +ss +" in "+ kk)
                    copy_filelist[count] = 0
                    filenames[count] = 0
                    count += 1

        # Move files
        for tt, pp in zip(copy_filelist, filenames):
            if tt != 0:
                if not os.path.exists(data_dir+ii+"/"+pp):
                    print("Moving file: "+tt)
                    print("To file: "+data_dir+ii+"/"+pp)
                    copyfile(tt, data_dir+ii+"/"+pp)

if __name__ == '__main__':
    main()