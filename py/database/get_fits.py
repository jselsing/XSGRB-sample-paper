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
    data_dir = "../specdb_input/"

    # Loop through data directories and copy files to local directories
    for ii, kk in zip(burst_names, dirs):
        # if not "GRB160203A" in kk:
        #     continue
        # Create directory to contain data
        if not os.path.exists(data_dir+ii):
            os.makedirs(data_dir+ii)
        # Files to be moved to local dir
        print(kk)
        filelist = [ll for ll in glob.glob(kk+"/reduced_data/*/*/*/*") if ("too" in ll.lower() or "rrm" in ll.lower()) and "IDP" in ll]
        OBs = list(set([ll.split("/")[-4] for ll in filelist]))
        lst = [ll for ll in glob.glob(kk+"/*") if "combine" in ll]
        if len(lst) != 0:
            OBs.append("com")

        arms = list(set([ll.split("/")[-3] for ll in filelist]))

        copy_filelist, filenames, count = [0]*len(OBs)*len(arms), [0]*len(OBs)*len(arms), 0

        for dd in OBs:
            for ss in arms:
                if dd == "com":
                    copy_filelist[count] = [ll for ll in glob.glob(kk+"/reduced_data/OB1/"+ss+"/*/*") if ("too" in ll.lower() or "rrm" in ll.lower()) and "IDP" in ll][0]
                else:
                    copy_filelist[count] = [ll for ll in glob.glob(kk+"/reduced_data/"+dd+"/"+ss+"/*/*") if ("too" in ll.lower() or "rrm" in ll.lower()) and "IDP" in ll][0]
                filenames[count] = dd+ss+"_IDP.fits"
                count += 1

        # print(copy_filelist, filenames)
        # Move files
        for tt, pp in zip(copy_filelist, filenames):
            if not os.path.exists(data_dir+ii+"/"+pp):
                print("Moving file: "+tt)
                print("To file: "+data_dir+ii+"/"+pp)
                copyfile(tt, data_dir+ii+"/"+pp)


if __name__ == '__main__':
    main()