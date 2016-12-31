#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

import numpy as np
import glob

def main():
    files = [ii for ii in glob.glob("../data/*/*") if "spec" in ii]

    for ii in files:
        file = np.genfromtxt(ii)
        data = zip(file[:, 0], file[:, 1]*1e-17, file[:, 2]*1e-17)
        np.savetxt(ii, data)
        # print(file[:, 0])
    # print(files)
    # pass

if __name__ == '__main__':
    main()