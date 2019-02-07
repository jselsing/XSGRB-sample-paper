#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

import datetime


def timediff(t1, t2):

    time1 = datetime.datetime(int(t1[0:4]), int(t1[5:7]), int(t1[8:10]), hour=int(round(float(t1[11:13]))), minute=int(round(float(t1[14:16]))), second=int(round(float(t1[17:19]))))
    time2 = datetime.datetime(int(t2[0:4]), int(t2[5:7]), int(t2[8:10]), hour=int(round(float(t2[11:13]))), minute=int(round(float(t2[14:16]))), second=int(round(float(t2[17:19]))))

    dt = time2 - time1
    return dt.total_seconds()/(60)

if __name__ == '__main__':

    t_trigger = "2019:01:14T20:57:03"
    t_obs = "2019:01;18T00:59:27.051"

    print(timediff(t_trigger, t_obs))
