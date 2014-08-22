#!/bin/env python

import os
import string

import slice3plot

runs = ["base",
        "CG",
        "gammae",
        "no_trace_grav",
        "ppm2",
        "temp_fix_1",
        "temp_fix_2",
        "trans_eos",
        "trans_rhoe"]

prefix = "plt"

var = "temperature"

for r in runs:

    # find the last plotfile output to that directory
    plot_num = -1
    for file in os.listdir(r):
        if os.path.isdir(r+"/"+file) and file.startswith(prefix):
            index = string.rfind(file, prefix)
            plot_num = max(int(file[index+len(prefix):]), plot_num)

    file = r + "/" + prefix + "{:05}".format(plot_num)

    slice3plot.doit(file, var, 0, r, r+".png")

