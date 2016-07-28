#!/usr/bin/env python

import yt

ds = yt.load("smallplt16623")

p = yt.ProjectionPlot(ds, "z", "density")
p.save()

