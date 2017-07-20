#!/usr/bin/env python

import matplotlib
matplotlib.use('agg')

# this example comes from 
# http://yt-project.org/doc/visualizing/volume_rendering.html

import math
import sys
import pylab
import numpy as np

import yt
import yt.visualization.volume_rendering.api as vr

def doit(plotfile, fname):

    ds = yt.load(plotfile)

    if fname == "vz":
        field = ('gas', 'velocity_z')
        use_log = False
        
        vals = [-1.e7, -5.e6, 5.e6, 1.e7]
        sigma = 5.e5

        fmt = None

    elif fname == "magvel":
        field = ('gas', 'velocity_magnitude')
        use_log = True
        
        vals = [1.e5, 3.16e5, 1.e6, 3.16e6, 1.e7]
        sigma = 0.1

    elif fname == "enucdot":
        field = ('boxlib', 'enucdot')
        use_log = True

        vals = [1.e16, 3.162e16, 1.e17, 3.162e17, 1.e18]
        vals = list(10.0**numpy.array([16.5, 17.0, 17.5, 18.0, 18.5]))
        sigma = 0.05

        fmt = "%.3g"
    

    dd = ds.all_data()

    mi = min(vals)
    ma = max(vals)
    
    if use_log:
        mi, ma = np.log10(mi), np.log10(ma)
        
    # this is hackish, but there seems to be no better way to set whether
    # you are rendering logs
    ds._get_field_info("density").take_log = use_log

    print mi, ma

    # Instantiate the ColorTransferfunction.
    tf =  vr.ColorTransferFunction((mi, ma))
    cm = "gist_rainbow"

    for v in vals:
        if (use_log):
            print v, math.log10(v)
            tf.sample_colormap(math.log10(v), sigma**2, colormap=cm) #, alpha=0.2)
        else:
            tf.sample_colormap(v, sigma**2, colormap=cm) #, alpha=0.2)


    # an attempt to get around the periodic assumption -- suggested by Nathan
    #root_dds = pf.domain_width/pf.domain_dimensions
    #half_w = pf.domain_width/2.# - root_dds
    #half_w[2] -= root_dds[2]
    #reg = pf.region(pf.domain_center, pf.domain_center-half_w, pf.domain_center+half_w)


    # alternate attempt
    ds.periodicity = (True, True, True)

    # Create a camera object
    # Set up the camera parameters: center, looking direction, width, resolution
    c = (ds.domain_right_edge + ds.domain_left_edge)/2.0
    L = np.array([0.0, -1.0, -1.0])

    north_vector=[0.0,0.0,1.0]



    im, sc = yt.volume_render(ds, [field])

    source = sc.get_source(0)
    source.set_transfer_function(tf)
    source.transfer_function.grey_opacity=True

    # note: this needs to come AFTER the set_transfer_function(), or
    # else we get an AttributeError
    #sc.annotate_domain(ds)

    cam = sc.camera

    cam.resolution = (720,720)
    cam.set_width(3.0*ds.domain_width)
    cam.focus = c

    # we cannot seem to switch views and get an image...
    #cam.switch_view(north_vector=north_vector)
    #cam.switch_view(normal_vector=L)
    #cam.switch_orientation(L, north_vector)

    # make an image
    im = sc.render("test.png")


if __name__ == "__main__":

    # Choose a field
    fname = "vz"
    plotfile = ""


    try: plotfile = sys.argv[1]
    except: sys.exit("ERROR: no plotfile specified")

    try: fname = sys.argv[2]
    except: pass

    doit(plotfile, fname)


        
