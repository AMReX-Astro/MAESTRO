#!/usr/bin/env python

# this is for the wdconvect problem

import matplotlib
matplotlib.use('agg')

# this example comes from 
# http://yt-project.org/doc/visualizing/volume_rendering.html

import math
import sys
import numpy as np
import os
import pylab

import yt
#import yt.visualization.volume_rendering.api as vr
import yt.visualization.volume_rendering.old_camera as vr


def doit(plotfile):

    ds = yt.load(plotfile)

    cm = "coolwarm"

    field = ('boxlib', 'radial_velocity')
        
    vals = [-5.e6, -2.5e6, -1.25e6, 1.25e6, 2.5e6, 5.e6]
    sigma = 3.e5
        
    # select a sphere
    center = (0, 0, 0)
    R = (5.e8, 'cm')

    dd = ds.sphere(center, R)

    # Instantiate the ColorTransferfunction.
    tf =  yt.ColorTransferFunction((min(vals), max(vals)))

    # Set up the camera parameters: center, looking direction, width, resolution
    c = np.array([0.0, 0.0, 0.0])
    L = np.array([1.0, 1.0, 1.0])
    L = np.array([1.0, 1.0, 1.2])
    W = 1.5*ds.domain_width
    N = 720

    north=[0.0,0.0,1.0]

    for v in vals:
        tf.sample_colormap(v, sigma**2, colormap=cm) #, alpha=0.2)


    # alternate attempt
    ds.periodicity = (True, True, True)

    # Create a camera object
    cam = vr.Camera(c, L, W, N, transfer_function=tf, ds=ds, data_source=dd, 
                    no_ghost=False,
                    north_vector=north,
                    fields = [field], log_fields = [False])

    #cam.zoom(3)

    # make an image
    im = cam.snapshot()


    # add an axes triad -- note if we do this, we HAVE to do draw
    # domain, otherwise the image is blank (likely a bug)
    #cam.draw_coordinate_vectors(im)

    # add the domain box to the image:
    #nim = cam.draw_domain(im)
    nim = im
    # increase the contrast -- for some reason, the enhance default
    # to save_annotated doesn't do the trick (likely a bug)
    max_val = im[:,:,:3].std() * 6.0
    nim[:,:,:3] /= max_val

    f = pylab.figure()

    pylab.text(0.2, 0.85, "{:.3g} s".format(float(ds.current_time.d)),
               transform=f.transFigure, color="white")

    cam._render_figure = f

    cam.save_image(im, "test.png")
    
    # save annotated -- this added the transfer function values, 
    # but this messes up our image size defined above
    cam.save_annotated("{}.png".format(os.path.normpath(plotfile)), 
                       im, 
                       dpi=200, clear_fig=False)



if __name__ == "__main__":

    # Choose a field
    plotfile = ""


    try: plotfile = sys.argv[1]
    except: sys.exit("ERROR: no plotfile specified")

    doit(plotfile)


        
