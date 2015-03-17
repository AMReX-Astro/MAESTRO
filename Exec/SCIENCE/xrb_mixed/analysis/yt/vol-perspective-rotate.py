#!/usr/bin/env python

import matplotlib
matplotlib.use('agg')

# this example comes from 
# http://yt-project.org/doc/visualizing/volume_rendering.html

import math
import sys
import numpy as np
import os
import matplotlib.pyplot as plt

import yt
import yt.visualization.volume_rendering.api as vr

def doit(plotfile, fname, view):

    ds = yt.load(plotfile)

    cm = "gist_rainbow"

    # rotate?
    Nrot = 0


    if fname == "vz":
        field = ('gas', 'velocity_z')
        use_log = False
        
        vals = [-1.e7, -5.e6, 5.e6, 1.e7]
        sigma = 5.e5

        cm = "coolwarm"
        
    elif fname == "magvel":
        field = ('gas', 'velocity_magnitude')
        use_log = False
        
        vals = [1.e6, 2.e6, 4.e6, 8.e6, 1.6e7]
        sigma = 2.e5

    elif fname == "radvel":
        field = ('boxlib', 'radial_velocity')
        use_log = False
        
        vals = [-1.e7, -5.e6, -2.5e6, 2.5e6, 5.e6, 1.e7]
        sigma = 2.e5
        
        cm = "coolwarm"


    dd = ds.all_data()

    mi = min(vals)
    ma = max(vals)
    
    if use_log:
        mi, ma = np.log10(mi), np.log10(ma)
        

    # Instantiate the ColorTransferfunction.
    tf =  vr.ColorTransferFunction((mi, ma))

    # Set up the camera parameters: center, looking direction, width, resolution
    #c = np.array([0.0, 0.0, 0.0])
    center = 0.5*(ds.domain_left_edge + ds.domain_right_edge)

    xmax, ymax, zmax = ds.domain_right_edge

    # this shows something face on
    c = np.array([-4.0*xmax, center[1], center[2]])
    L = np.array([1.0, 0.0, 0.0])

    # this shows something face on
    c = np.array([-2.0*xmax, -2.0*ymax, center[2]])
    L = np.array([1.0, 1.0, 0.0])

    if view == "top":
        # this is a view from above
        c = np.array([-2.0*xmax, -2.0*ymax, center[2]+0.75*zmax])

        # the normal vector should be pointing back through the center
        L = center.d - c
        L = L/np.sqrt((L**2).sum())

    W = 3.0*ds.domain_width
    N = 1080

    north=[0.0,0.0,1.0]

    for v in vals:
        if (use_log):
            tf.sample_colormap(math.log10(v), sigma**2, colormap=cm) #, alpha=0.2)
        else:
            tf.sample_colormap(v, sigma**2, colormap=cm) #, alpha=0.2)


    # alternate attempt
    ds.periodicity = (True, True, True)

    
    # loop doing a Matrix transform on the orientation vectors about
    # the domain center to simulate rotation
    if not Nrot == 0:
        dtheta = 2.0*np.pi/Nrot
    
    for n in range(max(1,Nrot)):


        # Create a camera object
        if n == 0:
            cam = vr.PerspectiveCamera(c, L, W, N, transfer_function=tf, ds=ds, 
                                       no_ghost=False,
                                       north_vector=north,
                                       fields = [field], log_fields = [use_log])
    

        if not Nrot == 0:
            cam.yaw(dtheta, center)
        

        # make an image
        im = cam.snapshot()


        # add an axes triad -- note if we do this, we HAVE to do draw
        # domain, otherwise the image is blank (likely a bug)
        cam.draw_coordinate_vectors(im)

        # add the domain box to the image:
        nim = cam.draw_domain(im)

        # increase the contrast -- for some reason, the enhance default
        # to save_annotated doesn't do the trick (likely a bug)
        if fname == "vz":
            # this normalization comes from looking at im[:,:,:3].std() * 4.0 for
            # a 3-d wide XRB visualization near 0.02 s
            max_val = 0.0276241228025

        elif n == 0:
            max_val = im[:,:,:3].std() * 4.0

        nim[:,:,:3] /= max_val

        f = plt.figure()

        plt.text(0.2, 0.15, "{:.3g} s".format(float(ds.current_time.d)),
                   transform=f.transFigure, color="white")

        cam._render_figure = f
    
        # save annotated -- this added the transfer function values, 
        # but this messes up our image size defined above
        #cam.save_annotated("{}_{}_{:03d}.png".format(os.path.normpath(plotfile), fname, n), 
        #                   nim, 
        #                   dpi=145, clear_fig=False)

        if view == "top":
            ostr = fname + "_top"
        else:
            ostr = fname
            
        cam.save_annotated("{}_{}_HD{:03d}.png".format(os.path.normpath(plotfile), ostr, n), 
                           nim, 
                           dpi=218, clear_fig=False)
        
        plt.close()


if __name__ == "__main__":

    # Choose a field
    fname = "vz"
    plotfile = ""

    view = "normal"
    
    try: plotfile = sys.argv[1]
    except: sys.exit("ERROR: no plotfile specified")

    try: fname = sys.argv[2]
    except: pass

    try: view = sys.argv[3]
    except: pass
    
    doit(plotfile, fname, view)


        
