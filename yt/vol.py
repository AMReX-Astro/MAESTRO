#!/usr/bin/env python

import matplotlib
matplotlib.use('agg')

# this example comes from 
# http://yt-project.org/doc/visualizing/volume_rendering.html

import math
import sys
import pylab

from yt.mods import *
import yt.visualization.volume_rendering.api as vr

def doit(plotfile, fname):

    ds = load(plotfile)

    if fname == "vz":
        field = ('gas', 'velocity_z')
        use_log = False
        
        vals = [-1.e7, -5.e6, 5.e6, 1.e7]
        sigma = 5.e5
        
    elif fname == "magvel":
        field = ('gas', 'velocity_magnitude')
        use_log = False
        
        vals = [1.e6, 2.5e6, 5.e6, 7.5e6, 1.e7]
        sigma = 1.e5


    dd = ds.h.all_data()

    ds.field_info[field].take_log = use_log


    mi = min(vals)
    ma = max(vals)
    
    if use_log:
        mi, ma = np.log10(mi), np.log10(ma)
        

    # Instantiate the ColorTransferfunction.
    tf =  vr.ColorTransferFunction((mi, ma))

    # Set up the camera parameters: center, looking direction, width, resolution
    c = np.array([0.0, 0.0, 0.0])
    L = np.array([1.0, 1.0, 1.0])
    L = np.array([1.0, 1.0, 1.2])
    W = 1.5*ds.domain_width
    N = 720

    north=[0.0,0.0,-1.0]

    cm = "gist_rainbow"

    for v in vals:
        if (use_log):
            tf.sample_colormap(math.log10(v), sigma**2, colormap=cm) #, alpha=0.2)
        else:
            tf.sample_colormap(v, sigma**2, colormap=cm) #, alpha=0.2)


    # alternate attempt
    ds.periodicity = (True, True, True)

    # Create a camera object
    cam = vr.Camera(c, L, W, N, transfer_function=tf, ds=ds, 
                    no_ghost=False, #data_source=reg,
                    north_vector=north,
                    fields = [field], log_fields = [use_log])

    #cam.rotate(3.0*np.pi/2., rot_vector=rot_vector)

    # make an image
    im = cam.snapshot()

    # add an axes triad
    cam.draw_coordinate_vectors(im)

    # colorbar
    cam.show_tf()

    # add the domain box to the image:
    nim = cam.draw_domain(im)

    # increase the contrast
    max_val = im[:,:,:3].std() * 4.0
    nim[:,:,:3] /= max_val

    #nim.write_png("xrb_vol_test.png")

    f = pylab.figure()

    pylab.text(0.2, 0.85, "{:.3g} s".format(float(ds.current_time.d)),
               transform=f.transFigure, color="white")

    cam._render_figure = f
    
    # save annotated -- this added the transfer function values, 
    # but this messes up our image size defined above
    cam.save_annotated("xrb_vol_{}.png".format(fname), nim, dpi=145, clear_fig=False)



if __name__ == "__main__":

    # Choose a field
    fname = "vz"
    plotfile = ""


    try: plotfile = sys.argv[1]
    except: sys.exit("ERROR: no plotfile specified")

    try: fname = sys.argv[2]
    except: pass

    doit(plotfile, fname)


        
