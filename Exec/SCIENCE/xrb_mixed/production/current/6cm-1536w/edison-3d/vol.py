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

    cm = "gist_rainbow"


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



    dd = ds.all_data()

    mi = min(vals)
    ma = max(vals)
    
    if use_log:
        mi, ma = np.log10(mi), np.log10(ma)
        

    # Instantiate the ColorTransferfunction.
    tf =  vr.ColorTransferFunction((mi, ma))

    # Set up the camera parameters: center, looking direction, width, resolution
    c = (ds.domain_right_edge + ds.domain_left_edge)/2.0
    L = np.array([1.0, 1.0, 1.0])
    L = np.array([1.0, 1.0, 1.2])
    W = 2.0*ds.domain_width
    N = 720

    north=[0.0,0.0,1.0]

    for v in vals:
        if (use_log):
            tf.sample_colormap(math.log10(v), sigma**2, colormap=cm) #, alpha=0.2)
        else:
            tf.sample_colormap(v, sigma**2, colormap=cm) #, alpha=0.2)


    # alternate attempt
    ds.periodicity = (True, True, True)

    # Create a camera object
    cam = vr.Camera(c, L, W, N, transfer_function=tf, ds=ds, 
                    no_ghost=False,
                    north_vector=north,
                    fields = [field], log_fields = [use_log])

    # make an image
    im = cam.snapshot()


    # add an axes triad -- note if we do this, we HAVE to do draw
    # domain, otherwise the image is blank (likely a bug)
    cam.draw_coordinate_vectors(im)

    # add the domain box to the image:
    nim = cam.draw_domain(im)

    # increase the contrast -- for some reason, the enhance default
    # to save_annotated doesn't do the trick (likely a bug)
    max_val = im[:,:,:3].std() * 4.0
    nim[:,:,:3] /= max_val

    f = pylab.figure()

    pylab.text(0.2, 0.15, "{:.3g} s".format(float(ds.current_time.d)),
               transform=f.transFigure, color="white")

    cam._render_figure = f
    
    # save annotated -- this added the transfer function values, 
    # but this messes up our image size defined above
    cam.save_annotated("{}_{}.png".format(os.path.normpath(plotfile), fname), 
                       nim, 
                       dpi=145, clear_fig=False)



if __name__ == "__main__":

    # Choose a field
    fname = "vz"
    plotfile = ""


    try: plotfile = sys.argv[1]
    except: sys.exit("ERROR: no plotfile specified")

    try: fname = sys.argv[2]
    except: pass

    doit(plotfile, fname)


        
