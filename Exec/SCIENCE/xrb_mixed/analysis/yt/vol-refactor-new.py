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
from yt.visualization.volume_rendering.api import Scene, Camera, VolumeSource, \
                                                  ColorTransferFunction


def doit(plotfile, fname):

    ds = yt.load(plotfile)

    cm = "gist_rainbow"                                                         

    if fname == "vz":
        field = ('gas', 'velocity_z')
        use_log = False
        
        vals = [-1.e7, -5.e6, 5.e6, 1.e7]
        sigma = 5.e5

        fmt = None

        cm = "coolwarm"                                                         


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
    


    # this is hackish, but there seems to be no better way to set whether
    # you are rendering logs
    ds._get_field_info(field).take_log = use_log

    # hack periodicity
    ds.periodicity = (True, True, True)

    mi = min(vals)
    ma = max(vals)
    
    if use_log:
        mi, ma = np.log10(mi), np.log10(ma)


    print mi, ma

        
    # setup the scene and camera
    sc = Scene()
    cam = Camera(ds, lens_type="perspective")

    # Set up the camera parameters: center, looking direction, width, resolution
    center = (ds.domain_right_edge + ds.domain_left_edge)/2.0
    xmax, ymax, zmax = ds.domain_right_edge                                     

    # this shows something face on                                              
    c = np.array([-8.0*xmax, center[1], center[2]])  

    # the normal vector should be pointing back through the center              
    L = center.d - c                                                            
    L = L/np.sqrt((L**2).sum())    

    north_vector=[0.0,0.0,1.0]

    cam.position = ds.arr(c)
    cam.switch_orientation(normal_vector=L, north_vector=north_vector)

    cam.set_width(ds.domain_width*4)

    cam.resolution = (720,720)

    # create the volume source
    vol = VolumeSource(ds, field=field)
    
    # Instantiate the ColorTransferfunction.
    tf = vol.transfer_function
    tf = ColorTransferFunction((mi, ma))
    #tf.grey_opacity=True                                  


    for v in vals:
        if use_log:
            tf.sample_colormap(math.log10(v), sigma**2, colormap=cm) #, alpha=0.2)
        else:
            print v
            tf.sample_colormap(v, sigma**2, colormap=cm) #, alpha=0.2)

    sc.camera = cam
    sc.add_source(vol)
    sc.render("test_perspective.png", clip_ratio=6.0)


if __name__ == "__main__":

    # Choose a field
    fname = "vz"
    plotfile = ""


    try: plotfile = sys.argv[1]
    except: sys.exit("ERROR: no plotfile specified")

    try: fname = sys.argv[2]
    except: pass

    doit(plotfile, fname)


        
