#!/usr/bin/env python

import matplotlib
matplotlib.use('agg')

import sys

import yt
import numpy as np
from yt.visualization.volume_rendering.api import \
    Scene, \
    VolumeSource


# this is for the wdconvect problem

def doit(plotfile):

    ds = yt.load(plotfile)
    ds.periodicity = (True, True, True)

    cm = "coolwarm"

    field = ('gas', 'velocity_z')
    ds._get_field_info(field).take_log = False
        
    sc = Scene()


    # add a volume: select a sphere
    #center = (0, 0, 0)
    #R = (5.e8, 'cm')

    #dd = ds.sphere(center, R)

    vol = VolumeSource(ds, field=field)
    sc.add_source(vol)


    # transfer function
    vals = [-1.e7, -5.e6, 5.e6, 1.e7]
    sigma = 5.e5

    tf =  yt.ColorTransferFunction((min(vals), max(vals)))

    tf.clear()
    cm = "coolwarm"
    for v in vals:
        tf.sample_colormap(v, sigma**2, colormap=cm) #, alpha=0.2)

    sc.get_source(0).transfer_function = tf

    cam = sc.add_camera(ds, lens_type="perspective")        
    cam.resolution = (1920, 1080)

    center = 0.5*(ds.domain_left_edge + ds.domain_right_edge)

    cam.position = [2.5*ds.domain_right_edge[0],
                    2.5*ds.domain_right_edge[1],
                    center[2]+0.25*ds.domain_right_edge[2]]
    
    # look toward the center -- we are dealing with an octant
    normal = (center - cam.position)
    normal /= np.sqrt(normal.dot(normal))

    cam.switch_orientation(normal_vector=normal,
                           north_vector=[0., 0., 1.])
    cam.set_width(ds.domain_width)

    sc.camera = cam
    #sc.annotate_axes(alpha=0.05)
    #sc.annotate_domain(ds, color=np.array([0.05, 0.05, 0.05, 0.05]))
    #sc.annotate_grids(ds, alpha=0.05)

    sc.render()
    sc.save("{}_radvel".format(plotfile), sigma_clip=4.0)
    sc.save_annotated("{}_radvel_annotated.png".format(plotfile),
                      sigma_clip=4.0, 
                      text_annotate=[[(0.05, 0.05), 
                                      "t = {}".format(ds.current_time.d),
                                      dict(horizontalalignment="left")],
                                     [(0.5,0.95), 
                                      "Maestro simulation of convection in a mixed H/He XRB",
                                      dict(color="y", fontsize="24",
                                           horizontalalignment="center")]])

if __name__ == "__main__":

    # Choose a field
    plotfile = ""


    try: plotfile = sys.argv[1]
    except: sys.exit("ERROR: no plotfile specified")

    doit(plotfile)


        
