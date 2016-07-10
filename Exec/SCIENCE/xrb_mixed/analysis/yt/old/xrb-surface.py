import sys

import yt

import numpy as np
import pylab

from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

print sys.argv

ds = yt.load(sys.argv[1])

# a hack to remove one zone in each direction to prevent boundaries from coming into play
root_dds = ds.domain_width/ds.domain_dimensions                              
half_w = ds.domain_width/2. - root_dds                                    
half_w[:] -= root_dds[:]                                                   
reg = ds.region(ds.domain_center, ds.domain_center-half_w, ds.domain_center+half_w)  

ds.periodicity = (True, True, True)
dd = ds.all_data()

# figure
fig = pylab.figure()
ax = fig.gca(projection='3d')
ax.set_aspect('equal')



# we'll do several different isocontours
vals = [-1.e7, -5.e6, 5.e6, 1.e7]

# colors blue -> red as rgb
colors = [(255, 0, 0), (255, 100, 100), 
          (100, 100, 255), (150, 150, 255)]

for n, v in enumerate(vals):

    print "working on value: {}".format(v)

    v = vals[n]
    c = colors[n]
    c = np.array([c[0], c[1], c[2], 25.6])/255.0  # add alpha

    surface = ds.surface(reg, "velocity_z", v)
    
    p3dc = Poly3DCollection(surface.triangles, linewidth=0.0)

    p3dc.set_facecolors(c)
    ax.add_collection(p3dc)


# trick -- put a max bounding box so we get equal aspect
# see: http://stackoverflow.com/questions/13685386/matplotlib-equal-unit-length-with-equal-aspect-ratio-z-axis-is-not-equal-to
max_range = max(ds.domain_right_edge.d - ds.domain_left_edge.d)/2.0

center = 0.5*(ds.domain_right_edge.d - ds.domain_left_edge.d)

# put all the triangls into view
ax.set_xlim(center[0] - max_range, center[0] + max_range)
ax.set_ylim(center[1] - max_range, center[1] + max_range)
ax.set_zlim(center[2] - max_range, center[2] + max_range)
                  
pylab.savefig("xrb_surf.png")

