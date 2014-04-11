# this example comes from 
# http://yt-project.org/doc/visualizing/volume_rendering.html

from yt.mods import *
import math
import yt.visualization.volume_rendering.api
import yt.units as u

pf = load("xrb-3d/plt58827")

# Choose a field
field = ('gas', 'velocity_magnitude')

# Do you want the log of the field?
use_log = True

# Find the bounds in log space of for your field
dd = pf.h.all_data()
mi, ma = dd.quantities["Extrema"](field)

pf.field_info[field].take_log = use_log

if use_log:
    vals = [1.e5, 3.16e5, 1.e6, 3.16e6, 1.e7]
else:
    vals = np.mgrid[mi:ma:7j].d

mi = min(vals)
ma = max(vals)

if use_log:
    mi, ma = np.log10(mi), np.log10(ma)

print mi, ma

# Instantiate the ColorTransferfunction.
tf =  yt.visualization.volume_rendering.api.ColorTransferFunction((mi, ma))

# Set up the camera parameters: center, looking direction, width, resolution
c = (pf.domain_right_edge + pf.domain_left_edge)/2.0
L = np.array([-0.0, -1.0, -1.0])
W = 3.0*pf.domain_width
N = 800

# Now let's add some isocontours, and take a snapshot, saving the image
# to a file.
#tf.add_layers(25, 0.3, colormap = 'spectral')
cm = "gist_rainbow"


for v in vals:
    if (use_log):
        tf.sample_colormap(math.log10(v), 0.002, colormap=cm) #, alpha=0.2)
    else:
        tf.sample_colormap(v, 0.01, colormap=cm) #, alpha=0.2)


#tf.add_layers(10,0.01,colormap='RdBu_r')


# Create a camera object
cam = pf.h.camera(c, L, W, N, tf, 
                  no_ghost=False,
                  fields = [field], log_fields = [use_log])

# make an image
im = cam.snapshot()
max_val = im[:,:,:3].std() * 4.0
im[:,:,:3] /= max_val

im.write_png('test-1.png')

# add an axes triad
cam.draw_coordinate_vectors(im)

# colorbar
cam.show_tf()

# add the domain box to the image:
nim = cam.draw_domain(im)
nim[:,:,:3] /= max_val

# save
nim.write_png('test_rendering_with_domain.png')

# save annotated
cam.save_annotated("test.png", nim)
