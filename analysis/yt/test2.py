# this example comes from 
# http://yt-project.org/doc/visualizing/volume_rendering.html

from yt.mods import *

pf = load("plt23437")

# Choose a field
field = 'magvel'

# Do you want the log of the field?
use_log = True

# Find the bounds in log space of for your field
dd = pf.h.all_data()
mi, ma = dd.quantities["Extrema"](field)[0]

pf.field_info[field].take_log = use_log

print mi, ma

if use_log:
    mi,ma = np.log10(mi), np.log10(ma)

print mi, ma

# Instantiate the ColorTransferfunction.
tf = ColorTransferFunction((mi, ma))

# Set up the camera parameters: center, looking direction, width, resolution
c = (pf.domain_right_edge + pf.domain_left_edge)/2.0
L = np.array([-0.0, -1.0, -1.0])
W = 3.0*pf.domain_width
N = 800

# Now let's add some isocontours, and take a snapshot, saving the image
# to a file.
#tf.add_layers(25, 0.3, colormap = 'spectral')
tf.sample_colormap(3, 0.02)
tf.sample_colormap(4, 0.02)
tf.sample_colormap(5, 0.02)
tf.sample_colormap(6, 0.02)
tf.sample_colormap(7, 0.02)

# Create a camera object
cam = pf.h.camera(c, L, W, N, tf, 
                  no_ghost=False,
                  fields = [field], log_fields = [use_log])

im = cam.snapshot('test_rendering.png')

# To add the domain box to the image:
nim = cam.draw_domain(im)
nim.write_png('test_rendering_with_domain.png')

