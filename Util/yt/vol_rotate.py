# this example comes from 
# http://yt-project.org/doc/visualizing/volume_rendering.html

from yt.mods import *
import math
import yt.visualization.volume_rendering.api

pf = load("plt23437")

# Choose a field
fname = "vz"

if fname == "vz":
    field = ('gas', 'velocity_z')
    use_log = False

    vals = [-1.e7, -5.e6, 5.e6, 1.e7]
    sigma = 5.e5

elif fname == "magvel":
    field = ('gas', 'velocity_magnitude')
    use_log = True

    vals = [1.e5, 3.16e5, 1.e6, 3.16e6, 1.e7]
    sigma = 0.1


dd = pf.h.all_data()

pf.field_info[field].take_log = use_log


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
W = 2.0*pf.domain_width
N = 720

rot_vector=[0.0,0.0,1.0]

cm = "gist_rainbow"

for v in vals:
    if (use_log):
        tf.sample_colormap(math.log10(v), sigma**2, colormap=cm) #, alpha=0.2)
    else:
        tf.sample_colormap(v, sigma**2, colormap=cm) #, alpha=0.2)


# an attempt to get around the periodic assumption -- suggested by Nathan
#root_dds = pf.domain_width/pf.domain_dimensions
#half_w = pf.domain_width/2.# - root_dds
#half_w[2] -= root_dds[2]
#reg = pf.region(pf.domain_center, pf.domain_center-half_w, pf.domain_center+half_w)


# alternate attempt
pf.periodicity = (True, True, True)

N_frames = 360

for f in range(N_frames):

    # Create a camera object
    cam = pf.h.camera(c, L, W, N, tf, 
                      no_ghost=False, #data_source=reg,
                      fields = [field], log_fields = [use_log])

    cam.rotate(2*np.pi*f/N_frames, rot_vector=rot_vector)

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

    nim.write_png("xrb_vol_{:04}.png".format(f))

    # save annotated -- this added the transfer function values, 
    # but this messes up our image size defined above
    cam.save_annotated("test_{:04}.png".format(f), nim, dpi=145)
