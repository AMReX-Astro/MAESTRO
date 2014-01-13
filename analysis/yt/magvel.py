from yt.mods import *

pf = load("plt23437")
field = "z_vel"

dd = pf.h.all_data()
mi, ma = dd.quantities["Extrema"](field)[0]
mi -= 0.1 ; ma += 0.1 # To allow a bit of room at the edges

print mi, ma

tf = ColorTransferFunction((mi, ma))
tf.add_layers(8, w=0.01)
c = (pf.domain_right_edge + pf.domain_left_edge)/2.0
L = na.array([0.0, 1.0, -1.0])
W = 0.5 / pf["unitary"]

N = 512

cam = pf.h.camera(c, L, W, N, tf, fields=[field])
fn = "%s_image.png" % pf

cam.snapshot(fn)
