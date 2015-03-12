import yt
import numpy as np

Nx = 8
Ny = 8
Nz = 16

xmin = ymin = zmin = 0.0
xmax = ymax = 1.0
zmax = 2.0

center = np.array( [0.5*(xmin + xmax), 0.5*(ymin + ymax), 0.5*(zmin + zmax)] )

arr = np.zeros((Nx,Ny,Nz), dtype=np.float64)
arr[:,:,:] = 1.0

bbox = np.array([ [xmin, xmax], [ymin, ymax], [zmin, zmax] ])

# coordinates -- in the notation data[i, j, k]
dcoord = (xmax - xmin)/Nx
x = np.linspace(xmin+dcoord, xmax-dcoord, Nx, endpoint=True)
y = np.linspace(ymin+dcoord, ymax-dcoord, Ny, endpoint=True)
z = np.linspace(zmin+dcoord, zmax-dcoord, Nz, endpoint=True)

x3d, y3d, z3d = np.meshgrid(x, y, z, indexing="ij")


data = dict(Density = arr)
ds = yt.load_uniform_grid(data, arr.shape, bbox=bbox)

import yt.visualization.volume_rendering.api as vr

# looking from edge in midplane z
c = np.array([-5*xmax,-5*ymax,0.5*(zmin+zmax)])
L = np.array([1.0, 1.0, 0.0])
zoom = 4.0

# looking from below
c = np.array([-5*xmax,-5*ymax,-1.0*zmax])
L = np.array([1.0, 1.0, 0.5])
zoom = 1.0

W = zoom*ds.domain_width
N = 720

north = np.array([0, 0, 1])

tf = vr.ColorTransferFunction((0.1,1.0))
tf.sample_colormap(1.0, 0.01, colormap="coolwarm")

cam = vr.PerspectiveCamera(c, L, W, N, transfer_function=tf, 
                north_vector = north, 
                ds=ds, fields=[('gas', 'Density')], 
                log_fields=[False])

im = cam.snapshot()
cam.draw_coordinate_vectors(im)
nim = cam.draw_domain(im)

nim.write_png("test_{:03}.png".format(0))
