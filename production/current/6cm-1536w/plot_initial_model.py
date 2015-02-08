# plot the initial model, showing the location of various cutoffs

import numpy as np
import matplotlib.pyplot as plt

# read the initial model
model = np.loadtxt("toy_xrb.hi_dens.hse.tanh.delta_12.000cm.dx_6.000cm.CNO")

ir = 0    # position
idens = 1 # density
itemp = 2 # temperature
ipres = 3 # pressure


plt.semilogy(model[:,ir], model[:,idens])

sponge_start_factor = 25.0

dens = 5.e2
plt.semilogy(model[:,ir], dens*np.ones_like(model[:,ir]), ls=":", color="r")
plt.semilogy(model[:,ir], sponge_start_factor*dens*np.ones_like(model[:,ir]), ls="--", color="r")

dens = 1.e3
plt.semilogy(model[:,ir], dens*np.ones_like(model[:,ir]), ls=":", color="g")
plt.semilogy(model[:,ir], sponge_start_factor*dens*np.ones_like(model[:,ir]), ls="--", color="g")

dens = 5.e3
plt.semilogy(model[:,ir], dens*np.ones_like(model[:,ir]), ls=":", color="b")
plt.semilogy(model[:,ir], sponge_start_factor*dens*np.ones_like(model[:,ir]), ls="--", color="b")

plt.ylim(1.0, 1.e8)

plt.savefig("initial_model.png")

