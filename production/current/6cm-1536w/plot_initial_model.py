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

dens = [5.e2, 1.e3, 5.e3]
colors = ["r", "g", "b", "k", "0.5", "c"]

for n, d in enumerate(dens):
    plt.semilogy(model[:,ir],
                 d*np.ones_like(model[:,ir]),
                 ls=":", color=colors[n], label="{}".format(d))
    plt.semilogy(model[:,ir],
                 sponge_start_factor*d*np.ones_like(model[:,ir]),
                 ls="--", color=colors[n])

plt.ylim(1.0, 1.e8)

plt.xlabel("r [cm]")
plt.ylabel("density [g/cc]")

plt.legend(frameon=False)

plt.savefig("initial_model.png")

