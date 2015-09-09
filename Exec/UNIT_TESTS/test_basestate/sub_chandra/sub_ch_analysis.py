import numpy as np
import matplotlib.pyplot as plt

castro = np.loadtxt("castro_hse_adjust_subchandra.out")
maestro_orig = np.loadtxt("base.orig")
maestro_new = np.loadtxt("base.new")

plt.subplot(211)

plt.plot(maestro_orig[:,0], maestro_orig[:,1], color="0.5", label="initial conditions")
plt.plot(castro[:,0], castro[:,1], color="k", label="Castro")
plt.scatter(maestro_new[:,0], maestro_new[:,1], color="r", label="Maestro", marker="x", s=7)

plt.legend(frameon=False, fontsize="small", loc="best")

ax = plt.gca()
ax.xaxis.set_major_formatter(plt.ScalarFormatter(useMathText=True))
ax.yaxis.set_major_formatter(plt.ScalarFormatter(useMathText=True))
        
ax.set_yscale("log")

plt.xlim(0., 7.e8)
plt.ylim(100.,1.e8)

plt.xlabel("r [cm]")
plt.ylabel(r"$\rho$ [g/cc]")


plt.subplot(212)

plt.plot(maestro_orig[:,0], maestro_orig[:,2], color="0.5", label="initial conditions")
plt.plot(castro[:,0], castro[:,5], color="k", label="Castro")
plt.scatter(maestro_new[:,0], maestro_new[:,2], color="r", label="Maestro", marker="x", s=7)

plt.legend(frameon=False, fontsize="small", loc="best")

ax = plt.gca()
ax.xaxis.set_major_formatter(plt.ScalarFormatter(useMathText=True))
ax.yaxis.set_major_formatter(plt.ScalarFormatter(useMathText=True))

plt.xlabel("r [cm]")
plt.ylabel(r"$T$ [K]")

plt.xlim(0., 7.e8)
plt.ylim(0.0, 5.e8)

plt.suptitle("test_basestate for the sub_chandra initial conditions", fontsize="medium")

f = plt.gcf()
f.set_size_inches(6.0, 7.5)

plt.savefig("test_basestate_subch.png")
