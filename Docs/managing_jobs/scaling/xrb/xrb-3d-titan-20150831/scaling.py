import math
import numpy as np
import matplotlib.pyplot as plt 

def scaling():

    data = np.loadtxt("scaling.txt")

    N_mpi = data[:,0]
    N_omp = data[:,1]
    
    t_adv = data[:,2]
    t_MAC = data[:,3]
    t_nodal = data[:,4]
    t_react = data[:,5]
    t_misc = data[:,6]
    t_total = data[:,7]
    std_total = data[:,8]
    width = data[:,9]

    cores = N_mpi * N_omp

    
    # we ran 2 tests, one that was 384 zones wide, and one that was 768 zones wide
    idx_384 = width[:] == 384
    idx_768 = width[:] == 768

    plt.errorbar(cores[idx_384], t_total[idx_384], yerr=std_total[idx_384], marker="o", color="k", label=r"$384\times 384\times 768$")
    plt.errorbar(cores[idx_768], t_total[idx_768], yerr=std_total[idx_768], marker="^", color="r", label=r"$768\times 768\times 768$")

    ax = plt.gca()
    ax.set_xscale("log")
    ax.set_yscale("log")

    # ideal
    cm = np.min(cores[idx_384])
    cM = np.max(cores[idx_384])
    id = np.argmin(cores[idx_384])
    plt.loglog([cm, cM], t_total[idx_384][id]*cm/np.array([cm, cM]), ":", color="k")

    cm = np.min(cores[idx_768])
    cM = np.max(cores[idx_768])
    id = np.argmin(cores[idx_768])
    plt.loglog([cm, cM], t_total[idx_768][id]*cm/np.array([cm, cM]), ":", color="k")

    plt.text(600, 1.25, "Cray 8.4.0 compilers; 2015-08-31")
    plt.xlim(512, 131072)

    plt.legend(frameon=False)

    plt.xlabel("number of cores")
    plt.ylabel("average time to advance timestep")

    plt.title("OLCF Titan Scaling for 3-d XRB")

    plt.tight_layout()

    plt.savefig("titan_xrb_scaling.png")


    plt.loglog(cores[idx_384], t_react[idx_384], "o--", color="0.5", label=r"$384\times 384\times 768$ reactions")
    plt.loglog(cores[idx_768], t_react[idx_768], "^--", color="m", label=r"$384\times 384\times 768$ reactions")

    plt.savefig("titan_xrb_scaling_react.png")

if __name__== "__main__":
    scaling()

