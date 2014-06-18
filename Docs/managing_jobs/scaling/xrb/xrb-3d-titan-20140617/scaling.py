import math
import numpy as np
import pylab 

def scaling():

    data = np.loadtxt("scaling.txt")

    problem_size = 384*384*768

    print data.shape

    N_mpi = data[:,0]
    N_omp = data[:,1]
    
    t_adv = data[:,2]
    t_MAC = data[:,3]
    t_nodal = data[:,4]
    t_react = data[:,5]
    t_misc = data[:,6]
    t_total = data[:,7]
    jmode = data[:,8]

    cores = N_mpi * N_omp

    # we ran in 3 batches, with 4, 8, and 16 threads -- separate them out
    # we also ran in 2 modes on titan -- with 8 cores per compute node and 16
    idx_omp1_j1 = np.logical_and(N_omp[:] == 1, jmode[:] == 1)
    idx_omp4_j1 = np.logical_and(N_omp[:] == 4, jmode[:] == 1)
    idx_omp8_j1 = np.logical_and(N_omp[:] == 8, jmode[:] == 1)
    idx_omp16_j1 = np.logical_and(N_omp[:] == 16, jmode[:] == 1)

    idx_omp1_j2 = np.logical_and(N_omp[:] == 1, jmode[:] == 2)
    idx_omp4_j2 = np.logical_and(N_omp[:] == 4, jmode[:] == 2)
    idx_omp8_j2 = np.logical_and(N_omp[:] == 8, jmode[:] == 2)
    idx_omp16_j2 = np.logical_and(N_omp[:] == 16, jmode[:] == 2)

    pylab.loglog(cores[idx_omp1_j1], t_total[idx_omp1_j1], "o-", color="k", label="MPI")
    pylab.loglog(cores[idx_omp4_j1], t_total[idx_omp4_j1], "o-", color="b", label="MPI + 4 OpenMP threads")
    pylab.loglog(cores[idx_omp8_j1], t_total[idx_omp8_j1], "o-", color="r", label="MPI + 8 OpenMP threads")
    pylab.loglog(cores[idx_omp16_j1], t_total[idx_omp16_j1], "o-", color="g", label="MPI + 16 OpenMP threads")

    pylab.loglog(cores[idx_omp1_j2], t_total[idx_omp1_j2], "^-", color="k")
    pylab.loglog(cores[idx_omp4_j2], t_total[idx_omp4_j2], "^-", color="b")
    pylab.loglog(cores[idx_omp8_j2], t_total[idx_omp8_j2], "^-", color="r")
    pylab.loglog(cores[idx_omp16_j2], t_total[idx_omp16_j2], "^-", color="g")

    # ideal
    cm = np.min(cores)
    cM = np.max(cores)
    id = np.argmin(cores)
    pylab.loglog([cm, cM], t_total[id]*cm/np.array([cm, cM]), ":", color="k", label="ideal scaling")

    pylab.legend(frameon=False)

    pylab.xlabel("number of cores")
    pylab.ylabel("average time to advance timestep")

    pylab.title("OLCF Titan Scaling for 3-d XRB (384 x 384 x 768 zones)")

    pylab.ylim(1.,200.)

    pylab.savefig("xrb_titan_scaling_by_parallel.png")

    # we also ran with 3 different grid sizes: 32^3, 48^3, and 64^3.
    # this is essentially N_mpi
    pylab.clf()

    grids = np.unique(N_mpi)
    jtype = [1, 2]
    for j in jtype:
        colors = ["g", "b", "r"]
        for g in grids:
            idx = np.logical_and(N_mpi[:] == g, jmode[:] == j)
            gsize = int(round((problem_size/g)**(1./3.)))
            if j == 1:
                pylab.loglog(cores[idx], t_total[idx], "o-", color=colors.pop(),
                             label=r"${}^3$ grid".format(gsize))
            else:
                pylab.loglog(cores[idx], t_total[idx], "^-", color=colors.pop())

    # ideal
    cm = np.min(cores)
    cM = np.max(cores)
    id = np.argmin(cores)

    pylab.loglog([cm, cM], t_total[id]*cm/np.array([cm, cM]), ":", color="k", label="ideal scaling")

    pylab.legend(frameon=False)

    pylab.xlabel("number of cores")
    pylab.ylabel("average time to advance timestep")

    pylab.title("OLCF Titan Scaling for 3-d XRB (384 x 384 x 768 zones)")

    pylab.ylim(1.,200.)
    pylab.savefig("xrb_titan_scaling_by_grid.png")



if __name__== "__main__":
    scaling()

