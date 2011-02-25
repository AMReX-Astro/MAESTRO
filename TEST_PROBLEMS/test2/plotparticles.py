#!/usr/bin/env python

import math
import numpy
import sys
import parseparticles
import pylab


#-----------------------------------------------------------------------------
def main(files):

    # this returns a dict whose keys are a unique identifier (based on 
    # id and CPU) and values are the actual particle objects
    particlesDict = parseparticles.parseParticleFile(files)

    # get just the particle objects
    particles = particlesDict.values()


    # print out some info and test the routines
    print "number of unique particles = ", len(particles)


    # make a plot of the particle paths
    pylab.clf()

    n = 0
    while (n < len(particles)):

        # get numpy arrays containing the time and coordinate
        # information for particle 0
        coords, time = particles[n].flatten()

        pylab.scatter([coords[0,0]], [coords[1,0]], marker="x")
        pylab.plot(coords[0,:], coords[1,:])

        n += 1

    pylab.xlabel("x")
    pylab.ylabel("y")

    a = pylab.gca()
    a.set_aspect("equal")

    pylab.savefig("particle_paths.png")


    # make an animation -- note: this assumes that all particles exist
    # at all timesteps
    nstep = 0
    while (nstep < len(particles[0].history)):
        
        pylab.clf()
        
        n = 0
        while (n < len(particles)):

            # plot the position of the current particle (n) at the
            # current step (nstep).  Color it by the abundance of
            # X(Mg24).
            XMg24_min = 1.e-6
            XMg24_max = 1.e-5

            pylab.scatter([particles[n].history[nstep].xyz[0]],
                          [particles[n].history[nstep].xyz[1]] ,
                          marker="o", s=1.0,
                          c=math.log10(particles[n].history[nstep].data["magnesium-24"]/
                                       particles[n].history[nstep].data["density"]),
                          vmin=math.log10(XMg24_min),
                          vmax=math.log10(XMg24_max),
                          edgecolor="None")

            n += 1

        # axis labels
        pylab.xlabel("x")
        pylab.ylabel("y")

        # color bar
        cb = pylab.colorbar(orientation="horizontal")
        cb.set_label(r"log$_{10}$[X($^{24}$Mg)]")

        pylab.axis([2.e7,2.2e8,4.e7,1.2e8])
        a = pylab.gca()
        a.set_aspect("equal")

        a.xaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))
        a.yaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))


        f = pylab.gcf()
        f.set_size_inches(8.0,6.0)

        pylab.savefig("particles_%04d.png" % (nstep) )

        nstep += 1


#-----------------------------------------------------------------------------
if __name__== "__main__":

    if (len(sys.argv) == 1):
        print "ERROR: no particle data files specified\n"
        sys.exit(2)

    main(sys.argv[1:])




        
