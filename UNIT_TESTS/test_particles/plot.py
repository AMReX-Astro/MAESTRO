#!/usr/bin/env python

import math
import numpy
import sys
import parseparticles
import pylab


#-----------------------------------------------------------------------------
def main(files):

    # domain information
    xcenter = 0.5
    ycenter = 0.5


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

    pylab.axis([0,1,0,1])

    a = pylab.gca()
    a.set_aspect("equal")

    pylab.savefig("particle_paths.png")


    # compute the relative change in the radius of the particle from
    # start to end -- this is the error in the integration.
    pylab.clf()

    # assume that all particles were around for all timesteps
    nstep = len(particles[0].history)

    error = numpy.zeros(len(particles), dtype=numpy.float64)
    
    n = 0
    while (n < len(particles)):

        # compute the relative change in radius
        r_start = math.sqrt( (particles[n].history[0].xyz[0] - xcenter)**2 +
                             (particles[n].history[0].xyz[1] - ycenter)**2)

        r_end = math.sqrt( (particles[n].history[nstep-1].xyz[0] - xcenter)**2 +
                           (particles[n].history[nstep-1].xyz[1] - ycenter)**2)

        error[n] = math.fabs(r_start - r_end)/r_start

        print "particle: (%d, %d), r init = %f, rel error = %g" % \
            (particles[n].pid, particles[n].originCPU, r_start, error[n])
        
        n += 1


    print " "
    print "maximum error = ", max(error)
    print "L2 norm error = ", math.sqrt(numpy.dot(error,error)/len(particles))


#-----------------------------------------------------------------------------
if __name__== "__main__":

    if (len(sys.argv) == 1):
        print "ERROR: no particle data files specified\n"
        sys.exit(2)

    main(sys.argv[1:])




        
