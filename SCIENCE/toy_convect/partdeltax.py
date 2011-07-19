#!/usr/bin/env python

# compute the maximum delta x of a particle over a timestep -- for
# debugging purposes.  Try to account for periodic BCs in x

import math
import numpy
import sys
import parseparticles
import pylab


xmin = 0.0
xmax = 7.5e7
dx = (xmax - xmin)/float(512)


#-----------------------------------------------------------------------------
def main(files):

    # this returns a dict whose keys are a unique identifier (based on 
    # id and CPU) and values are the actual particle objects
    particlesDict = parseparticles.parseParticleFile(files)

    # get just the particle objects
    particles = particlesDict.values()


    # print out some info and test the routines
    print "number of unique particles = ", len(particles)


    if (not particles[0].dim == 2):
        print "ERROR: ony 2-d supported\n"
        sys.exit(-1)



    # maxdx is the maximum distance moved by a particle in a single
    # timestep
    maxdx = numpy.zeros((len(particles)), dtype=numpy.float64)
    maxdx[:] = -1.e33


    # dump out the data, particle-by-particle
    n = 0
    while (n < len(particles)):

        # get numpy arrays containing the time and coordinate
        # information for particle 0
        coords, time = particles[n].flatten()

        i = 1
        while (i < len(particles[n].history)):
            
            xo = coords[0,i-1]
            xn = coords[0,i]

            yo = coords[1,i-1]
            yn = coords[1,i]

            # compute distance -- check for periodic BCs
            if (xo > xmax-dx and xn < xmin+dx):
                # we flowed through right BC
                d = math.sqrt((xn - (xo - xmax))**2 + (yn - yo)**2)
                
            elif (xo < xmin+dx and xn > xmax-dx):
                # we flowed through the left BC
                d = math.sqrt((xn - (xo + xmax))**2 + (yn - yo)**2)

            else:
                # compute the distance
                d = math.sqrt((xn - xo)**2 + (yn - yo)**2)
            
            # store maximum
            maxdx[n] = max(maxdx[n],d)

            i += 1

        n += 1



    # output to a flie
    of = open("particle_maxdx.out", 'w')
    n = 0
    of.write("# particle,   max{d},     max{d}/dx\n")
    while (n < len(particles)):
        of.write("%d %f %f\n" % (n, maxdx[n], maxdx[n]/dx))
        n += 1
        
    of.close()



#-----------------------------------------------------------------------------
if __name__== "__main__":

    if (len(sys.argv) == 1):
        print "ERROR: no particle data files specified\n"
        sys.exit(2)

    main(sys.argv[1:])




        
