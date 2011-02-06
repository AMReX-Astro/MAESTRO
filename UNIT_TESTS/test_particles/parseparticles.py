#!/usr/bin/env python

import sys
import os
import string
import math
import numpy
import pylab

#-----------------------------------------------------------------------------
# a particleInstance object represents a single instance of a
# particle, and will carry the position, time, and associated data.
class particleInstance:

    def __init__(self, x, y, z, t):
        self.x = x
        self.y = y
        self.z = z
        self.t = t
        self.data = {}


    def addData(self, name, value):
        self.data[name] = value


    def __str__(self):
        string = "particle pos:  (%g, %g, %g) \n" % (self.x, self.y, self.z) + \
                 "         time: %g \n" % (self.t) + \
                 "         data: \n"
        for key in self.data.keys():
            string += "            %s: %g \n" % (key, self.data[key])
        
        return string

    
    # we will sort based on the particle instance time
    def value(self):
        return self.t

    def __cmp__(self, other):
        return cmp(self.value(), other.value())



#-----------------------------------------------------------------------------
# a particle object stores a history of a single particle, each
# element of which is a particleInstance object.
class particle:
    
    def __init__(self, pid, originCPU, dim):
        
        # a MAESTRO particle is identified by 2 numbers, the pid and
        # the CPU that it was created on.  Together, these uniquely
        # identify the particle.
        self.pid = pid
        self.originCPU = originCPU
        self.dim = dim

        # the history list will store instances of the particle at
        # different times.
        self.history = []

        # keep track of the number of particle history instances we've
        # stored
        self.numInstances = 0
    

    def addInstance(self, x=0.0, y=0.0, z=0.0, t=0.0, 
                    dataNames=[], dataValues=[]):
        
        # add this particle instance to the particle history
        self.history.append(particleInstance(x,y,z,t))

        # sanity check
        if (not (len(dataNames) == len(dataValues)) ):
            print "ERROR: len(names) != len(values)"
            sys.exit(2)
            
        # add the data associate with the particle to the particle instance
        n = 0
        while (n < len(dataNames)):
            self.history[self.numInstances].addData(dataNames[n], dataValues[n])
            n += 1

        self.numInstances += 1


    # return numpy arrays containing the particle data
    def flatten(self):

        # sort the history by time -- the particleInstance objects use
        # time as their value for comparison, so sorted() should handle
        # this nicely for us.
        shistory = sorted(self.history)

        coords = numpy.zeros((self.dim, self.numInstances), dtype=numpy.float64)
        time   = numpy.zeros((self.numInstances), dtype=numpy.float64)

        n = 0
        while (n < len(shistory)):
            coords[0,n] = shistory[n].x
            if (self.dim >= 2):
                coords[1,n] = shistory[n].y
            if (self.dim == 3): 
                coords[2,n] = shistory[n].z
            
            time[n] = shistory[n].t

            n += 1

        return coords, time


    def __str__(self):
        string = "particle ID:           %d \n" % (self.pid) + \
                 "particle origin CPU:   %d \n" % (self.originCPU) + \
                 "dimensionality:        %d \n" % (self.dim) + \
                 "num. stored instances: %d \n" % (self.numInstances)

        return string



#-----------------------------------------------------------------------------
# loop through all of the elements of the particleList list, looking
# for the particle object that has the desired pid and originCPU.  If
# none is found, return -1.
def findParticleIndex(particleList, pid, originCPU):
    index = -1

    n = 0
    while (n < len(particleList)):
        if (particleList[n].pid == pid and particleList[n].originCPU == originCPU):
            index = n
            exit
            
        n += 1

    return index



#-----------------------------------------------------------------------------
# read in all the particle data from a file and add each particle instance
# to the particleList list, grouping all of the history for a unique particle
# (identified by the pid and originCPU) together in a particle object in
# the list.
def parseParticleFile(MaestroParticleFile, particleList):

    # read the file line by line
    mf = open(MaestroParticleFile, "r")

    for line in mf:

        # skip blank lines
        if (line.lstrip() == ""):
            continue

        # look for a header information
        if (line.startswith("#")):
            fields = string.split(line[1:])

            # make sure we know what we are doing -- the first 2
            # fields should be the particle ID and origin CPU
            if (fields[0] == "part-ID" and fields[1] == "origin-CPU"):
                ipid = 0
                ioCPU = 1

            else:
                print "ERROR: particle file columns not in expected order"
                sys.exit(2)


            # the next fields should be x, y, and z, depending on the
            # dimensionality
            if (fields[2] == "x" and fields[3] == "y" and fields[4] == "z"):
                dim = 3
                ix = 2; iy = 3; iz = 4
                
            elif (fields[2] == "x" and fields[3] == "y"):
                dim = 2
                ix = 2; iy = 3

            elif (fields[2] == "x"):
                dim = 1
                ix = 2

            else:
                print "ERROR: particle file columns not in expected order"
                sys.exit(2)


            # then comes time
            if (fields[2 + dim] == "time"):
                it = 2 + dim
                
            else:
                print "ERROR: particle file columns not in expected order"
                sys.exit(2)
            

            # everything else is associated data
            if (len(fields) > 3 + dim):
                idata = 3 + dim
                ndata = len(fields) - idata
                dataNames = fields[idata:]
            else:
                ndata = 0

            # done with the header
            continue


        # normal particle data -- first find the required particle data
        fields = string.split(line)

        pid = int(fields[ipid])
        originCPU = int(fields[ioCPU])

        if (dim == 1):
            x = float(fields[ix])
            y = -1
            z = -1

        elif (dim == 2):
            x = float(fields[ix])
            y = float(fields[iy])
            z = -1

        else:
            x = float(fields[ix])
            y = float(fields[iy])
            z = float(fields[iz])


        time = float(fields[it])

        dataValues = []
        values = fields[idata:]
        n = 0
        while (n < len(values)):
            dataValues.append(float(values[n]))
            n += 1

        # add the particle instance to the right particle object in
        # the particleList list
        index = findParticleIndex(particleList, pid, originCPU)

        if (index == -1):
            particleList.append(particle(pid, originCPU, dim))
            index = findParticleIndex(particleList, pid, originCPU)

        particleList[index].addInstance(x=x, y=y, z=z, t=time, 
                                        dataNames=dataNames, 
                                        dataValues=dataValues)
    


#-----------------------------------------------------------------------------
def main(files):

    # domain information
    xcenter = 0.5
    ycenter = 0.5


    # the master list of particles -- each element will store a particle
    # object containing the complete history of that particle
    particles = []

    # loop over the files
    for f in files:

        # make sure that the file exists
        if (not os.path.isfile(f)):
            print "ERROR: file %s not found\n" % (f)
    
        # parse the file, creating particle objects for each unique particle
        # found and storing the history
        parseParticleFile(f, particles)


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
        r_start = math.sqrt( (particles[n].history[0].x - xcenter)**2 +
                             (particles[n].history[0].y - ycenter)**2)

        r_end = math.sqrt( (particles[n].history[nstep-1].x - xcenter)**2 +
                           (particles[n].history[nstep-1].y - ycenter)**2)

        error[n] = math.fabs(r_start - r_end)/r_start

        print "particle: (%d, %d), r init = %f, rel error = %g" % \
            (particles[n].pid, particles[n].originCPU, r_start, error[n])
        
        n += 1



#-----------------------------------------------------------------------------
if __name__== "__main__":

    if (len(sys.argv) == 1):
        print "ERROR: no particle data files specified\n"
        sys.exit(2)

    main(sys.argv[1:])




        
