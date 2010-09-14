# a simple script to plot the initial density and location of cutoff density and
# sponge parameters

anelastic_cutoff = 100.0
sponge_center_density = 100.0
sponge_start_factor = 2.5



import math
import numpy
import pylab
import string

def initial_model():

    initialModel = "glasner_T5_Ginvsq_Snone.hse"

    # read in the initial model and store the data in 
    # initialData
    mf = open(initialModel, "r")
    
    numLines = 0
    numFields = -1
    for line in mf:
        if (not line.startswith("#")):
            numLines +=1
            
            if (numFields == -1):
                numFields = len(string.split(line))

    mf.close()

    initialData = numpy.zeros( (numLines, numFields), numpy.float64)

    mf = open(initialModel, "r")
    
    indexLine = 0
    for line in mf:
        if (not line.startswith("#")):
            initialData[indexLine,:] = string.split(line)
            indexLine += 1
            
    mf.close()

    # find the r coordinate of the anelastic cutoff
    n = 0
    while (n < numLines):
        if (initialData[n,1] <= anelastic_cutoff):
            cutoffR = initialData[n,0]
            break

        n += 1


    # find the r coordinate of the middle of the sponge
    n = 0
    while (n < numLines):
        if (initialData[n,1] <= sponge_center_density):
            r_md = initialData[n,0]
            break

        n += 1


    # find the r coordinate of the start of the sponge
    sponge_start = sponge_start_factor*sponge_center_density
    n = 0
    while (n < numLines):
        if (initialData[n,1] <= sponge_start):
            r_sp = initialData[n,0]
            break

        n += 1
        
    spongeR = r_sp

    # compute the r coordinate of the top of the sponge
    r_tp = 2*r_md - r_sp


    print "cutoffR = ", cutoffR
    print "r_sp = ", r_sp
    print "r_md = ", r_md
    print "r_tp = ", r_tp


    # make the plots
    rMin = 0
    rMax = 3.e8

    # density
    densMin = 1.e-4
    densMax = 1.e10
        
    sp = pylab.subplot(211)

    sp.set_yscale('log')

    pylab.plot(initialData[:,0], initialData[:,1], color="k")

    ax = pylab.gca()
    ax.xaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))
    #ax.yaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))

    # draw in the anelastic cutoff and sponge
    pylab.plot([cutoffR, cutoffR], [densMin, densMax], color="0.5", linestyle=":")
    pylab.plot([spongeR, spongeR], [densMin, densMax], color="0.5", linestyle="--")
    pylab.plot([r_md, r_md ], [densMin, densMax], color="0.5", linestyle="--")
    pylab.plot([r_tp, r_tp ], [densMin, densMax], color="0.5", linestyle="--")

    pylab.xlabel("r (cm)")
    pylab.ylabel(r"density (g cm$^{-3}$)")

    pylab.xlim(4.1e8,5.1e8)
    pylab.ylim(1.e-3,1.e7)

    # sponge plot
    sp = pylab.subplot(212)

    r = initialData[:,0].copy()
    f_damp = numpy.zeros( (len(r) ), numpy.float64)

    n = 0
    while (n < len(r)):

        if (r[n] < r_sp):
            f_damp[n] = 0
        elif (r[n] >= r_sp and r[n] < r_tp):
            f_damp[n] = 0.5*(1 - math.cos(math.pi*( (r[n] - r_sp)/(r_tp - r_sp) )) )
        else:
            f_damp[n] = 1.0

        n += 1





    pylab.plot(r, f_damp, color='k')
    pylab.xlabel("r (cm)")
    pylab.ylabel(r"$f_\mathrm{damp}$")

    ax = pylab.gca()
    ax.xaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))
    #ax.yaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))

    pylab.xlim(4.1e8,5.1e8)
    pylab.ylim(-0.1,1.1)



    pylab.subplots_adjust(hspace=0.3)

    f = pylab.gcf()
    f.set_size_inches(6.0,10.0)

    pylab.savefig("initial_model_paper.png")



if __name__== "__main__":
    initial_model()

