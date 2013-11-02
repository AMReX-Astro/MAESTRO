import matplotlib.pyplot as plt
import numpy
import glob

labelDict={'max{enuc}': r"$\epsilon_{\rm max}$ (erg g$^{-1}$ s$^{-1}$)",
           'max{T}': r"$T_{\rm max}$ (K)",
           'max{Machno}': r"$M_{\rm max}$"}

fns = glob.glob("xrb*diag.out")
names = [fn.split('_')[1] for fn in fns]

for ifn,fn in enumerate(fns):
    print "working on %s" % fn
    with open(fn) as fh:
        for i,line in enumerate(fh):
            if i == 5: # skip the header info and get the labels
                cols = line.split()[1:]
                indx = [i for i in range(len(cols)) if cols[i].find("loc") < 0]
                cols = [c for c in cols if c.find("loc") < 0]

    print "Found these columns: %s" % ' '.join(cols)

    assert len(indx) == len(cols), \
        "indx and cols don't match: %s\n%s" % (' '.join(indx), ' '.join(cols))
    assert "time" in cols, "Couldn't find a time column!"
    nplts = len(indx)
    dat = numpy.loadtxt(fn,usecols=indx)
    for iplt in range(1,nplts):
        # only plot the quantities in labelDict
        if labelDict.has_key(cols[iplt]):
            plt.plot(dat[:,0],dat[:,iplt])
            ax = plt.gca()
            ax.set_xlabel("%s (s)" % cols[0])
            ax.set_ylabel(labelDict[cols[iplt]])
            fsave = "%s_%s" % (names[ifn],
                               cols[iplt].replace('{','_').strip('}'))
            plt.savefig("%s.eps" % fsave,bbox_inches='tight')
            plt.savefig("%s.pdf" % fsave,bbox_inches='tight')
            plt.clf()
            plt.cla()

