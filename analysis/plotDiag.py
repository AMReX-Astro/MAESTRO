import matplotlib.pyplot as plt
import numpy
import glob
import sys

labelDict={'max{enuc}': r"$\epsilon_{\rm max}$ (erg g$^{-1}$ s$^{-1}$)",
           'max{T}': r"$T_{\rm max}$ (K)",
           'max{Machno}': r"$M_{\rm max}$"}

smoothType = None 
smoothSize = 49     # needs to be odd

# parse the args
if "--smooth-type" in sys.argv:
    smoothType = sys.argv[sys.argv.index("--smooth-type")+1]
if "--smooth-size" in sys.argv:
    smoothSize = int(sys.argv[sys.argv.index("--smooth-size")+1])

if smoothType:
    print "Doing smoothing with a %s window of size %d" % (smoothType,
                                                           smoothSize)

# smoothing functions
# taken from wiki.scipy.org/Cookbook/SignalSmooth
def smooth(x,window_len=11,window=None):
    """smooth the data using a window with requested size.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd 
                    integer
        window: the type of window from 'flat', 'hanning', 'hamming', 
                'bartlett', 'blackman'
                flat window will produce a moving average smoothing.

    NOTE: length(output) != length(input), to correct this: 
          return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """
    if not window: return x
    if x.ndim != 1: raise ValueError, "smooth only accepts 1 dimension arrays."
    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."
    if window_len < 3: return x
    if window not in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is not one of 'flat', 'hanning', 'hamming', 'bartlett', or 'blackman'"
    s=numpy.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    if window == 'flat': # moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')
    # this has size len(s)-window_len-1 with the first and last window_len/2
    # elements added from padding and convolving
    y=numpy.convolve(w/w.sum(),s,mode='valid')
    # lets get the right size back; note Cookbook recipe has typo here
#    return y
    return y[(window_len/2):-(window_len/2)]


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
            # apply smoothing
            dat[:,iplt] = smooth(dat[:,iplt],smoothSize,smoothType)
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

