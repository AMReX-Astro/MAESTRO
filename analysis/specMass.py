import numpy
import matplotlib
matplotlib.rcParams['font.size'] =  18
matplotlib.rcParams["lines.linewidth"] = 2
matplotlib.rcParams["legend.columnspacing"] = 0.5
matplotlib.rcParams["legend.labelspacing"] = 0.1
matplotlib.rcParams["legend.fontsize"] = "medium"
matplotlib.rcParams["legend.frameon"] = False
import matplotlib.pyplot as plt
import subprocess
import glob
import StringIO
import sys

# comparison function for sorting plt's with 5 vs 6 digits
def fnCmp(fn1,fn2):
    # strip the 'plt' off and compare
    fn1=int(fn1[3:])
    fn2=int(fn2[3:])
    return cmp(fn1,fn2)

datFile = "masses.txt"
fmass = "./fspec_total_mass.Linux.gfortran.exe"

lts = ['b-','g-','r-','c-','m-',
       'b--','g--','r--','c--','m--']

justPlot = False
append = False
every = 1

justPlot = '--just-plot' in sys.argv
append = '--append' in sys.argv
everyFlag = "--every" in sys.argv
if everyFlag:
    every = int(sys.argv[sys.argv.index("--every")+1])
    print "Only summing every %s files" % every

if justPlot:
    print "Just doing the plots"
    masses = numpy.loadtxt(datFile)
    labels = open(datFile,'r').readline().split()[1:]
    # make sure this is the label list; might be second header line
    if "time" not in labels:
        labels = open(datFile,'r').readline().split()[1:]
    assert "time" in labels, "Can't find headers!"
else:
    fns = glob.glob("plt*")
    fns.sort(cmp=fnCmp)
    fns = fns[::every]
    print "Working these files:\n %s" % ' '.join(fns)

    cmd = "%s %s" % (fmass, ' '.join(fns))

    mass = subprocess.Popen(cmd,shell=True,
                            stdout=subprocess.PIPE).communicate()[0]
    labels = mass.split('\n')[0].split()[1:]
    lastFile = "# last file used: %s\n" fns[-1]
    labels = [lastFile,] + labels

    fh = StringIO.StringIO(mass)
    masses = numpy.loadtxt(fh)
    if append:
        print "Just appending the new values to the old %s file" % datFile
        oldMasses = numpy.loadtxt(datFile)
        masses = numpy.vstack((oldMasses,masses))
    numpy.savetxt(datFile,masses,header=' '.join(labels))

# normalize to starting value
print 'Normalizing...'
masses[:,1:] /= masses[0,1:]

print "Plotting..."
for i,curve in enumerate(labels[1:]):
    plt.plot(masses[:,0],masses[:,i+1],lts[i],label=curve)

ax = plt.gca()
ax.semilogy()
ax.set_xlabel(r"$t$ (s)")
ax.set_ylabel(r"$M_X/M_{X,0}$")
plt.legend(loc="upper center",ncol=4)

plt.savefig("masses.eps",bbox_inches='tight')
plt.savefig("masses.pdf",bbox_inches='tight')
