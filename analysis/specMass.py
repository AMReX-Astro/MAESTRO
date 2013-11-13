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
import os.path

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

specLabels={"C12": r"$^{12}\mathrm{C}$",
            "O14": r"$^{14}\mathrm{O}$",
            "O15": r"$^{15}\mathrm{O}$",
            "O16": r"$^{16}\mathrm{O}$",
            "F17": r"$^{17}\mathrm{F}$",
            "Mg22": r"$^{22}\mathrm{Mg}$",
            "S30": r"$^{30}\mathrm{S}$",
            "Ni56": r"$^{56}\mathrm{Ni}$",
            "He4": r"$^{4}\mathrm{He}$",
            "H1": r"$^{1}\mathrm{H}$"}

justPlot = False
append = False
every = 1

# parse the cmdline options
justPlot = '--just-plot' in sys.argv
append = '--append' in sys.argv
everyFlag = "--every" in sys.argv
if everyFlag:
    every = int(sys.argv[sys.argv.index("--every")+1])
    print "Only summing every %s files" % every

if justPlot:
    print "Just doing the plots"
    masses = numpy.loadtxt(datFile)
    fh = open(datFile,'r')
    labels = fh.readline().split()[1:]
    # make sure this is the label list
    if "time" not in labels:
        labels = fh.readline().split()[1:]
    fh.close()
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
    lastFile = "last file used: %s\n" % fns[-1]
    header = [lastFile,] + labels

    fh = StringIO.StringIO(mass)
    masses = numpy.loadtxt(fh)
    if append:
        print "Just appending the new values to the old %s file" % datFile
        oldMasses = numpy.loadtxt(datFile)
        masses = numpy.vstack((oldMasses,masses))
    if os.path.exists(datFile):
        raw_input("Getting ready to overwrite %s; " % datFile
                  + "continue if this is ok, otherwise C-c and quit.")
    numpy.savetxt(datFile,masses,header=' '.join(header))

# normalize to starting value
print 'Normalizing...'
startVals = masses[0,1:].copy()
masses[:,1:] /= startVals

print "Plotting..."
for i,curve in enumerate(labels[1:]):
    plt.plot(masses[:,0],masses[:,i+1],lts[i],label=specLabels[curve])
# thin line at 1
plt.axhline(y=1.0,lw=1.0,ls=':',color='0.5')

ax = plt.gca()
ax.semilogy()
ax.set_xlabel(r"$t$ (s)")
ax.set_ylabel(r"$M_X/M_{X,0}$")
ax.set_xlim(0,0.175)
plt.legend(loc="upper center",ncol=4, frameon=False)

plt.savefig("masses.eps",bbox_inches='tight')
plt.savefig("masses.pdf",bbox_inches='tight')
