#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('modelfile', type=str, help='Name of MAESTRO initial model file.')
args = parser.parse_args()

Msun = 1.988435e33 # grams

def getlinedata(f):
    l = f.readline()
    ls = l.strip().split()
    a = [float(x) for x in ls]
    na = len(a)
    return na, a

def plotvsk(data, key, xlab):
    # key should be the key into data to use
    # as the independent variable in the plot.
    # xlab is the string for the xlabel.
    # Plot each variable versus data[key]
    for c in colnames:
        d = data[c]
        m = data[key]
        plt.subplots()
        plt.plot(m,d,'bo')
        plt.xlabel(xlab)
        plt.ylim(np.amin(d), np.amax(d))
        plt.ylabel(c)
        plt.tight_layout()
        plt.savefig(c.replace(' ','_')+'_{}.eps'.format(key))
        plt.clf()

if __name__=="__main__":
    f = open(args.modelfile, 'r')
    nstr = f.readline().strip()
    npts  = int(nstr.split()[-1])
    nvarstr = f.readline().strip()
    nvar = int(nvarstr.split()[-1])

    colnames = ['radius']

    for i in range(nvar):
        cstr = f.readline().strip()
        cname = cstr[2:]
        colnames.append(cname)

    ncols = len(colnames)

    data = {}
    for c in colnames:
        data[c] = []
    for i in range(npts):
        nd, d = getlinedata(f)
        while nd < ncols:
            nd2, d2 = getlinedata(f)
            nd = nd + nd2
            for x in d2:
                d.append(x)
        for c, x in zip(colnames, d):
            data[c].append(x)

    for c in colnames:
        data[c] = np.array(data[c])

    print('Found columns:')
    for k in data.keys():
        print(k)

    # Add mass coordinate to dataset
    dx = data['radius'][1]-data['radius'][0]
    data['zone_xhi'] = data['radius']+0.5*dx
    data['zone_xlo'] = data['radius']-0.5*dx
    data['zone_mass'] = data['density']*(data['zone_xhi']**3 - data['zone_xlo']**3)*4.0*np.pi/3.0
    data['mass'] = np.cumsum(data['zone_mass'])/Msun

    # Plot each variable versus the mass coordinate
    plotvsk(data, 'mass', 'Mass ($M_\\odot$)')

    # Plot each variable versus the radius coordinate
    plotvsk(data, 'radius', 'radius ($cm$)')

    # Plot each variable versus the density coordinate
    plotvsk(data, 'density', 'density ($g/cm^3$)')

    # Plot each variable versus the pressure coordinate
    plotvsk(data, 'pressure', 'pressure ($erg/cm^3$)')


