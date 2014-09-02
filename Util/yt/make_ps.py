import yt
import pylab
import powerspectrum
import numpy as np

ds = yt.load("plt164582")

zmin = 1300 
zmax = 3550

k, Ek = powerspectrum.powerspectrum(ds, nindex_rho=1./3., zrange=(zmin, zmax))

n = 0
while n < len(k):
    print n, k[n], Ek[n]
    n += 1


index = np.argmax(Ek)

print 'index = ', index

kmax = k[index]
Emax = Ek[index]

print 'kmax = ', kmax
print 'Emax = ', Emax


pylab.loglog(k, Ek)
pylab.loglog(k, Emax*(k/kmax)**(-5./3.), ls=":", color="0.5")

pylab.xlabel(r"$k$")
pylab.ylabel(r"$E(k)$")

pylab.savefig("spectrum.eps")
