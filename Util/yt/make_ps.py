from yt.mods import *
import pylab
import powerspectrum
import numpy as na

pf = load("wd_576_base_plt216711")

k, Ek = powerspectrum.powerspectrum(pf, nindex_rho=1./3.)

n = 0
while n < len(k):
    print n, k[n], Ek[n]
    n += 1


index = na.argmax(Ek)

print 'index = ', index

kmax = k[index]
Emax = Ek[index]

print 'kmax = ', kmax
print 'Emax = ', Emax


pylab.loglog(k, Ek)
pylab.loglog(k, Emax*(k/kmax)**(-5./3.), ls=":", color="0.5")

pylab.savefig("spectrum.png")
