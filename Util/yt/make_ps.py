from yt.mods import *
import pylab
import powerspectrum
import numpy as na

pf = load("plt23437")

k, Ek = powerspectrum.powerspectrum(pf, nindex_rho=1./3.)

index = Ek == na.max(Ek)

kmax = k[index]
Emax = Ek[index]

pylab.loglog(k, Ek)
pylab.loglog(k, Emax*(k/kmax)**(-5./3.), ls=":", color="0.5")

pylab.savefig("spectrum.png")
