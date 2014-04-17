from yt.mods import *
import pylab
import powerspectrum

pf = load("plt23437")

k, Ek = powerspectrum.powerspectrum(pf, nindex_rho=1./3.)

pylab.loglog(k, Ek)

pylab.savefig("spectrum.png")
