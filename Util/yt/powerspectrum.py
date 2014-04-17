from yt.mods import *
import numpy as na
import pylab


def powerspectrum(pf, nindex_rho=0.0):
    """Compute a density-weighted velocity power spectrum.  In particular,
    this computes:

             1           1  ^      ^*
     E(k) =  - integral  -  V(k) . V(k) dS
             N           2 

                 n
    where V = rho  U is the density-weighted velocity field, and N is the
    volume.  
      
    Parameters
    ----------
    pf : yt plotfile object

    nindex_rho: (optional) power to use for the density-weighting

    Returns
    -------
    k : ndarray
       The physical wavenumbers of the binned power spectrum

    E_spectrum : ndarray
       The power, E(k), at a given wavenumber
    """

    irho = ('gas', 'density')
    iu = ('gas', 'velocity_x')
    iv = ('gas', 'velocity_y')
    iw = ('gas', 'velocity_z')

    # get the data on a uniform grid at the highest resolution
    max_level = pf.index.max_level

    cube = pf.covering_grid(max_level, left_edge=pf.domain_left_edge, 
                            dims=pf.domain_dimensions,
                            fields=[irho, iu, iv, iw])


    # only triply-periodic is supported at the moment, although it is 
    # not at all clear that BCs should come into play for a sampled
    # covering grid
    pf.periodicity = (True, True, True)

    rho = cube[irho].d
    u = cube[iu].d
    v = cube[iv].d
    w = cube[iw].d
    
    nx, ny, nz = rho.shape

    # do the FFTs -- note that since our data is real, there will be
    # too much information here.  By default, fftn will put the
    # positive frequency terms in the first half of all axes -- that's
    # what we want to keep
    ruk = na.fft.fftn(rho**nindex_rho * u)[0:nx/2+1,0:ny/2+1,0:nz/2+1]
    rvk = na.fft.fftn(rho**nindex_rho * v)[0:nx/2+1,0:ny/2+1,0:nz/2+1]
    rwk = na.fft.fftn(rho**nindex_rho * w)[0:nx/2+1,0:ny/2+1,0:nz/2+1]

    # normalize -- the 8 here is because we are a real-valued function,
    # so only one octant is unique.
    ruk = 8.0*ruk/(nx*ny*nz)
    rvk = 8.0*rvk/(nx*ny*nz)
    rwk = 8.0*rwk/(nx*ny*nz)

    # wavenumbers -- unfortunately, yt uses an older version of NumPy,
    # so we don't have access to the rfftfreq function.  The last
    # element is negative, because of the symmetry, but should be
    # positive (rfftfreq would get this right).
    L = (pf.domain_right_edge - pf.domain_left_edge).d

    kx = na.fft.fftfreq(nx)[0:nx/2+1]
    kx[-1] *= -1
    kx = kx*nx/L[0]

    ky = na.fft.fftfreq(ny)[0:ny/2+1]
    ky[-1] *= -1
    ky = ky*ny/L[1]

    kz = na.fft.fftfreq(rho.shape[2])[0:nz/2+1]
    kz[-1] *= -1
    kz = kz*nz/L[2]

    # density-weighted kinetic energy density
    Kk = 0.5*(abs(ruk)**2 + abs(rvk)**2 + abs(rwk)**2)

    kmin = 0
    kmax = na.sqrt(na.max(kx)**2 + na.max(ky)**2 + na.max(kz)**2)

    N = na.floor(na.sqrt(nx**2 + ny**2 + nz**2))

    # bins holds the edges
    bins = na.linspace(kmin, kmax, N+1)

    kx3d, ky3d, kz3d = na.meshgrid(kx, ky, kz, indexing="ij")

    k = numpy.sqrt(kx3d**2 + ky3d**2 + kz3d**2)

    # bin the radii -- digitize returns an array with the same shape
    # as the input array but with elements of the array specifying
    # which bin that location belongs to.  The value of whichbin will
    # be 1 if we are located in the bin defined by bins[0] to bins[1].
    # This means that there will be no 0s
    whichbin = na.digitize(k.flat, bins)

    # bincount counts the number of occurrences of each non-negative
    # integer value in whichbin.  Each entry in ncount gives the number of
    # occurrences of it in whichbin.  The length of ncount is set by the
    # maximum value in whichbin
    ncount = na.bincount(whichbin)

    E_spectrum = na.zeros(len(ncount)-1, dtype=na.float64)

    n = 1
    while n < len(ncount):
        E_spectrum[n-1] = na.sum(Kk.flat[whichbin==n])/ncount[n]
        n += 1

    k = bins[1:]
    E_spectrum = E_spectrum[0:len(k)]

    return k, E_spectrum

