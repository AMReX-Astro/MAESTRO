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

    ref = int(na.product(pf.ref_factors[0:max_level]))

    # allocate for our uniformly-gridded result
    dims = pf.domain_dimensions*ref
    nx, ny, nz = dims

    Kk = na.zeros( (nx/2+1, ny/2+1, nz/2+1), dtype=na.float32)

    print "doing ux"
    Kk += 0.5*fft_comp(pf, irho, iu, nindex_rho, max_level, dims)

    print "doing uy"
    Kk += 0.5*fft_comp(pf, irho, iv, nindex_rho, max_level, dims)

    print "doing uz"
    Kk += 0.5*fft_comp(pf, irho, iw, nindex_rho, max_level, dims)

    # wavenumbers
    L = (pf.domain_right_edge - pf.domain_left_edge).d

    kx = na.fft.fftfreq(nx)[0:nx/2+1]
    kx[-1] *= -1
    kx = kx*nx/L[0]

    ky = na.fft.fftfreq(ny)[0:ny/2+1]
    ky[-1] *= -1
    ky = ky*ny/L[1]

    kz = na.fft.fftfreq(nz)[0:nz/2+1]
    kz[-1] *= -1
    kz = kz*nz/L[2]

    kmin = 0
    kmax = na.sqrt(na.max(kx)**2 + na.max(ky)**2 + na.max(kz)**2)

    N = na.floor(na.sqrt(nx**2 + ny**2 + nz**2))

    # bins holds the edges
    bins = na.linspace(kmin, kmax, N+1)

    kx3d, ky3d, kz3d = na.meshgrid(kx, ky, kz, indexing="ij", dtype=na.float32)

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

    E_spectrum = na.zeros(len(ncount)-1, dtype=na.float32)

    n = 1
    while n < len(ncount):
        if ncount[n] == 0: break
        E_spectrum[n-1] = na.sum(Kk.flat[whichbin==n]) #/ncount[n]
        n += 1

    k = bins[1:n]
    E_spectrum = E_spectrum[0:len(k)]

    return k, E_spectrum



def fft_comp(pf, irho, iu, nindex_rho, level, dims ):

    print "covering grid"
    cube = pf.covering_grid(level, left_edge=pf.domain_left_edge, 
                            dims=dims,
                            fields=[irho, iu])

    print "accessing rho cube"
    rho = cube[irho].d    

    print rho.dtype
    nx, ny, nz = rho.shape

    # do the FFTs -- note that since our data is real, there will be
    # too much information here.  By default, fftn will put the
    # positive frequency terms in the first half of all axes -- that's
    # what we want to keep

    # normalize -- the 8 here is because we are a real-valued function,
    # so only one octant is unique.

    print "accessing cube"
    u = cube[iu].d

    print "ffting"
    ru = na.fft.fftn(rho**nindex_rho * u)[0:nx/2+1,0:ny/2+1,0:nz/2+1]

    print "normalizing"
    ru = 8.0*ru/(nx*ny*nz)

    print "returning"
    return abs(ru)**2
