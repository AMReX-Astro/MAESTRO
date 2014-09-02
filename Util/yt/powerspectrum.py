import yt
import numpy as np
import pylab

def powerspectrum(ds, nindex_rho=0.0, zrange=None):
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
    ds : yt dataset object

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
    max_level = ds.index.max_level

    ref = int(np.product(ds.ref_factors[0:max_level]))

    # allocate for our uniformly-gridded result
    dims = ds.domain_dimensions*ref


    if not zrange == None:
        zmin, zmax = zrange

        # figure out the range in zones in z that we want to keep
        dz = ds.domain_width[2]/dims[2]
        z = (np.arange(dims[2])+0.5)*dz + ds.domain_left_edge[2]

        izlo = np.nonzero(zmin > z)[0]
        izlo = izlo[len(izlo)-1]
        
        izhi = np.nonzero(zmax > z)[0]
        izhi = izhi[len(izhi)-1]

        # make the range even -- for FFT sake
        if not (izhi-izlo) % 2 == 0:
            izhi += 1

        low = (ds.domain_left_edge[0], ds.domain_left_edge[1], izlo)
        delta = (dims[0], dims[1], izhi-izlo)

    else:

        low = ds.domain_left_edge
        delta = dims


    nx, ny, nz = delta

    Kk = np.zeros( (nx/2+1, ny/2+1, nz/2+1), dtype=np.float32)

    print "doing ux"
    Kk += 0.5*fft_comp(ds, irho, iu, nindex_rho, max_level, low, delta)

    print "doing uy"
    Kk += 0.5*fft_comp(ds, irho, iv, nindex_rho, max_level, low, delta)

    print "doing uz"
    Kk += 0.5*fft_comp(ds, irho, iw, nindex_rho, max_level, low, delta)

    # wavenumbers
    L = (ds.domain_right_edge - ds.domain_left_edge).d

    if not zrange == None:
        L[2] = z[izhi] - z[izlo]

    kx = np.fft.fftfreq(nx)[0:nx/2+1]
    kx[-1] *= -1
    kx = kx*nx/L[0]

    ky = np.fft.fftfreq(ny)[0:ny/2+1]
    ky[-1] *= -1
    ky = ky*ny/L[1]

    kz = np.fft.fftfreq(nz)[0:nz/2+1]
    kz[-1] *= -1
    kz = kz*nz/L[2]

    kmin = 0
    kmax = np.sqrt(np.max(kx)**2 + np.max(ky)**2 + np.max(kz)**2)

    N = np.floor(np.sqrt(nx**2 + ny**2 + nz**2))

    # bins holds the edges
    bins = np.linspace(kmin, kmax, N+1)

    kx3d, ky3d, kz3d = np.meshgrid(kx, ky, kz, indexing="ij", dtype=np.float32)

    k = np.sqrt(kx3d**2 + ky3d**2 + kz3d**2)

    # bin the radii -- digitize returns an array with the same shape
    # as the input array but with elements of the array specifying
    # which bin that location belongs to.  The value of whichbin will
    # be 1 if we are located in the bin defined by bins[0] to bins[1].
    # This means that there will be no 0s
    whichbin = np.digitize(k.flat, bins)

    # bincount counts the number of occurrences of each non-negative
    # integer value in whichbin.  Each entry in ncount gives the number of
    # occurrences of it in whichbin.  The length of ncount is set by the
    # maximum value in whichbin
    ncount = np.bincount(whichbin)

    E_spectrum = np.zeros(len(ncount)-1, dtype=np.float32)

    n = 1
    while n < len(ncount):
        if ncount[n] == 0: break
        E_spectrum[n-1] = np.sum(Kk.flat[whichbin==n]) 
        n += 1

    k = bins[1:n]
    E_spectrum = E_spectrum[0:len(k)]

    return k, E_spectrum



def fft_comp(ds, irho, iu, nindex_rho, level, low, delta ):

    print "covering grid"
    cube = ds.covering_grid(level, left_edge=low,
                            dims=delta,
                            fields=[irho, iu])

    print "accessing rho cube"
    rho = cube[irho].d    

    print rho.shape
    nx, ny, nz = rho.shape

    # do the FFTs -- note that since our data is real, there will be
    # too much information here.  By default, fftn will put the
    # positive frequency terms in the first half of all axes -- that's
    # what we want to keep

    # normalize -- the 8 here is because we are a real-valued function,
    # so only one octant is unique.

    u = cube[iu].d

    ru = np.fft.fftn(rho**nindex_rho * u)[0:nx/2+1,0:ny/2+1,0:nz/2+1]

    ru = 8.0*ru/(nx*ny*nz)

    return abs(ru)**2
