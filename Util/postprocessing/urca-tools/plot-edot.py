#!/usr/bin/env python
"""
Use yt to plot the integrated energy generation rate over
concentric spheres of varying radius from the center
of the star outward.

Energy generation rate is SUM(enucdot * density * cell volume)
and has units of erg/s

Donald E. Willcox
"""
import yt
from yt import derived_field
import numpy as np
import matplotlib

matplotlib.use('Agg')

from matplotlib import pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('infile', type=str, help='Name of input plotfile.')
parser.add_argument('-maxr', '--max_radius', type=float, help='Maximum radius (km).')
parser.add_argument('-minr', '--min_radius', type=float, default=10,
                    help='Minimum radius (km). Defaults to 10 km.')
parser.add_argument('-n', '--num_radial_steps', type=int, default=1000,
                    help='Number of radial steps to take between --max_radius and --min_radius. Default is 1000.')
args = parser.parse_args()

@derived_field(name="edot", units="erg/s")
def _edot(field, data):
    return data[('boxlib', 'enucdot')] * data[('boxlib', 'density')] * data[('gas', 'cell_volume')]

def get_max_radius(ds):
    # Given dataset ds, get the maximum radius of an inscribing sphere
    maxradius = min(0.5*(ds.domain_width)).in_units('km')
    return maxradius

def get_rads_edots(ds):
    rads = np.linspace(args.min_radius, get_max_radius(ds).to_ndarray(), num=args.num_radial_steps)
    edots = []
    for r in rads:
        s = ds.sphere(ds.domain_center, (r, 'km'))
        edots.append(s.sum('edot'))
    return rads, edots

dsname = args.infile
ds = yt.load(dsname)

rads, edots = get_rads_edots(ds)

plt.plot(rads, edots, color='black')
plt.xlabel('radius (km)')
plt.ylabel('Integrated $\\dot{E}(r)$ (erg/s)')
if args.max_radius:
    plt.xlim([0, args.max_radius])
plt.tight_layout()
plt.savefig('{}_edotvr.png'.format(dsname), dpi=600)
