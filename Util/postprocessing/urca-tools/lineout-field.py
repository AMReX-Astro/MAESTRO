#!/usr/bin/env python
import matplotlib.pyplot as plt
import yt
from yt import derived_field
import numpy as np
from yt_urca_fields import UrcaShellFields, DatasetHelpers
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('infile', type=str, help='Name of input plotfile.')
parser.add_argument('-f', '--field', type=str, default='density',
                    help='Name of the field to plot. Eg. "tfromp". Default is density.')
parser.add_argument('-axis', '--axis', type=str, default='x',
                    help='Axis along which to take the lineout. Default is "x". Can be "x", "y", or "z".')
parser.add_argument('-w', '--width', type=float,
                    help='Width of lineout (cm). Default is the full width of the domain.')
parser.add_argument('-flo', '--flo', type=float,
                    help='Lower bound of field in lineout plot.')
parser.add_argument('-fhi', '--fhi', type=float,
                    help='Upper bound of field in lineout plot.')
parser.add_argument('-log', '--logscale', action='store_true', help='If supplied, use a log scale for the field.')
args = parser.parse_args()


if __name__ == "__main__":
    ds = yt.load(args.infile)

    ushell_fields = UrcaShellFields()
    ushell_fields.setup(ds)
    
    field, field_short_name = DatasetHelpers.get_field(ds, args.field)
    assert(field)

    c = ds.domain_center

    axmap = {'x': 0, 'y': 1, 'z': 2}
    axis_str = args.axis.lower()
    ax = axmap[axis_str]

    transverse_indices = [0,1,2]
    transverse_indices.pop(ax)

    # cut through the transverse axis such that the ray intersects the center of the domain
    ray = ds.ortho_ray(ax, (c[transverse_indices[0]], c[transverse_indices[1]]))

    # Sort the ray values by axis coordinate so there are no discontinuities
    srt = np.argsort(ray[axis_str])

    plt.subplot(111)
    if args.logscale:
        plt.semilogy(np.array(ray[axis_str][srt]), np.array(ray[field][srt]))
    else:
        plt.plot(np.array(ray[axis_str][srt]), np.array(ray[field][srt]))
    if args.width:
        center_axis = c[ax].in_units('cm').d
        lower = center_axis - 0.5*args.width
        upper = center_axis + 0.5*args.width
        plt.gca().set_xlim(left=lower, right=upper)
    plt.gca().set_ylim(bottom=args.flo, top=args.fhi)
    plt.xlabel(axis_str)
    plt.ylabel(field_short_name)
    plotname = "{}.lineout.{}.{}.png".format(args.infile, axis_str, field_short_name)
    print('Saving lineout plot: {}'.format(plotname))
    plt.savefig(plotname)
