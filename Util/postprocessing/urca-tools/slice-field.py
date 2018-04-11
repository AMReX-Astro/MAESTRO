#!/usr/bin/env python
"""
Use yt to slice a boxlib plotfile supplied through the domain center.

Donald E. Willcox
"""
import yt
from yt import derived_field
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('infile', type=str, help='Name of input plotfile.')
parser.add_argument('-f', '--field', type=str,
                    help='Name of the field to plot. Eg. "(boxlib, tfromp)". Default is to slice all fields.')
parser.add_argument('-axis', '--axis', type=str, default='x',
                    help='Axis across which to slice at the center of the domain. Default is "x".')
parser.add_argument('-w', '--width', type=float,
                    help='Width of slice (cm). Default is domain width.')
parser.add_argument('-log', '--logscale', action='store_true', help='If supplied, use a log scale for the field.')
parser.add_argument('-symlog', '--symlog', type=float, help='If supplied, use symlog scaling, which is linear near zero, to accomodate positive and negative values of the field. Pass the value of the field at which to linearize the colorbar.')
parser.add_argument('-ctr', '--center', type=float, nargs='+', help='Centers the plot on the coordinates provided (x, y, z).')
parser.add_argument('-min', '--field_min', type=float, help='Minimim field value for colormap.')
parser.add_argument('-max', '--field_max', type=float, help='Maximum field value for colormap.')
parser.add_argument('-cmap', '--colormap', type=str, default='viridis',
                    help='Name of colormap to use. Default is "viridis".')
parser.add_argument('-res', '--resolution', type=int, default=2048,
                    help='Resolution to use in each direction in pixels. Default is 2048.')
parser.add_argument('-dc', '--drawcells', action='store_true', help='If supplied, draw the cell edges.')
parser.add_argument('-dg', '--drawgrids', action='store_true', help='If supplied, draw the grids.')
parser.add_argument('-octant', '--octant', action='store_true', help='Sets slice view appropriately for octant dataset.')
parser.add_argument('-natorg', '--native_origin', action='store_true', help='Use the native origin location for the axes.')
parser.add_argument('-ls', '--list_fields', action='store_true', help='If supplied, do nothing except list the available fields.')
args = parser.parse_args()

def slicefield(ds, field, field_short_name):
    if not args.width:
        width = max(ds.domain_width)
    else:
        width = yt.YTQuantity(args.width, 'cm')

    center_loc = None

    if args.octant:
        if args.center and len(args.center) == 3:
            center_loc = ds.arr(args.center, 'cm')
        else:
            dcenter = width.in_units('cm').v/2.0
            center_loc = ds.arr([dcenter, dcenter, dcenter], 'cm')
        s = yt.SlicePlot(ds, args.axis, field, center=center_loc, width=width, origin="native")
    else:
        if args.center and len(args.center) == 3:
            center_loc = ds.arr(args.center, 'cm')
        else:
            center_loc = 'c'
        if args.native_origin:
            s = yt.SlicePlot(ds, args.axis, field, center=center_loc, width=width, origin="native")
        else:
            s = yt.SlicePlot(ds, args.axis, field, center=center_loc, width=width)

    # Colormaps and Scaling
    maxv = ds.all_data().max(field)
    minv = ds.all_data().min(field)
    if (minv < 0.0 and maxv > 0.0 and
        not (args.field_min and args.field_max)):
        # Use symlog scaling and a two-tone colormap
        pos_maxv = np.ceil(np.log10(maxv))
        neg_maxv = np.ceil(np.log10(minv))
        dlog = abs(np.log10(maxv)) + abs(np.log10(abs(minv)))
        logmaxv = max(pos_maxv, neg_maxv)
        linmaxv = max(maxv, -minv)
        s.set_cmap(field, 'PiYG')
        if args.logscale:
            if args.symlog:
                s.set_log(field, args.logscale, linthresh=args.symlog)
            else:
                s.set_log(field, args.logscale, linthresh=1.0e3)
        else:
            s.set_log(field, args.logscale)
        if dlog >= 2.0:
            s.set_zlim(field, -10.0**logmaxv, 10.0**logmaxv)
        else:
            s.set_zlim(field, -linmaxv, linmaxv)
    else:
        s.set_log(field, args.logscale)    
        s.set_cmap(field, args.colormap)
        if args.field_min:
            zmin = args.field_min
        else:
            zmin = 'min'
        if args.field_max:
            zmax = args.field_max
        else:
            zmax = 'max'
        s.set_zlim(field, zmin, zmax)

    # Annotations
    s.annotate_scale()
    if args.drawcells:
        s.annotate_cell_edges()
    if args.drawgrids:
        s.annotate_grids()

    # Sizing and saving
    s.set_buff_size(args.resolution)
    s.save('{}.slice.{}.{}.png'.format(args.infile, args.axis, field_short_name))

if __name__=="__main__":
    ds = yt.load(args.infile)

    if args.list_fields:
        print('Native fields in dataset:')
        for f in ds.field_list:
            print(f)
        print('Derived fields from dataset:')
        for f in ds.derived_field_list:
            print(f)
        exit()
    
    # Check axis input
    axes_list = ['x', 'y', 'z']
    if (args.axis != 'x' and
        args.axis != 'y' and
        args.axis != 'z'):
        print('Improper axis argument.')
        exit()

    if args.field:
        field = None
        field_short_name = None
        for f in ds.field_list + ds.derived_field_list:
            if f[1] == args.field:
                field_short_name = f[1]
                field = f
                break
        if not field:
            print('Field {} not present.'.format(args.field))
            exit()
        slicefield(ds, field, field_short_name)
    else:
        for f in ds.field_list:
            slicefield(ds, f, f[1])
