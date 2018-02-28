#!/usr/bin/env python
import yt
from yt import derived_field
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('infile', type=str, help='Name of input plotfile.')
parser.add_argument('-w', '--width', type=float,
                    help='Width of slice (cm). Default is domain width.')
parser.add_argument('-axis', '--axis', type=str, default='x', help='Axis along which to slice. Should be "x", "y", or "z".')
parser.add_argument('-dc', '--drawcells', action='store_true', help='If supplied, draw the cell edges.')
parser.add_argument('-dg', '--drawgrids', action='store_true', help='If supplied, draw the grids.')
parser.add_argument('-octant', '--octant', action='store_true', help='Sets slice view appropriately for octant dataset.')
parser.add_argument('-rho', '--rhocontours', type=float, nargs='+', help='Draws contours for the densities provided (g/cm^3).')
parser.add_argument('-res', '--resolution', type=int, default=2048, help='Resolution for output plot.')
parser.add_argument('-vel', '--velocityarrows', action='store_true', help='If supplied, annotate velocity arrows.')
args = parser.parse_args()

@derived_field(name='pos_radial_velocity', units='cm/s')
def _pos_radial_velocity(field, data):
    return np.maximum(data[('boxlib','radial_velocity')],
                      yt.YTQuantity(1.0e-99, 'cm/s'))
@derived_field(name='neg_radial_velocity', units='cm/s')
def _neg_radial_velocity(field, data):
    return np.maximum(-data[('boxlib','radial_velocity')],
                      yt.YTQuantity(1.0e-99, 'cm/s'))

# Check axis input
axes_list = ['x', 'y', 'z']
if (args.axis != 'x' and
    args.axis != 'y' and
    args.axis != 'z'):
    print('Improper axis argument.')
    exit()

ds = yt.load(args.infile)

if not args.width:
    width = max(ds.domain_width)
else:
    width = yt.YTQuantity(args.width, 'cm')
    
pos_maxv = np.ceil(np.log10(ds.all_data().max('pos_radial_velocity')))
neg_maxv = np.ceil(np.log10(ds.all_data().max('neg_radial_velocity')))
maxv = max(pos_maxv, neg_maxv)

if args.octant:
    dcenter = width.in_units('cm').v/2.0
    cpos    = ds.arr([dcenter, dcenter, dcenter], 'cm')
    s = yt.SlicePlot(ds, args.axis, ('boxlib', 'radial_velocity'),
                     center=cpos, width=width, origin='native')
else:
    s = yt.SlicePlot(ds, args.axis, ('boxlib', 'radial_velocity'),
                     center='c', width=width)

if args.rhocontours:
    for rhoc in args.rhocontours:
        rhounit = yt.YTQuantity(rhoc, 'g/(cm**3)')
        s.annotate_contour('density', ncont=1, clim=(rhounit, rhounit))

s.set_cmap(('boxlib', 'radial_velocity'), 'RdBu')
s.set_log(('boxlib', 'radial_velocity'), True, linthresh=1.0e3)
s.set_zlim(('boxlib', 'radial_velocity'), -10.0**maxv, 10.0**maxv)

if args.velocityarrows:
    s.annotate_velocity(normalize=True)

s.annotate_scale()

if args.drawcells:
    s.annotate_cell_edges()

if args.drawgrids:
    s.annotate_grids()

s.set_buff_size(args.resolution)
s.save('{}.slice.{}.vel_radial.png'.format(args.infile, args.axis))

if args.octant:
    dcenter = width.in_units('cm').v/2.0
    cpos    = ds.arr([dcenter, dcenter, dcenter], 'cm')
    s = yt.SlicePlot(ds, args.axis, ('boxlib', 'circum_velocity'),
                     center=cpos, width=width, origin="native")
else:
    s = yt.SlicePlot(ds, args.axis, ('boxlib', 'circum_velocity'),
                     center='c', width=width)

if args.rhocontours:
    for rhoc in args.rhocontours:
        rhounit = yt.YTQuantity(rhoc, 'g/(cm**3)')
        s.annotate_contour('density', ncont=1, clim=(rhounit, rhounit))

if args.velocityarrows:
    s.annotate_velocity(normalize=True)

s.annotate_scale()

if args.drawcells:
    s.annotate_cell_edges()

if args.drawgrids:
    s.annotate_grids()

s.set_buff_size(args.resolution)
s.save('{}.slice.{}.vel_circum.png'.format(args.infile, args.axis))
