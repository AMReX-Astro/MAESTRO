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
parser.add_argument('-symlog', '--symlog', action='store_true', help='If supplied, use symlog scaling, which is linear near zero, to accomodate positive and negative values of enucdot.')
parser.add_argument('-logmax', '--logmax', type=float,
                    help='Log of the +/- maximum enucdot value.')
parser.add_argument('-rho', '--rhocontours', type=float, nargs='+', help='Draws contours for the densities provided (g/cm^3).')
parser.add_argument('-dc', '--drawcells', action='store_true', help='If supplied, draw the cell edges.')
parser.add_argument('-dg', '--drawgrids', action='store_true', help='If supplied, draw the grids.')
parser.add_argument('-lic', '--lic_velocity', action='store_true', help='If supplied, draw line-integral convolution for in-plane velocity field.')
parser.add_argument('-octant', '--octant', action='store_true', help='Sets slice view appropriately for octant dataset.')
parser.add_argument('-res', '--resolution', type=int, default=2048, help='Resolution for output plot.')
args = parser.parse_args()

# Check axis input
axes_list = ['x', 'y', 'z']
vel_list = []
if (args.axis != 'x' and
    args.axis != 'y' and
    args.axis != 'z'):
    print('Improper axis argument.')
    exit()
else:
    for i, ax in enumerate(axes_list):
        if ax != args.axis:
            vel_list.append('{}_vel'.format(ax))

@derived_field(name='pos_enucdot', units='erg/(g*s)')
def _pos_radial_velocity(field, data):
    return np.maximum(data[('boxlib','enucdot')],
                      yt.YTQuantity(1.0e-99, 'erg/(g*s)'))
@derived_field(name='neg_enucdot', units='erg/(g*s)')
def _neg_radial_velocity(field, data):
    return np.maximum(-data[('boxlib','enucdot')],
                      yt.YTQuantity(1.0e-99, 'erg/(g*s)'))

ds = yt.load(args.infile)

if not args.width:
    width = max(ds.domain_width)
else:
    width = yt.YTQuantity(args.width, 'cm')

pos_maxv = np.ceil(np.log10(ds.all_data().max('pos_enucdot')))
neg_maxv = np.ceil(np.log10(ds.all_data().max('neg_enucdot')))
maxv = max(max(pos_maxv, neg_maxv), 3)

if args.octant:
    dcenter = width.in_units('cm').v/2.0
    cpos    = ds.arr([dcenter, dcenter, dcenter], 'cm')
    s = yt.SlicePlot(ds, args.axis, ('boxlib', 'enucdot'),
                     center=cpos, width=width, origin="native")
else:
    s = yt.SlicePlot(ds, args.axis, ('boxlib', 'enucdot'),
                     center='c', width=width)

if args.rhocontours:
    for rhoc in args.rhocontours:
        rhounit = yt.YTQuantity(rhoc, 'g/(cm**3)')
        s.annotate_contour('density', ncont=1, clim=(rhounit, rhounit))

if args.logmax:
    log_max_val = args.logmax
else:
    log_max_val = maxv
    
if args.symlog:
    s.set_cmap(('boxlib', 'enucdot'), 'PiYG')
    s.set_log(('boxlib', 'enucdot'),  True, linthresh=10.0**(log_max_val-6))
    s.set_zlim(('boxlib', 'enucdot'), -10.0**log_max_val, 10.0**log_max_val)
else:
    s.set_cmap(('boxlib', 'enucdot'), 'PiYG')
    s.set_zlim(('boxlib', 'enucdot'), -10.0**log_max_val, 10.0**log_max_val)

#s.annotate_velocity(normalize=True)
#s.annotate_streamlines(('boxlib', 'y_vel'), ('boxlib', 'z_vel'), density=5)

if args.lic_velocity:
    s.annotate_line_integral_convolution(('boxlib', vel_list[0]),
                                         ('boxlib', vel_list[1]),
                                         lim=(0.5,0.65))

s.annotate_scale()

if args.drawcells:
    s.annotate_cell_edges()

if args.drawgrids:
    s.annotate_grids()

s.set_buff_size(args.resolution)
s.save('{}.slice.{}.enucdot.png'.format(args.infile, args.axis))
