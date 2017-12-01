#!/usr/bin/env python
import yt
from yt import derived_field
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('infile', type=str, help='Name of input plotfile.')
parser.add_argument('-w', '--width', type=float,
                    help='Width of slice (cm). Default is domain width.')
parser.add_argument('-symlog', '--symlog', action='store_true', help='If supplied, use symlog scaling, which is linear near zero, to accomodate positive and negative values of enucdot.')
parser.add_argument('-dc', '--drawcells', action='store_true', help='If supplied, draw the cell edges.')
parser.add_argument('-dg', '--drawgrids', action='store_true', help='If supplied, draw the grids.')
parser.add_argument('-octant', '--octant', action='store_true', help='Sets slice view appropriately for octant dataset.')
args = parser.parse_args()

@derived_field(name='pos_enucdot', units='erg/(g*s)')
def _pos_radial_velocity(field, data):
    return np.maximum(data[('boxlib','enucdot')], yt.YTQuantity(1.0e-99, 'erg/(g*s)'))
@derived_field(name='neg_enucdot', units='erg/(g*s)')
def _neg_radial_velocity(field, data):
    return np.maximum(-data[('boxlib','enucdot')], yt.YTQuantity(1.0e-99, 'erg/(g*s)'))
    
ds = yt.load(args.infile)

if not args.width:
    width = max(ds.domain_width)
else:
    width = yt.YTQuantity(args.width, 'cm')

pos_maxv = np.ceil(np.log10(ds.all_data().max('pos_enucdot')))
neg_maxv = np.ceil(np.log10(ds.all_data().max('neg_enucdot')))
maxv = max(pos_maxv, neg_maxv)

if args.octant:
    dcenter = width.in_units('cm').v/2.0
    cpos    = ds.arr([dcenter, dcenter, dcenter], 'cm')
    s = yt.SlicePlot(ds, 'x', ('boxlib', 'enucdot'),
                     center=cpos, width=width, origin="native")
else:
    s = yt.SlicePlot(ds, 'x', ('boxlib', 'enucdot'),
                     center='c', width=width, origin="native")

if args.symlog:
    s.set_cmap(('boxlib', 'enucdot'), 'PiYG')
    s.set_log(('boxlib', 'enucdot'), True, linthresh=1.0e3)
    s.set_zlim(('boxlib', 'enucdot'), -10.0**maxv, 10.0**maxv)
else:
    s.set_cmap(('boxlib', 'enucdot'), 'Greens')    
#s.annotate_velocity(normalize=True)
#s.annotate_streamlines(('boxlib', 'y_vel'), ('boxlib', 'z_vel'), density=5)
s.annotate_line_integral_convolution(('boxlib', 'y_vel'), ('boxlib', 'z_vel'), lim=(0.5,0.65))
s.annotate_scale()

if args.drawcells:
    s.annotate_cell_edges()

if args.drawgrids:
    s.annotate_grids()

s.set_buff_size(2048)
s.save('{}.slice.enucdot.png'.format(args.infile))
