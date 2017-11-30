import yt
from yt import derived_field
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('infile', type=str, help='Name of input plotfile.')
parser.add_argument('-w', '--width', type=float,
                    help='Width of slice (cm). Default is domain width.')
parser.add_argument('-dc', '--drawcells', action='store_true', help='If supplied, draw the cell edges.')
parser.add_argument('-dg', '--drawgrids', action='store_true', help='If supplied, draw the grids.')
parser.add_argument('-octant', '--octant', action='store_true', help='Sets slice view appropriately for octant dataset.')
args = parser.parse_args()

@derived_field(name='pos_radial_velocity', units='cm/s')
def _pos_radial_velocity(field, data):
    return np.maximum(data[('boxlib','radial_velocity')], yt.YTQuantity(1.0e-99, 'cm/s'))
@derived_field(name='neg_radial_velocity', units='cm/s')
def _neg_radial_velocity(field, data):
    return np.maximum(-data[('boxlib','radial_velocity')], yt.YTQuantity(1.0e-99, 'cm/s'))

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
    s = yt.SlicePlot(ds, 'x', ('boxlib', 'radial_velocity'),
                     center=cpos, width=width, origin="native")
else:
    s = yt.SlicePlot(ds, 'x', ('boxlib', 'radial_velocity'),
                     center='c', width=width, origin="native")

s.set_cmap(('boxlib', 'radial_velocity'), 'RdBu')
s.set_log(('boxlib', 'radial_velocity'), True, linthresh=1.0e3)
s.set_zlim(('boxlib', 'radial_velocity'), -10.0**maxv, 10.0**maxv)
s.annotate_scale()

if args.drawcells:
    s.annotate_cell_edges()

if args.drawgrids:
    s.annotate_grids()

s.set_buff_size(2048)
s.save('{}.slice.vel_radial.png'.format(args.infile))

if args.octant:
    dcenter = width.in_units('cm').v/2.0
    cpos    = ds.arr([dcenter, dcenter, dcenter], 'cm')
    s = yt.SlicePlot(ds, 'x', ('boxlib', 'circum_velocity'),
                     center=cpos, width=width, origin="native")
else:
    s = yt.SlicePlot(ds, 'x', ('boxlib', 'circum_velocity'),
                     center='c', width=width, origin="native")

s.annotate_velocity(normalize=True)
s.annotate_scale()

if args.drawcells:
    s.annotate_cell_edges()

if args.drawgrids:
    s.annotate_grids()

s.set_buff_size(2048)
s.save('{}.slice.vel_circum.png'.format(args.infile))
