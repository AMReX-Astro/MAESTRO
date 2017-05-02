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
args = parser.parse_args()

@derived_field(name='pos_radial_velocity', units='cm/s')
def _pos_radial_velocity(field, data):
    return np.maximum(data[('boxlib','radial_velocity')], 1.0e-99)
@derived_field(name='neg_radial_velocity', units='cm/s')
def _neg_radial_velocity(field, data):
    return np.maximum(-data[('boxlib','radial_velocity')], 1.0e-99)

ds = yt.load(args.infile)

if not args.width:
    width = max(ds.domain_width)
else:
    width = (args.width, 'cm')
    
pos_maxv = np.ceil(np.log10(ds.all_data().max('pos_radial_velocity')))
neg_maxv = np.ceil(np.log10(ds.all_data().max('neg_radial_velocity')))
maxv = max(pos_maxv, neg_maxv)
s = yt.SlicePlot(ds, 'x', ('boxlib', 'radial_velocity'), center='c', width=width)
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

s = yt.SlicePlot(ds, 'x', ('boxlib', 'circum_velocity'), center='c', width=width)
s.annotate_velocity(normalize=True)
s.annotate_scale()

if args.drawcells:
    s.annotate_cell_edges()

if args.drawgrids:
    s.annotate_grids()

s.set_buff_size(2048)
s.save('{}.slice.vel_circum.png'.format(args.infile))
