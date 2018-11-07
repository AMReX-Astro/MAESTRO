#!/usr/bin/env python
"""
Use yt to get the sum of a field in a boxlib plotfile.

Donald E. Willcox
"""
import yt
from yt import derived_field
import numpy as np
import argparse
from yt_urca_fields import UrcaShellFields

parser = argparse.ArgumentParser()
parser.add_argument('infile', type=str, help='Name of input plotfile.')
parser.add_argument('-f', '--field', type=str, default='cell_volume',
                    help='Name of the field to sum. Eg. "tfromp". Default is cell volume.')
parser.add_argument('-r', '--radius', type=str,
                    help='Radius in centimeters of a centered sphere over which to take the sum. Default is the entire domain.')
args = parser.parse_args()

def get_field(ds):
    field = None
    field_short_name = None
    for f in ds.field_list + ds.derived_field_list:
        if f[1] == args.field:
            field_short_name = f[1]
            field = f
            return field, field_short_name
    if not field:
        print('Field {} not present.'.format(args.field))
        return None, None

if __name__=="__main__":
    ds = yt.load(args.infile)

    # Add Urca fields
    ushell_fields = UrcaShellFields()
    ushell_fields.setup(ds)

    field, field_short_name = get_field(ds)
    assert(field)

    region = None
    if args.radius:
        region = ds.sphere('c', (args.radius, 'cm'))
    else:
        region = ds.all_data()

    fsum = region.quantities.total_quantity(field)

    if args.radius:
        print('The sum of {} within a sphere of radius {} cm is: {}'.format(args.field, args.radius, fsum))
    else:
        print('The domain sum of {} is: {}'.format(args.field, fsum))
