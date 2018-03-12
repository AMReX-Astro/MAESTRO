#!/usr/bin/env python
"""
Use yt to profile a given field averaged versus a given axis.

Donald E. Willcox
"""
import yt
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('infiles', type=str, nargs='+', help='Name of input plotfile.')
parser.add_argument('-labels', '--labels', type=str, nargs='+',
                    help='Labels to use corresponding to each input file.')
parser.add_argument('-f', '--field', type=str, default='enucdot',
                    help='Name of the field to profile. Eg. "enucdot".')
parser.add_argument('-f_log', '--field_log', action='store_true',
                    help='If supplied, use log scaling for field.')
parser.add_argument('-ax', '--axis', type=str, default='radius',
                    help='Axis (field) along which to profile the desired field. Default is "x".')
parser.add_argument('-ax_log', '--axis_log', action='store_true',
                    help='If supplied, use log scaling for axis.')
parser.add_argument('-nbins', '--nbins', type=int, default=64,
                    help='Number of bins.')
# parser.add_argument('-f_min', '--field_min', type=float,
#                     help='Minimum to use for plotting the field.')
# parser.add_argument('-f_max', '--field_max', type=float,
#                     help='Maximum to use for plotting the field.')
parser.add_argument('-o', '--outfile', type=str, default='profile.png',
                    help='Name of output file to save. Defaults to "profile.png"')
args = parser.parse_args()

# def doit():
#     ds = yt.load(args.infile)
#     ad = ds.all_data()

#     plt = yt.ProfilePlot(ad, args.axis, args.field,
#                          n_bins=args.nbins,
#                          x_log=args.axis_log,
#                          y_log={args.field:args.field_log})

#     if args.field_min and args.field_max:
#         plt.set_ylim(args.field, ymin=args.field_min, ymax=args.field_max)
#     elif args.field_min:
#         plt.set_ylim(args.field, ymin=args.field_min)
#     elif args.field_max:
#         plt.set_ylim(args.field, ymax=args.field_max)

#     plt_name = '{}.{}-vs-{}.png'.format(args.infile, args.field, args.axis)
#     plt.save(plt_name)

def create_profile(infile):
    ds = yt.load(infile)
    ad = ds.all_data()
    profile = yt.create_profile(ad, args.axis, fields=[args.field],
                                weight_field=None, n_bins=args.nbins,
                                logs={args.axis: args.axis_log,
                                      args.field: args.field_log})
    return profile

def doit(labels):
    profiles = [create_profile(f) for f in args.infiles]
    plot = yt.ProfilePlot.from_profiles(profiles, labels=labels)
    plot.save(args.outfile)

if __name__=='__main__':
    if args.labels:
        if len(args.labels) != len(args.infiles):
            print('Enter the same number of labels as input files.')
            exit()
        else:
            labels = args.labels
    else:
        labels = args.infiles
    doit(labels)
