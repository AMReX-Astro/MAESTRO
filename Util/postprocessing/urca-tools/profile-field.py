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
                    help='Number of bins. Defaults to 64.')
parser.add_argument('-f_min', '--field_min', type=float,
                    help='Minimum to use for plotting the field.')
parser.add_argument('-f_max', '--field_max', type=float,
                    help='Maximum to use for plotting the field.')
parser.add_argument('-ax_min', '--axis_min', type=float,
                    help='Minimum axis coordinate to use for plotting the field.')
parser.add_argument('-ax_max', '--axis_max', type=float,
                    help='Maximum axis coordinate to use for plotting the field.')
parser.add_argument('-o', '--outfile', type=str,
                    help='Name of output file to save. Defaults to "args.infile[0].[field]-vs-[axis].png"')
args = parser.parse_args()

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

    # Set vertical limits
    if args.field_min and args.field_max:
        plot.set_ylim(args.field, ymin=args.field_min, ymax=args.field_max)
    elif args.field_min:
        plot.set_ylim(args.field, ymin=args.field_min)
    elif args.field_max:
        plot.set_ylim(args.field, ymax=args.field_max)

    # Set horizontal limits
    if args.axis_min and args.axis_max:
        plot.set_xlim(args.field, xmin=args.axis_min, xmax=args.axis_max)
    elif args.axis_min:
        plot.set_xlim(args.field, xmin=args.axis_min)
    elif args.axis_max:
        plot.set_xlim(args.field, xmax=args.axis_max)

    # Save
    if args.outfile:
        plt_name = args.outfile
    else:
        plt_name = '{}.{}-vs-{}.png'.format(args.infile, args.field, args.axis)
    plot.save(plt_name)

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
