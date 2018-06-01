#!/usr/bin/env python
import yt
from yt.units import dimensions
from yt import derived_field
import yt.visualization.volume_rendering.api as vr
from yt.visualization.volume_rendering.transfer_function_helper import TransferFunctionHelper
from yt.visualization.volume_rendering.api import Scene, VolumeSource, Camera, ColorTransferFunction
import numpy as np
import argparse
from yt_urca_fields import UrcaShellFields, DatasetHelpers

parser = argparse.ArgumentParser()
parser.add_argument('infile', type=str, help='Name of input plotfile.')
parser.add_argument('-f', '--field', type=str, help='Name of the field to render, e.g. tfromp.')
parser.add_argument('-rup', '--rup', type=float, default=1.0e8, help='Maximum radius (cm). Default 1.0e8.')
parser.add_argument('-zoom', '--zoom', type=float, default=1.0, help='Camera zoom factor. Default 1.0.')
parser.add_argument('-cpos', '--camera_position', type=float, nargs=3, help='3-D Camera position in fractions of maximum radius (--rup).')
parser.add_argument('-cnorth', '--camera_north', type=float, nargs=3, help='Camera north vector (direction of up).')
parser.add_argument('-cmap', '--colormap', type=str, default='octarine', help='Name of the yt/matplotlib colormap. Defaults to octarine.')
parser.add_argument('-min', '--minimum', type=float, default=0.0, help='Minimum field value for transfer function.')
parser.add_argument('-max', '--maximum', type=float, default=1.0, help='Maximum field value for transfer function.')
parser.add_argument('-sig', '--sigma', type=float, default=0.08, help='Transfer function width parameter. (Default is 0.08).')
parser.add_argument('-n', '--num_layers', type=int, default=5, help='Number of layers. (Default is 5).')
parser.add_argument('-dd', '--drawdomain', action='store_true', help='If supplied, draw the boundaries of the domain.')
parser.add_argument('-dg', '--drawgrids', action='store_true', help='If supplied, draw the grids.')
parser.add_argument('-da', '--drawaxes', action='store_true', help='If supplied, draw an axes triad.')
parser.add_argument('-alpha_ones', '--alpha_ones', action='store_true', help='If supplied, set the transfer function values to ones.')
parser.add_argument('-res', '--resolution', type=int, default=2048, help='Resolution for output plot.')
parser.add_argument('-dry', '--dry_run', action='store_true', help='Plot only the transfer functions and quit.')
args = parser.parse_args()

# Open Dataset
ds = yt.load(args.infile)
core = ds.sphere(ds.domain_center, (args.rup, 'cm'))

# Load urca-specific fields
ushell_fields = UrcaShellFields()
ushell_fields.setup(ds)

# Get the field from the dataset
field, field_short_name = DatasetHelpers.get_field(ds, args.field)
assert(field)

# Create Scene
sc = Scene()

# Create Sources
so = VolumeSource(core, field)

bounds = np.array([args.minimum, args.maximum])
log_min = np.log10(args.minimum)
log_max = np.log10(args.maximum)
sigma  = args.sigma

nlayers = args.num_layers
if args.alpha_ones:
    alphavec = np.ones(nlayers)
else:
    alphavec = np.logspace(-3, 0, num=nlayers, endpoint=True)

tfh = TransferFunctionHelper(ds)
tfh.set_field(field)
print(field)
tfh.set_log(True)
tfh.grey_opacity = False
tfh.set_bounds(bounds)
tfh.build_transfer_function()
tfh.tf.add_layers(nlayers, colormap=args.colormap, w=sigma**2, alpha=alphavec)
tfh.plot('{}_tfun_{}.png'.format(args.infile, field_short_name))
so.transfer_function = tfh.tf

if args.dry_run:
    exit()

# Add sources to scene
sc.add_source(so)

# Add camera to scene
sc.add_camera()

# Set camera properties
sc.camera.focus = ds.domain_center
sc.camera.resolution = args.resolution
sc.camera.north_vector = yt.YTArray(args.camera_north, 'cm')
sc.camera.position = ds.domain_center + yt.YTArray(args.camera_position, 'cm') * args.rup
sc.camera.zoom(args.zoom)

# Annotate domain - draw boundaries
if args.drawdomain:
    sc.annotate_domain(ds, color=[1, 1, 1, 0.01])

# Annotate by drawing grids
if args.drawgrids:
    sc.annotate_grids(ds, alpha=0.01)

# Annotate by drawing axes triad
if args.drawaxes:
    sc.annotate_axes(alpha=0.01) 

# Render
sc.render()
sc.save('{}_rendering_{}.png'.format(args.infile, field_short_name), sigma_clip=3)
