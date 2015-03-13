#!/usr/bin/env python

# a simple script to plot 2-d XRB data of a bunch of variables all together in
# the same figure.  Note: the variables and limits are hardcoded here.
#
# 2012-01-10 M. Zingale
import matplotlib
matplotlib.use('Agg')
import fsnapshot
import numpy
import pylab
import os
import sys
import getopt
import math
import string
from mpl_toolkits.axes_grid1 import ImageGrid


# data
var1 = "tfromh"
m1 = 5.e7
M1 = 8.e8
l1 = 1

var2 = "vort"
m2 = -2.e6
M2 = 2.e6
l2 = 0

var3 = "magvel"
m3 = 1.e4
M3 = 1.e7
l3 = 1

var4 = "enucdot"
m4 = 1.e14
M4 = 1.e19
l4 = 1
epsEnuc = 1.e-10

var5 = "X(Mg22)"
m5 = 1.e-8
M5 = 1.e-3
l5 = 1 


#==============================================================================
# do_plot
#==============================================================================
def do_plot(plotfile,
            eps, dpi, 
            xmax_pass, ymax_pass, zmax_pass):

    global var1, m1, M1, l1, var2, m2, M2, l2, var3, m3, M3, l3, var4, m4, M4, l4, var5, m5, M5, l5, epsEnuc

    pylab.rc("font", size=9)


    #--------------------------------------------------------------------------
    # construct the output file name
    #--------------------------------------------------------------------------
    outFile = os.path.normpath(plotfile) + "_all"

    if (not eps):
        outFile += ".png"

    else:
        outFile += ".eps"


    #--------------------------------------------------------------------------
    # read in the meta-data from the plotfile
    #--------------------------------------------------------------------------
    (nx, ny, nz) = fsnapshot.fplotfile_get_size(plotfile)

    time = fsnapshot.fplotfile_get_time(plotfile)

    (xmin, xmax, ymin, ymax, zmin, zmax) = \
        fsnapshot.fplotfile_get_limits(plotfile)

    dx = (xmax - xmin)/nx
    x = xmin + numpy.arange( (nx), dtype=numpy.float64 )*dx

    dy = (ymax - ymin)/ny
    y = ymin + numpy.arange( (ny), dtype=numpy.float64 )*dy

    if (nz > 0):
        sys.exit("ERROR: 3-d not supported")


    extent = [xmin, xmax, ymin, ymax]

    if (not xmax_pass == None):
        extent[1] = xmax_pass

    if (not ymax_pass == None):
        extent[3] = ymax_pass

    ix = nx
    if (not xmax_pass == None):
        ix = int((xmax_pass - xmin)/dx)

    iy = ny
    if (not ymax_pass == None):
        iy = int((ymax_pass - ymin)/dy)


    sparseX = 0
    if (ny >= 3*nx):
        sparseX = 1

    # read in the data

    # 1
    data1 = numpy.zeros( (nx, ny), dtype=numpy.float64)

    (data1, err) = fsnapshot.fplotfile_get_data_2d(plotfile, var1, data1)
    if (not err == 0):
        sys.exit(2)

    data1 = numpy.transpose(data1)

    if (l1):
        data1 = numpy.log10(data1)
        m1 = math.log10(m1)
        M1 = math.log10(M1)

    # 2
    data2 = numpy.zeros( (nx, ny), dtype=numpy.float64)

    (data2, err) = fsnapshot.fplotfile_get_data_2d(plotfile, var2, data2)
    if (not err == 0):
        sys.exit(2)

    data2 = numpy.transpose(data2)

    if (l2):
        data1 = numpy.log10(data2)
        m2 = math.log10(m2)
        M2 = math.log10(M2)

    # 3
    data3 = numpy.zeros( (nx, ny), dtype=numpy.float64)

    (data3, err) = fsnapshot.fplotfile_get_data_2d(plotfile, var3, data3)
    if (not err == 0):
        sys.exit(2)

    data3 = numpy.transpose(data3)

    if (l3):
        data3 = numpy.log10(data3)
        m3 = math.log10(m3)
        M3 = math.log10(M3)

    # 4
    data4 = numpy.zeros( (nx, ny), dtype=numpy.float64)

    (data4, err) = fsnapshot.fplotfile_get_data_2d(plotfile, var4, data4)
    if (not err == 0):
        sys.exit(2)

    data4 = numpy.transpose(data4)

    # floor values
    data4[data4 < epsEnuc] = epsEnuc

    if (l4):
        data4 = numpy.log10(data4)
        m4 = math.log10(m4)
        M4 = math.log10(M4)

    # 5
    data5 = numpy.zeros( (nx, ny), dtype=numpy.float64)

    (data5, err) = fsnapshot.fplotfile_get_data_2d(plotfile, var5, data5)
    if (not err == 0):
        sys.exit(2)

    data5 = numpy.transpose(data5)

    if (l5):
        data5 = numpy.log10(data5)
        m5 = math.log10(m5)
        M5 = math.log10(M5)


    # use an imagegrid to place all the axes close together.
    # see http://matplotlib.org/examples/axes_grid/demo_axes_grid2.html
    F = pylab.figure(1, (12, 8))
    F.clf()

    formatter = matplotlib.ticker.ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((-3,3))

    grid = ImageGrid(F, 111, # similar to subplot(111)
                     nrows_ncols = (1, 5),
                     direction="row",
                     axes_pad = 0.2 ,
                     add_all=True,
                     label_mode = "L",
                     share_all = True,
                     cbar_location="top",
                     cbar_mode="each",
                     cbar_size="3%",
                     cbar_pad="5%",
                     )

    # variable 1
    ax = grid[0]
    im=ax.imshow(data1[0:iy,0:ix], origin='lower', extent=extent, vmin=m1, vmax=M1)
    ax.set_title(var1)
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    ax.xaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))
    ax.yaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))

    cb = ax.cax.colorbar(im, format=formatter)
    #ax.cax.set_title(var1)
        
    # variable 2
    ax = grid[1]
    im=ax.imshow(data2[0:iy,0:ix], origin='lower', extent=extent, vmin=m2, vmax=M2)
    ax.set_title(var2)
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    ax.xaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))
    ax.yaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))
    ax.set_xticks([xmin, xmin + (xmax-xmin)/3.0, xmin + 2.0*(xmax-xmin)/3.0])

    ax.cax.colorbar(im, format=formatter)

    # variable 3
    ax = grid[2]
    im=ax.imshow(data3[0:iy,0:ix], origin='lower', extent=extent, vmin=m3, vmax=M3)
    ax.set_title(var3)
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    ax.xaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))
    ax.yaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))

    ax.cax.colorbar(im, format=formatter)


    # variable 4
    ax = grid[3]
    im=ax.imshow(data4[0:iy,0:ix], origin='lower', extent=extent, vmin=m4, vmax=M4)
    ax.set_title(var4)
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    ax.xaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))
    ax.yaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))

    ax.cax.colorbar(im, format=formatter)


    # variable 5
    ax = grid[4]
    im=ax.imshow(data5[0:iy,0:ix], origin='lower', extent=extent, vmin=m5, vmax=M5)
    ax.set_title(var5)
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    ax.xaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))
    ax.yaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))

    ax.cax.colorbar(im, format=formatter)


    # write the time
    pylab.text(0.0,0.05, "t = %g s" % (time), transform = F.transFigure)


    #--------------------------------------------------------------------------
    # save the figure
    #--------------------------------------------------------------------------
    try: pylab.tight_layout()  # requires matplotlib > 1.1
    except:
        pass

    if (not eps):
        pylab.savefig(outFile, bbox_inches='tight', dpi=dpi, pad_inches=0.33)
    else:
        pylab.savefig(outFile, bbox_inches='tight', pad_inches=0.33)



#==============================================================================
# usage
#==============================================================================
def usage():
    usageStr = """
    ./plotsinglevar.py [options] plotfile component [component2]

    Make a simple colormap plot of variable "component" from the
    BoxLib plotfile "plotfile".  Support for 2-d and 3-d.

    If two component names are specified after the plotfile name then
    two side-by-side plots will be made, one for each specified
    component.

    Options:

       -o outfile   save the plot to the file outfile

       -m value     set the minimum data range for the plot to value

       -M value     set the maximum data range for the plot to value

       --log        plot the logarithm (base-10) of the data
                    (note: if the data < 0 anywhere, then the abs
                     is taken first)

       --eps        make an EPS plot instead of a PNG

       --dpi value  (PNG only) make the plot with the dpi specified by
                    value

       --origin     (3-d only) slice through the origin (0,0,0) instead
                    of the center of the domain.

       --annotate string (2-d only) add an annotation under the time


    Note: this script requires the fsnapshot.so library, compiled with
    f2py using the GNUmakefile in data_processing/python_plotfile/

    """              
    print usageStr



#==============================================================================
# main
#==============================================================================
if __name__== "__main__":

    eps = 0
    dpi = 100
    xmax = None
    ymax = None
    zmax = None

    try: opts, next = getopt.getopt(sys.argv[1:], "X:Y:Z:", 
                                    ["eps","dpi="])
    except getopt.GetoptError:
        print "invalid calling sequence"
        usage()
        sys.exit(2) 
               

    for o, a in opts:

        if o == "-X":
            try: xmax = float(a)
            except ValueError:
                print "invalid value for -X"
                sys.exit(2)            

        if o == "-Y":
            try: ymax = float(a)
            except ValueError:
                print "invalid value for -Y"
                sys.exit(2)            

        if o == "-Z":
            try: zmax = float(a)
            except ValueError:
                print "invalid value for -Z"
                sys.exit(2)            
 
        if o == "--eps":
            eps = 1

        if o == "--dpi":
            try: dpi = int(a)
            except ValueError:
                print "invalid value for --dpi"
                sys.exit(2)



    try: plotfile = next[0]
    except IndexError:
        print "ERROR: plotfile not specified"
        usage()
        sys.exit(2)

    do_plot(plotfile, 
            eps, dpi, 
            xmax, ymax, zmax)
