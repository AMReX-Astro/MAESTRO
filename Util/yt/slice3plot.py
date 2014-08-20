#!/bin/env python

import argparse

import yt

import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

def doit(file, var, log):

    # load the data
    ds = yt.load(file)
    dd = ds.all_data()

    center = 0.5*(ds.domain_left_edge + ds.domain_right_edge)

    # see http://yt-project.org/docs/3.0/cookbook/complex_plots.html#multipanel-with-axes-labels

    fig = plt.figure()

    grid = AxesGrid(fig, (0.1, 0.1, 0.85, 0.85),
                    nrows_ncols = (1, 3),
                    axes_pad = 1.1,
                    label_mode = "all",
                    share_all = False,
                    cbar_location = "right",
                    cbar_mode = "each",
                    cbar_size = "3%",
                    cbar_pad = "0%")

    formatter = matplotlib.ticker.ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((-3,3))

    for i, d in enumerate(["x", "y", "z"]):

        p = yt.SlicePlot(ds, d, var, center=(center[0], center[1], center[2]),
                         origin="native", fontsize=10)
        p.set_log(var, log)

        plot = p.plots[var]
        plot.figure = fig
        plot.axes = grid[i].axes
        plot.cax = grid.cbar_axes[i]
        
        cb = plot.cb
        cb.formatter = formatter
        #cb.formatter.set_scientific(True)
        #cb.formatter.set_powerlimits((-3,3))
        #plot.cax.yaxis.set_major_formatter(formatter)

        #cb.update_ticks()

        p._setup_plots()

        ax = plot.axes
        ax.xaxis.set_major_formatter(formatter)
        ax.yaxis.set_major_formatter(formatter)

        ax.xaxis.offsetText.set_fontsize("small")
        ax.yaxis.offsetText.set_fontsize("small")

        #if i > 0:
        #    ax.yaxis.offsetText.set_visible(False)

        # cl = plt.getp(ax, 'ymajorticklabels')
        # plt.setp(cl, fontsize=10)
        # cl = plt.getp(ax, 'xmajorticklabels')
        # plt.setp(cl, fontsize=10)


        

    fig.set_size_inches(12.80, 7.20)

    plt.savefig("test.png")

    


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--log", help="plot the log of the variable", action="store_true")
    
    parser.add_argument("file", help="the name of the file to read", type=str)
    parser.add_argument("var", help="the name of the variable to plot", type=str)

    args = parser.parse_args()

    doit(args.file, args.var, args.log)

