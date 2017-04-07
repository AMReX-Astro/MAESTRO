from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt

# Angle-averaged profile Class
class AngleAveragedProfile(object):
    def __init__(self, filename=None):
        self.init_vars()
        if filename:
            self.read_from_file(filename)
            
    def init_vars(self):
        self.header = {}
        self.data = {}
        self.data_keys = []
        self.filename = ''
        
    def read_from_file(self, filename):
        # Clear my variables
        self.init_vars()
        self.filename = filename
        # Given a profile filename, read the profile
        f = open(filename, 'r')
        # Get the header global values
        num_sep = 0 # Number of '-----' separator lines encountered
        readColumnLabels = False # True if next line has column labels
        for line in f:
            ls = line.strip()
            if ls[0] == '#' or readColumnLabels:
                if readColumnLabels:
                    cl = []
                    if ls[0] == '#':
                        ls = ls[1:]
                    for ci in ls.split('  '):
                        ci = ci.strip()
                        # Need to check because splitting by '  '
                        # can yield ' ' which strip to ''.
                        if ci:
                            # If field is enclosed in [ ] brackets, remove them
                            if ci[0] == '[' and ci[-1] == ']':
                                ci = ci[1:-1].strip()
                            cl.append(ci) 
                    self.data_keys = cl
                    print(self.data_keys)
                    for ci in cl:
                        self.data[ci] = []
                    readColumnLabels = False
                    num_sep = 0
                else:
                    # See if there's an equals sign and get value
                    if '=' in ls:
                        k, v = ls[2:].split('=', 1)
                        k = k.strip()
                        v = [float(vf) for vf in v.strip().split()]
                        if len(v) > 1:
                            v = np.array(v)
                        else:
                            v = v[0]
                        self.header[k] = v
                    elif '----------------------------' in ls:
                        num_sep += 1
                        if num_sep == 3 and not readColumnLabels:
                            readColumnLabels = True
            else:
                # Read data line
                ld = [float(di) for di in ls.split()]
                for k, v in zip(self.data_keys, ld):
                    self.data[k].append(v)
        f.close()
        # Turn data into numpy arrays
        for k in self.data.keys():
            self.data[k] = np.array(self.data[k])

    def gen_tick_spacing(self):
        # Generate possible tick spacings
        initvals = [0.25, 0.5, 1.0, 5.0]
        n = 0
        if initvals:
            for v in initvals:
                yield v
        else:
            while(True):
                n += 1
                yield float(10*n)

    def plot_var(self, var, fmt='png', rup=None):
        # Plot the variable corresponding to the data key var
        # Independent axis is radius r
        # Plots are log scale on the dependent axis
        if var not in self.data.keys():
            return
        fig = plt.figure()
        ax = fig.add_subplot(111)
        idxup = -1
        if rup:
            ax.set_xlim([0, rup])
            # Get the lowest index where radius > rup
            idxup = np.where(self.data['r'] > rup)[0][0]
        neg_idx = np.where(self.data[var][:idxup] < 0.0)
        pos_idx = np.where(self.data[var][:idxup] > 0.0)
        ax.set_xlabel('r')
        # find smallest non-zero log10 magnitude in quantity to plot
        try:
            neg_min = np.log10(np.amin(np.absolute(self.data[var][neg_idx])))
        except:
            neg_min = None
        try:
            pos_min = np.log10(np.amin(np.absolute(self.data[var][pos_idx])))
        except:
            pos_min = None
        if pos_min and neg_min:
            lwb = min(neg_min, pos_min)
        else:
            if pos_min:
                lwb = pos_min
            elif neg_min:
                lwb = neg_min
            else:
                lwb = None
        upb = np.log10(np.amax(np.absolute(self.data[var][:idxup])))
        if (not lwb) or upb-lwb <= 1.0:
            # plot quantity on linear axis
            # plot linear scale magnitudes
            ax.plot(self.data['r'][:idxup], self.data[var][:idxup], color='green')
            # plot positive points in blue
            ax.plot(self.data['r'][:idxup][pos_idx], self.data[var][:idxup][pos_idx],
                    linestyle='None', marker='^', color='blue', markersize=8, alpha=0.5)
            # plot negative points in red
            ax.plot(self.data['r'][:idxup][neg_idx], self.data[var][:idxup][neg_idx],
                    linestyle='None', marker='v', color='red', markersize=8, alpha=0.5)
            ax.set_ylabel('$\mathrm{' + var.replace('_','\_') + '}$')                        
        else:
            # plot quantity on log10 axis
            # plot log scale magnitudes
            ax.plot(self.data['r'][:idxup], np.log10(np.absolute(self.data[var][:idxup])), color='green')
            # plot positive points in blue
            ax.plot(self.data['r'][:idxup][pos_idx], np.log10(self.data[var][:idxup][pos_idx]),
                    linestyle='None', marker='^', color='blue', markersize=8, alpha=0.5)
            # plot negative points in red
            ax.plot(self.data['r'][:idxup][neg_idx], np.log10(np.absolute(self.data[var][:idxup][neg_idx])),
                    linestyle='None', marker='v', color='red', markersize=8, alpha=0.5)
            ax.set_ylabel('$\mathrm{Log_{10} \ | ' + var.replace('_','\_') + ' |}$')
            upb = np.ceil(upb*10.0)/10.0
            lwb = np.floor(lwb*10.0)/10.0
            yticks = None
            for tspac in self.gen_tick_spacing():
                nticks = int(np.floor((upb-lwb)/tspac) + 1)
                eps = upb - (lwb + tspac*(nticks-2))
                if nticks <= 10 and eps > 0.5*tspac:
                    yticks = np.array([lwb + tspac*(j) for j in range(nticks-1)] + [upb])
                    break
            ax.set_yticks(yticks)
            ax.set_ylim((lwb, upb))
        # List the time above the plot
        tart = ax.text(1.0, 1.01, 'time = {}'.format(self.header['time']),
                       transform=ax.transAxes,
                       verticalalignment='bottom',
                       horizontalalignment='right')
        outname = '.'.join([self.filename, var.replace(' ','-'), fmt])
        if fmt=='png':
            plt.savefig(outname, bbox_extra_artists=(tart,), dpi=300)
        else:
            plt.savefig(outname, bbox_extra_artists=(tart,))
        plt.close(fig)
            
    def plot_all_vars(self, fmt='png', rup=None):
        # Plot all variables in the profile and save
        # fmt is the suffix passed to savefig
        for var in self.data.keys():
            if var != 'r':
                self.plot_var(var, fmt, rup)
