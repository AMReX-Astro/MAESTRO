#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument('infiles', type=str, nargs='+',
                    help='Name of diagnostic files to plot.')
parser.add_argument('-tmin', '--tmin', type=float,
                    help='Minimum time to plot.')
args = parser.parse_args()

class RunDiagnostics(object):
    def __init__(self, date=None, time=None,
                 directory=None, job=None, fields=[], data={}):
        # Holds the run diagnostics for one run
        self.output_date = date
        self.output_time = time
        self.output_dir = directory
        self.job_name = job
        self.fields = fields
        self.data = data
        for f in self.fields:
            if not f in self.data.keys():
                self.data[f] = []

    def __eq__(self, other):
        # Two RunDiagnostics are equal if their metadata is equal.
        # THEY CAN HAVE DIFFERENT DATA
        equal = True
        equal = equal and self.output_date == other.output_date
        equal = equal and self.output_time == other.output_time
        equal = equal and self.output_dir == other.output_dir
        equal = equal and self.job_name == other.job_name
        return equal

    def __add__(self, other):
        # Adds two RunDiagnostics together by checking if other has fields self does not
        # and other things being equal, returning a new RunDiagnostics with the fields of both
        # time values should be equal, however.
        # Return None if the RunDiagnostics cannot be added (i.e. times are not equal)
        if self == other:
            fields = []
            data = {}
            for k in self.data.keys():
                if not k in fields:
                    data[k] = self.data[k]
                    fields.append(k)
            for k in other.data.keys():
                if not k in fields:
                    data[k] = other.data[k]
                    fields.append(k)
            return RunDiagnostics(date=self.output_date,
                                  time=self.output_time,
                                  directory=self.output_dir,
                                  job=self.job_name,
                                  fields=fields,
                                  data=data)
        else:
            # if the times are different but the fields are the same, then append
            if set(self.data.keys()) == set(other.data.keys()) and self.job_name == other.job_name:
                fields = self.data.keys()
                data = {}
                for f in fields:
                    data[f] = []
                # determine which has the earlier time data
                self_time_zero = self.data['time'][0]
                other_time_zero = other.data['time'][0]
                if self_time_zero < other_time_zero:
                    first = self
                    second = other
                    second_start_time = other_time_zero
                else:
                    first = other
                    second = self
                    second_start_time = self_time_zero
                # fill data from first up to second_start_time
                for i in range(len(first.data['time'])):
                    if first.data['time'][i] >= second_start_time:
                        break
                    for f in fields:
                        data[f].append(first.data[f][i])
                # continue filling data from second
                print(len(second.data['time']))
                for i in range(len(second.data['time'])):
                    print(i)
                    for f in fields:
                        print(len(second.data[f]))
                        data[f].append(second.data[f][i])
                # Return a new RunDiagnostics with this new data
                return RunDiagnostics(job=self.job_name, fields=fields, data=data)
            else:
                return None

    def read(self, line):
        # Read this line data
        ls = line.split()
        assert(len(ls)==len(self.fields))
        for f, lsi in zip(self.fields, ls):
            lsf = float(lsi)
            self.data[f].append(lsf)

    def plot_field(self, t, y, ylabel):
        fig, ax = plt.subplots()
        ax.plot(np.log10(t), y)
        ax.set_xlabel('Log time (s)')
        ax.set_ylabel('{}'.format(ylabel))
        plt.savefig('{}.diag.{}.png'.format(args.infile, ylabel.replace(' ', '_')))
        plt.clf()

    def plot(self):
        # Plot each column data vs. time
        for k in self.data.keys():
            self.data[k] = np.array(self.data[k])
        ilow = 0
        if args.tmin:
            for i, t in enumerate(self.data['time']):
                if t >= args.tmin:
                    ilow = i
                    break
        for c in self.data.keys():
            if c != 'time':
                print('plotting {}'.format(c))
                self.plot_field(self.data['time'][ilow:], self.data[c][ilow:], c)

class FileDiagnostics(object):
    re_date = re.compile('# output date: (.*)')
    re_time = re.compile('# output time: (.*)')
    re_dir  = re.compile('# output dir: (.*)')
    re_job  = re.compile('# job name: (.*)')
    re_fields = re.compile('#                time')
    
    def __init__(self, filename=None, diagnostics=[]):
        # Holds the Diagnostics in a given file
        self.diagnostics = diagnostics
        if filename:
            self.read(filename)

    def __add__(self, other):
        # Adds two FileDiagnostics together
        if self.diagnostics and other.diagnostics:
            assert(len(self.diagnostics)==len(other.diagnostics))
            new_diagnostics = []
            for ds, do in zip(self.diagnostics, other.diagnostics):
                new = ds + do
                assert(new)
                new_diagnostics.append(new)
            return FileDiagnostics(diagnostics=new_diagnostics)
        elif self.diagnostics:
            return self
        elif other.diagnostics:
            return other
        else:
            return FileDiagnostics()

    def combined_diagnostics(self):
        # Returns the combined data from all the RunDiagnostics,
        # properly interleaved.
        if len(self.diagnostics) == 0:
            return None
        elif len(self.diagnostics) == 1:
            return self.diagnostics[0]
        else:
            sum_diag = self.diagnostics[0]
            for diag in self.diagnostics[1:]:
                sum_diag += diag
            return sum_diag

    def read(self, filename):
        # Read diagnostics from filename
        fin = open(filename, 'r')
        flines = []
        for line in fin:
            ls = line.strip()
            if ls:
                flines.append(ls)
        fin.close()

        diagnostics = None
        while True:
            if len(flines) == 0:
                break

            line = flines.pop(0)
            m = FileDiagnostics.re_date.match(line)
            if m:
                # Get diagnostics object info
                date = m.group(1)
                time = FileDiagnostics.re_time.match(flines.pop(0)).group(1)
                jdir = FileDiagnostics.re_dir.match(flines.pop(0)).group(1)
                job  = FileDiagnostics.re_job.match(flines.pop(0)).group(1)
                fields = []
                field_string = flines.pop(0).split('  ')
                for fs in field_string:
                    fs = fs.strip()
                    if fs and fs != '#':
                        fields.append(fs)

                # Check to see if we already have this diagnostics object
                # Otherwise, create a diagnostics object for the following data.
                if diagnostics:
                    self.diagnostics.append(diagnostics)
                diagnostics = RunDiagnostics(date, time, jdir, job, fields)
            else:
                # Read data into diagnostics object
                diagnostics.read(line)

class SimulationDiagnostics(object):
    def __init__(self, files=[]):
        # Holds the diagnostics for a given set of simulation diagnostics files
        self.files = files
        self.file_diagnostics = [FileDiagnostics(f) for f in self.files]
        self.combine()

    def combine(self):
        fd = FileDiagnostics()
        for f in self.file_diagnostics:
            fd += f
        self.diagnostics = fd.combined_diagnostics()

    def plot(self):
        # Plot the diagnostics if they are present
        if self.diagnostics:
            self.diagnostics.plot()

if __name__ == "__main__":
    sdiag = SimulationDiagnostics(args.infiles)
    sdiag.plot()
