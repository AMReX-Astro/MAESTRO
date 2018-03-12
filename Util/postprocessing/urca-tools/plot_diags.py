#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('infile', type=str,
                    help='Name of diagnostic file to plot.')
parser.add_argument('-tmin', '--tmin', type=float,
                    help='Minimum time to plot.')
args = parser.parse_args()

fin = open(args.infile, 'r')

# Throw away 4 lines
for _ in range(4):
    fin.readline()

# Get the column header
chead = fin.readline().strip().split('  ')[1:]
columns = []
for ch in chead:
    if ch.strip():
        columns.append(ch)

# Get the columns
data = {}
for c in columns:
    data[c] = []

for line in fin:
    ls = line.strip().split()
    for li, c in zip(ls, columns):
        data[c].append(float(li))
fin.close()

# Plotting function vs. time
def plot(t, y, ylabel):
    fig, ax = plt.subplots()
    ax.plot(np.log10(t), y)
    ax.set_xlabel('Log time (s)')
    ax.set_ylabel('{}'.format(ylabel))
    plt.savefig('{}.diag.{}.png'.format(args.infile, ylabel.replace(' ', '_')))
    plt.clf()

# Plot each column data vs. time
for k in data.keys():
    data[k] = np.array(data[k])

ilow = 0
if args.tmin:
    for i, t in enumerate(data['time']):
        if t >= args.tmin:
            ilow = i
            break

for c in data.keys():
    if c != 'time':
        plot(data['time'][ilow:], data[c][ilow:], c)

