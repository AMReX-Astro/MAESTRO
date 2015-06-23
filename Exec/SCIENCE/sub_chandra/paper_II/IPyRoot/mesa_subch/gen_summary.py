#!/usr/bin/env python

#Title:         Summary Generator (gen_summary.py)
#Author:        Adam Jacobs
#Creation Date: 05/29/2015

#Usage: In the directory of the tXrhoYY.out files, just execute: ./gen_summary.py

#Description: Generate a summary of the results stored in the tXrhoYY.out files.
#

#TODO: 
# 1) 
# 2) 

import sys
import operator 
from glob import glob

def summarize_series(fglob, outfile):
   """Write out a file summarizing a series of output files from sub-Ch MESA calcs"""
   with open(outfile, mode='w') as of:
      #Iterate over files
      flist = glob(fglob)
      flist = sorted(flist)
      lgrho = [] #list of log(rho) values, parallel to the list of rxnmap maps
      rxns = []  #list of maps of form 'rxn_name' --> energy release [erg/g/s]
      for f in flist:
         rxnmap = {}
         currxn = ''
         eps_nuc = ''
         for line in open(f,mode='r'):
            if not currxn and line.count('reaction name') == 1:
               i1 = line.index('<') + 1
               i2 = line.index('>')
               currxn = line[i1:i2]
            elif currxn and line.count('eps_nuc') == 1:
               eps_nuc = float(line.partition(':')[2].strip())
               rxnmap[currxn] = eps_nuc
               currxn = ''
            elif line.count('log rho') == 1:
               lgrho.append(line.partition('rho')[2].strip())
         srtmap = sorted(rxnmap.items(), key=operator.itemgetter(1), reverse=True) #sort on values
         rxns.append(srtmap)

      #Write header
      of.write('log(rho): ' + (' {:3.3s}                  |'*len(lgrho)).format(*lgrho) + '\n')

      #Write rows of data for each logrho, include top ten rxns
      start =  '           '
      for i in range(10):
         of.write(start)
         for tup in rxns:
            of.write('{:23s}'.format(tup[i][0]))
         of.write('\n')
         of.write(start)
         for tup in rxns:
            of.write('{:<23.8e}'.format(tup[i][1]))
         of.write('\n\n')
         
if __name__== "__main__":
   summarize_series('t2*.out', '200MKSeries.out')
   summarize_series('t3*.out', '300MKSeries.out')
   
