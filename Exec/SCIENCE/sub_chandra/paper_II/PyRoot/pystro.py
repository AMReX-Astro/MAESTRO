#!/usr/bin/env python

#Title:         Maestro Python Module (pystro.py)
#Author:        Adam Jacobs
#Creation Date: June 26, 2015

#Usage: load as a module

#Description: This module provides various abstractions to facilitate carrying out
#             Maestro simulations.
#             For more on Maestro or to download it, see http://boxlib-codes.github.io/MAESTRO  
#

#TODO: 
# 1) 
# 2) 

#Notes on conventions, programming style:
#  1) All classes are in CamelCase with the first letter capitalized.  Class methods
#     are also in CamelCase with the first letter lower-case, e.g. myMethod().
#  2) Non-class functions and all variables are lowercase with underscores (_) 
#     acting as delimiters when needed/wanted.
#  3) An underscore prefix, e.g. _variable, means the named item is intended to be
#     private (Python doesn't provide for easily enforcing data hiding, it's all by 
#     convention).
#  4) Names that are in ALL_CAPS are intended as constants.

###########################################
### Global Imports, Data, and Constants ###
###########################################
from __future__ import print_function
import sys

###############
### Classes ###
###############
class Archiver(object):
   """A class abstracting the data and actions needed to archive Maestro output to HPSS."""
   from supercomputer import HPSS

   ##Constructor##
   def __init__(self, hpss_obj, prefix_tup=('plt','chk','inputs3d')):
      """self        --> implicitly passed reference to this instance of SCSimulation
         hpss_obj    --> A subclass of HPSS, represents the HPSS system to use
         prefix_tup  --> 3-tuple of (pltfile prefix, chkfile prefix, and inputs file prefix)"""
      if len(prefix_tup) != 3:
         raise ValueError, 'prefix_tup should have exactly three entries'
      if not isinstance(hpss_obj, HPSS):
         raise ValueError, 'hpss_obj is not an instance of HPSS!'

      self.plt_prefix    = prefix_tup[0]
      self.chk_prefix    = prefix_tup[1]
      self.inputs_prefix = prefix_tup[2]
      self.myhpss = hpss_obj
      self.pltdir = 'plotfiles'
      self.chkdir = 'checkfiles'
      self.procf  = 'processed.out'
      self.diag_interval = 4 #TODO: add this to init args

   ##Public Methods##
   def checkForFiles(self):
      """Check for new files to archive. Return True if new files found, False otherwise."""
     
      newdiag = self.checkForDiagFiles()
      newplt = self.checkForMaeFiles(self.plt_prefix, self.pltdir, 1)
      newchk = self.checkForMaeFiles(self.chk_prefix, self.chkdir, 2)

      return newdiag or newplt or newchk

   def checkForDiagFiles(self):
      """Check for new diag files to archive. Return True if new files found, False otherwise."""
      from glob import glob
      from datetime import datetime, timedelta

      #TODO: give user control of this
      DIAG_INTERVAL = timedelta(0, 4.*3600.) # be sure to archive diag files every 4 hours

      #Check if the diag file archiving interval has passed, 
      #if so signal that we have new files to archive
      dars = glob('diag_files*.tar')
      if len(dars) > 0:
         dars.sort()
         most_recent = dars[len(dars)-1]
         dar_time = self._extractTimestamp(most_recent)
         n = datetime.now()
         diag_delta = n-dar_time
         if diag_delta > DIAG_INTERVAL:
            return True
      else:
         #If there's no tar in the directory and there are diag files,
         #we need to archive
         #TODO: change this so it doesn't assume form of diag file? Let user set it?
         diags = glob('*diag.out')
         return len(diags) > 0

   def checkForMaeFiles(self, prefix, proc_dir, new_thresh=2):
      """Check for new Maestro "files."  These are actually directories of the form prefix#####, where ##### may have 5 or 6
      characters and represents the timestep of a simulation.  The directories contain a file HEADER with metadata describing
      the stored state data for either checkpointing or plotting.

      proc_dir is the directory where processed files are moved.  It is assumed to contain a file 'processed.out,' created by
      this script, that lets the script know if the file it's currently considering has already been processed.
      
      Return True if more than new_thresh new files are found, False otherwise."""

      maelist = self._getMaeFileList(prefix, proc_dir, new_thresh)
      return len(maelist) > 0

   def archiveNewFiles(self):
      """Archive newly generated files to HPSS"""

      #Archive any diag files
      if self.checkForDiagFiles():
         self._archiveDiags()

      #Archive checkpoint files
      if self.checkForMaeFiles(self.chk_prefix, self.chkdir, 2):
         self._archiveMaeFiles(self.chk_prefix, self.chkdir, 2)
      
      #Archive plot files
      if self.checkForMaeFiles(self.plt_prefix, self.pltdir, 1):
         self._archiveMaeFiles(self.plt_prefix, self.pltdir, 1)

   ##Private Methods##
   @staticmethod
   def _extractTimestamp(diag_archive):
      """Extract the timestamp from the tar file of diag files, return it as a datetime object."""
      from re import search
      from datetime import datetime

      darregex = '[0-9]{8}_([0-9]{4})'
      m  = search(darregex, diag_archive)
      ts_str = diag_archive[m.start():m.end()]
      #format of ts_str is YYYYMMDD_HHMM
      year  = int(ts_str[:4])
      month = int(ts_str[4:6]) 
      day   = int(ts_str[6:8]) 
      hour  = int(ts_str[9:11]) 
      min   = int(ts_str[11:13]) 
      ts    = datetime(year, month, day, hour, min)

      return ts

   def _archiveDiags(self):
      """Tar up diagnostic files, archive them on HPSS"""
      from subprocess import Popen, PIPE, STDOUT
      from datetime import datetime
      from glob import glob
      import shlex
      from warnings import warn

      #Generate string of diag files to tar up
      #ASSUMPTION/TODO: We assume diag files have the form *diag.out, 
      #should make this user-controlled.
      dglob = glob('*diag.out')
      if len(dglob) < 1:
         warn("You called _archiveDiags when there are no diag files to archive")
         return
      diag_files_str = ' '
      for s in dglob:
         diag_files_str = diag_files_str + s + ' '
        
      #Generate filename for the tar file
      #format of tstamp is YYYYMMDD_HHMM
      tarchive_pre = 'diag_files_'
      tstamp = datetime.today().strftime("%Y%m%d_%H%M")
      tarchive_post = '.tar'
      tarchive = tarchive_pre + tstamp + tarchive_post

      #Construct command
      command_line = 'tar cf ' + tarchive + ' ' + diag_files_str
      args = shlex.split(command_line)

      #Tar the files
      tarproc = Popen(args, stdout=PIPE, stderr=PIPE)
      (sout, serr) = tarproc.communicate()
      retcode = tarproc.poll()
      if retcode != 0:
         print 'tar stdout: ', sout
         print 'tar stderr: ', serr
         raise RuntimeError, 'tar of diag files failed!'
     
      #Archive the tar
      self.myhpss.sendToHPSS(tarchive)

   def _archiveMaeFiles(self, prefix, proc_dir, new_thresh=2):
      """Archive Maestro files prefix##### (chk and plt) to HPSS, move them to proc_dir, mark them processed.
      
      new_thresh is the number of files to leave untouched.  For example, if it's 2 then the 2 files with the 
      largest timestamp will not be archived."""
      from subprocess import Popen, PIPE, STDOUT
      from datetime import datetime
      from glob import glob
      import shlex
      from warnings import warn
      from os.path import join, isfile, isdir

      #Get a list of Maestro files to archive
      maelist = self._getMaeFileList(prefix, proc_dir, new_thresh)

      for f in maelist:
         #Archive the file
         self.myhpss.sendToHPSS(f)

         #Move it and mark as processed
         os.rename(f, join(proc_dir, f))
         proc_file = join(proc_dir, self.procf)
         print 'pf: ', proc_file
         with open(proc_file, 'a') as pf:
            print 'writing', f + '\n'
            pf.write(f + '\n')
        

   def _getMaeFileList(self, prefix, proc_dir, new_thresh=2):
      """Return a list of Maestro files (prefix#####) that are not marked as processed in proc_dir.
      Ignore the new_thresh(an integer) most recent files."""
      from re import match
      from os.path import join, isfile, isdir

      #Build set of already processed files
      proc_set = set()
      proc_file = join(proc_dir, self.procf)
      if isfile(proc_file):
         with open(proc_file) as f:
            for line in f:
               proc_set.add(line.strip())

      #Look for new plt and chk files, 
      #we want to archive if there are more than 2
      maelist = []
      for f in os.listdir('.'):
         #Is this a mae file with the right prefix?
         maeregex = '.*(' + prefix + ')([0-9]{5,6})$'
         if isdir(f):
            if match(maeregex, f) and f not in proc_set:
               if isfile(join(f,'HEADER')):  #Only add fully written files, i.e. those with HEADER
                  maelist.append(f)
      maelist.sort()
      return maelist[:-new_thresh]  #Will return empty list if len <= 2



#################
### Functions ###
#################


#################
### Execution ###
#################
#This is only for testing.  pystro.py is intended to be used as a module.
if __name__== "__main__":
   if(len(sys.argv) <= 1):
      #TODO: Add any desired testing
      pass
