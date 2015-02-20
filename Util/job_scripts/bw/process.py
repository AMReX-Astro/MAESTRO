#!/usr/bin/env python2

# Quickstart: check out the main execution region at the bottom for variables you'll want to
#  customize.  Basic usage is to execute (assuming bash)
#      python2 process.py >> process.out 2>&1 < /dev/null &
#      PID=$!
#      trap 'kill -s TERM $PID' EXIT TERM HUP XCPU KILL
#  in your job submission script.  Output will be archived until the job's time limit is reached.
#  You may also want to add
#      rm -f process.pid
#  to the end of your job script.  If there's no process.pid file, the script terminates itself.
#
#  TODO/WARNING: The script is currently only operational for Blue Waters, but is designed to be easily
#  extendable to other computing centers with whatever HPSS setup they have.
#
#
# process.py is a Python script for processing output of a Maestro run.
# It archives chk, plt, and diag files to HPSS.
#
# Please see the bottom of this script for a set of customizable variables
#
# This script made use of nhoffman's useful pyscript template:
# repo: https://gist.github.com/3006600.git 
# web:  https://gist.github.com/nhoffman/3006600
#
# ASSUMPTION: This script will be executed in a scratch directory where the simulation is running.

import os
import sys
import abc

#Global variables
pidfile_name = 'process.pid'  #Lock file to prevent multiple instances of this script
# Use this to turn on global debugging statements
#DEBUG = True
DEBUG = False

class HPSS(object):
   """An abstract class representing a high performance storage system."""
   __metaclass__ = abc.ABCMeta

   ##Constructor##
   @abc.abstractmethod
   def __init__(self, hpss_base_dir):
      """self          --> implicitly passed reference to this instance of HPSS
         hpss_base_dir --> the base HPSS directory all other files and directories will be relative to"""

      self.loc_dir = os.getcwd()
      self.loc_base = os.path.basename(self.loc_dir) 
      self.hpss_base_dir = hpss_base_dir

   ##Public Methods##
   @abc.abstractmethod
   def sendToHPSS(self, item):
      """Send self.loc_dir/item (a file or directory) 
      to the HPSS directory self.hpss_base_dir/self.loc_base/item"""
      raise NotImplementedError, "You must implement this abstract method!"

   @abc.abstractmethod
   def getFromHPSS(self, item):
      """Retrieve self.hpss_base_dir/self.loc_base/item (a file or directory) 
      and save it to self.loc_dir/item"""
      raise NotImplementedError, "You must implement this abstract method!"


   ##Private Methods##


class BWHPSS(HPSS):
   """Concrete subclass of HPSS representing an implementation of Blue Waters' HPSS system."""

   ##Constructor##
   def __init__(self, hpss_base_dir):
      """self          --> implicitly passed reference to this instance of HPSS
         hpss_base_dir --> the base HPSS directory all other files and directories 
                           will be relative to"""
      super(BWHPSS, self).__init__(hpss_base_dir)
      self.local_ep     = "ncsa#BlueWaters"
      self.storage_ep     = "ncsa#Nearline"
      self.go_host    = "cli.globusonline.org"
      self.sync_level = "3" # 0: Copy files that do not exist at the destination
                            # 1: Copy files if the destination size doesn't match
                            #      source size
                            # 2: Copy files if destination timestamp is older than
                            #      source
                            # 3: Copy files if source and destination checksums 
                            #      don't match
   
   ##Public Methods##
   def sendToHPSS(self, item):
      """Send self.loc_dir/item (a file or directory) 
      to the HPSS directory self.hpss_base_dir/self.loc_base/item"""

      #Generate a task id
      task_id = self._genID()

      #Execute the transfer
      self._transfer(task_id, item)

   
   def getFromHPSS(self, item):
      """Retrieve self.hpss_base_dir/dir/item (a file or directory) 
      and save it to self.loc_dir/item"""
      raise NotImplementedError, "You must implement this abstract method!"

   ##Private Methods##
   def _genID(self):
      """Use the Globus CLI to generate a task ID for a transfer"""
      #  NOTE: Why do we want to do this?  This is Globus' way of making transfers resilient to
      #  failure.  We transfer in a loop that keeps trying the transfer until successful.
      #  Having a unique id allows Globus to keep track of details.  If a transfer fails and you
      #  try again, Globus will notice you're attempting the same transfer and will clean up any 
      #  leftover state or meta data from the failed attempt.  If a transfer
      #  succeeds and you for some reason try again, it will do nothing and return successfully.
      from subprocess import Popen, PIPE, STDOUT
      import shlex

      TRY_MAX = 4
      tries = 0
      tid = None
      serr = ''

      while tries < TRY_MAX:
         #Create process to generate id
         #Command to generate id:   task_id=$(ssh $go_host transfer --generate-id)
         gid_command = 'ssh ' + self.go_host + ' transfer --generate-id'
         args = shlex.split(gid_command)
         if DEBUG:
            print 'executing: ', gid_command
         gid_proc = Popen(args, stdout=PIPE, stderr=PIPE)

         #Check the output
         (sout, serr) = gid_proc.communicate()
         retcode = gid_proc.poll()
         if retcode == 0:
            tid = sout.replace('\n','')
            break
         else:
            tries += 1

      if tries >= TRY_MAX:
         print '_genID() error output: ', serr
         raise RuntimeError, 'Globus failed to generate a task id! Make sure all endpoints are activated'

      return tid

   def _transfer(self, task_id, item):
      """Use the Globus CLI to transfer item using the given task"""
      from os.path import isdir, isfile
      from subprocess import Popen, PIPE, STDOUT
      import shlex
      from warnings import warn

      TRY_MAX = 4
      SUCCESS = 0
      SSH_ERROR = 255
      tries = 0
      
      #Is item a file or directory?
      dir_flag = ''
      if isdir(item):
         dir_flag = '-r'
         if not item.endswith('/'):
            item = item + r'/'  #Directories must end with a / according to Globus
      else:
         if not isfile(item):
            raise ValueError, 'The argument item is invalid! Must be local file or directory'

      #Command to transfer: ssh $go_host transfer --verify-checksum --taskid=$task_id --label=\"$task_label\" -s $sync_level   -- $src_dst [-r]
      #Build src/dst string for the transfer command
      src = self.local_ep + self.loc_dir + '/' + item
      dst = self.storage_ep + self.hpss_base_dir + '/' + self.loc_base + '/' + item
      src_dst = src + ' ' + dst
     
      #Make a label (remove any illegal characters [./])
      task_label = 'Archive ' + item
      task_label = task_label.replace('.', '-')
      task_label = task_label.replace('/', '')


      #Build the transfer command string
      trans_command = 'ssh ' + self.go_host + ' transfer --verify-checksum --taskid=' + task_id
      trans_command += r' --label=\"' + task_label + r'\" -s ' + self.sync_level 
      trans_command += ' -- ' + src_dst + ' ' + dir_flag

      #Execute the transfer
      if DEBUG:
         print 'calling: ', trans_command
      args = shlex.split(trans_command)
      while tries < TRY_MAX:
         trans_proc = Popen(args, stdout=PIPE, stderr=PIPE)
         (sout, serr) = trans_proc.communicate()
         retcode = trans_proc.poll()

         if retcode == SUCCESS:
            if DEBUG:
               print 'Transfer successfully established!'
            #Wait for the transfer to complete
            #TODO: Maybe we don't want to do this? For now I like it, makes sure transfers finish
            wait_command = 'ssh ' + self.go_host + ' wait ' + task_id
            args = shlex.split(wait_command)
            wait_proc = Popen(args)
            retcode = wait_proc.wait()
            if retcode == SUCCESS:
               return
            else:
               raise RuntimeError, 'Globus wait failed on task ' + task_id
         elif retcode == SSH_ERROR:
            tries += 1
            warn('ssh had a connection error, retry ' + str(tries) + ' of ' + str(TRY_MAX))
         else:
            print sout
            print serr
            raise RuntimeError, 'Globus transfer command failed fatally! Make sure endpoints are activated'

      if tries >= TRY_MAX:
         raise RuntimeError, 'Globus failed to transfer ' + item + ' after ' + str(tries) + ' tries!'


class Archiver(object):
   """A class abstracting the data and actions needed to archive Maestro output to HPSS."""

   ##Constructor##
   def __init__(self, hpss_obj, prefix_tup=('plt','chk')):
      """self        --> implicitly passed reference to this instance of SCSimulation
         hpss_obj    --> A subclass of HPSS, represents the HPSS system to use
         prefix_tup  --> 2-tuple of (pltfile prefix, chkfile prefix)
                         These are the common prefixes preceding the timestep shared by all
                         plt or chk files."""
      if len(prefix_tup) != 2:
         raise ValueError, 'prefix_tup should have exactly three entries'
      if not isinstance(hpss_obj, HPSS):
         raise ValueError, 'hpss_obj is not an instance of HPSS!'

      self.plt_prefix    = prefix_tup[0]
      self.chk_prefix    = prefix_tup[1]
      self.myhpss = hpss_obj
      self.pltdir = 'plotfiles'
      self.chkdir = 'checkfiles'
      self.procf  = 'processed.out'
      self.diag_interval = 4  #archive diagnostic files every 4 hours
      self.chk_interval  = 50 #archive checkpoint data for every 50 timesteps
      self.chk_thresh    = 2  #Do not archive the 2 most recent chk files
      self.plt_thresh    = 1  #Do not archive the most recent plt file

   ##Public Methods##
   def checkForFiles(self):
      """Check for new files to archive. Return True if new files found, False otherwise."""
     
      newdiag = self.checkForDiagFiles()
      newplt = self.checkForMaeFiles(self.plt_prefix, self.pltdir, self.plt_thresh)
      newchk = self.checkForMaeFiles(self.chk_prefix, self.chkdir, self.chk_thresh)

      return newdiag or newplt or newchk

   def checkForDiagFiles(self):
      """Check for new diag files to archive. Return True if new files found, False otherwise."""
      from glob import glob
      from datetime import datetime, timedelta

      DIAG_INTERVAL = timedelta(0, self.diag_interval*3600.) # be sure to archive diag files every 4 hours

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
      """Check for new Maestro "files."  These are actually directories of the form *prefix#####, where ##### may have 5 or 6
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
      if self.checkForMaeFiles(self.chk_prefix, self.chkdir, self.chk_thresh):
         self._archiveMaeFiles(self.chk_prefix, self.chkdir, self.chk_thresh, self.chk_interval)
      
      #Archive plot files
      if self.checkForMaeFiles(self.plt_prefix, self.pltdir, self.plt_thresh):
         self._archiveMaeFiles(self.plt_prefix, self.pltdir, self.plt_thresh)

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
      print 'archiving {0}'.format(tarchive)
      self.myhpss.sendToHPSS(tarchive)

   def _archiveMaeFiles(self, prefix, proc_dir, new_thresh=2, interval=1):
      """Archive Maestro files *prefix##### (chk and plt) to HPSS, move them to proc_dir, mark them processed.
      
      new_thresh is the number of files to leave untouched.  For example, if it's 2 then the 2 files with the 
      largest timestamp will not be archived.
      
      Optionally, an interval can be set.  An interval of 50 means only every 50th timestep, 
      based on the filename, will be archived to HPSS.  The other files will be moved out of the
      way to proc_dir, but will NOT be sent to HPSS.  This can help to keep storage use below
      quota."""
      from subprocess import Popen, PIPE, STDOUT
      from datetime import datetime
      from glob import glob
      import shlex
      from warnings import warn
      from os.path import join, isfile, isdir

      #Get a sorted list of Maestro files to archive, initialize processed file
      maelist = self._getMaeFileList(prefix, proc_dir, new_thresh)
      proc_file = join(proc_dir, self.procf)

      #Get the timestep of the last archived file if we're using an interval
      last_tstep = -interval - 1
      if interval > 1:
         if isfile(proc_file):
            with open(proc_file, 'r') as pf:
               lastline = ''
               for line in pf:
                  lastline = line
               if lastline.strip():
                  #The proc file wasn't empty, so extract last timestep
                  tokens = lastline.partition(prefix)
                  last_tstep = tokens[2]
                  last_tstep = int(last_tstep)

      #Go through all files, archive based on interval, 
      #move files out of the way into proc_dir
      for f in maelist:
         #Get the timestep
         tokens = f.partition(prefix)
         tstep = tokens[2]
         try:
            tstep = int(tstep)
         except ValueError:
            print 'file is in unexpected format! expected to end with timestep: ', f
            raise ValueError

         if DEBUG:
            print 'tstep, last_tstep, interval', tstep, last_tstep, interval
         if tstep - last_tstep >= interval:
            #Archive the file
            print 'archiving ', f
            self.myhpss.sendToHPSS(f)

            #Mark it as processed
            if DEBUG:
               print 'pf: ', proc_file
            with open(proc_file, 'a') as pf:
               if DEBUG:
                  print 'writing', f + '\n'
               pf.write(f + '\n')
           
            #Update the last archived timestep
            last_tstep = tstep

         #Move the file out of the way
         print f
         print proc_dir
         print 'moving {0} to {1}'.format(f, proc_dir)
         os.rename(f, join(proc_dir, f))
        
   def _getMaeFileList(self, prefix, proc_dir, new_thresh=2):
      """Return a list of Maestro files (*prefix#####) that are not marked as processed in proc_dir.
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
         if DEBUG:
            print 'checking ', f
         if isdir(f):
            if DEBUG:
               print 'is dir!'
            if match(maeregex, f) and f not in proc_set:
               if DEBUG:
                  print 'regex match!'
               if isfile(join(f,'Header')):  #Only add fully written files, i.e. those with HEADER
                  if DEBUG:
                     print 'found ', f
                  maelist.append(f)
      maelist.sort()
      return maelist[:-new_thresh]  #Will return empty list if len <= 2

def cleanup(signum, curstackframe):
   """A signal handler that does some cleanup when this process is killed."""
   #Get rid of the pidfile to indicate no instances of process.py are running
   from os.path import isfile
   print 'process.py: caught terminating signal {0}, shutting down'.format(signum)
   if isfile(pidfile_name):
      os.remove(pidfile_name)
   sys.exit(0)


if __name__ == '__main__':
   import argparse
   import time
   from datetime import datetime
   from signal import signal, SIGHUP, SIGTERM, SIGXCPU, SIGINT, SIGQUIT, SIGABRT

   #  NOTE: If this is run in a directory mysimulation/, then a file prefix##### would be
   #  archived to hpss_dir/mysimulation/prefix#####
   sleeptime = 60 #number of seconds to sleep before checking again for new files
   hpss_dir="/projects/sciteam/jni/ajacobs/sub_chandra" #base directory all archived files 
                                                        #will be sent to

   #We only want one instance of this script running in a directory
   #Use a lock file named pidfile_name to enforce this
   pidfile = None
   if os.path.isfile(pidfile_name):
      print pidfile_name, ' already exists! Make sure no instances of this script are running.'
      sys.exit(1)
   else:
      pidfile = open(pidfile_name, 'w')
      pidfile.write(str(os.getpid()))
      pidfile.close()

   print 'process.py starting up ', datetime.today().strftime("%Y-%m-%d %H:%M")

   #Setup signal handling.  This way we clean up before the script terminates
   signal(SIGHUP, cleanup); signal(SIGTERM, cleanup); signal(SIGXCPU, cleanup)
   signal(SIGINT, cleanup); signal(SIGQUIT, cleanup); signal(SIGABRT, cleanup)

   #Initialize 
   bwhpss = BWHPSS(hpss_dir)
   arch = Archiver(bwhpss)
   #Customize archiving options
   arch.pltdir = 'plotfiles'
   arch.chkdir = 'checkfiles'
   arch.procf  = 'processed.out'
   arch.diag_interval = 4  #archive diagnostic files every 4 hours
   arch.chk_interval  = 50 #archive checkpoint data for every 50 timesteps
   arch.chk_thresh    = 2  #Do not archive the 2 most recent chk files
   arch.plt_thresh    = 1  #Do not archive the most recent plt file

   #Execute the archiving algorithm repeatedly until this process is killed
   #or the lock file is removed.
   while os.path.isfile(pidfile_name):
      if arch.checkForFiles():
         #New files found, let's archive them!
         if DEBUG:
            print 'found new files!'
         arch.archiveNewFiles()
      else:
         if DEBUG:
            print 'found no new files :('
      if DEBUG:
         print 'sleep for {0}s'.format(sleeptime)
      time.sleep(sleeptime)

   print 'process.py: lock file {0} removed, shutting down...'.format(pidfile_name)

#TODO: change the script to just use 'python', not 'python2'. Add a check for version. 
#      As of now we assume python 2.6 because this is what a lot of supercomputers have by default
