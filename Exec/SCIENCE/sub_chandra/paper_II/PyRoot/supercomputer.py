#!/usr/bin/env python

#Title:         Supercomputer Python Module (supercomputer.py)
#Author:        Adam Jacobs
#Creation Date: June 26, 2015

#Usage: load as a module

#Description: This module contains code abstractions for a supercomputer.  The goal is
#             to have a set of abstractions other scripts can use.  This will allow
#             the same script to run on multiple supercomputers without modification.
#             Many of the classes, methods, etc will need to be implemented for each
#             supercomputing center one uses.

#TODO: 
# 1) 
# 2) 

###########################################
### Global Imports, Data, and Constants ###
###########################################
import abc  #AbstractBaseClass

###############
### Classes ###
###############
class TermColors(object):
   """A struct-like class containing terminal color constants.
   
   Example: red_string = TermColors.START_RED + "I'm red!" + TermColors.END """

   START_RED   = '\033[31m'
   START_GREEN = '\033[32m'
   START_BLUE  = '\033[34m'
   RESET       = '\033[0m'

class SupercomputerConfig(object):
   """A struct-like class for collecting the configuration of this supercomputer"""
   default_queue = ''    #Default PBS queue
   default_accnt = ''    #Default account for charging CPU use to
   projs_base = ''       #Base directory for projects
   scratch_base = ''     #Base scratch directory
   proc_script = ''      #Name of the processing script to be run in background of the simulation 
   run_script = ''       #Name of PBS batch run script
   npn    = -1           #NUMA nodes per physical node

class JobConfig(object):
   """A struct-like class for collecting the configuration of 
   a given job to be run on the supercomputer."""
   nodes  = -1      #Total number of nodes
   mpn    = -1      #MPI tasks per node
   omp    = -1      #OpenMP threads
   label  = ''      #Job's label
   queue  = ''      #queue for job
   accnt  = ''      #Account job is to be charged to
   time   = -1      #Wall time for job to run
   exe    = ''      #Name of executable
   inputs = ''      #Inputs file to pass to the executable

   @staticmethod
   def toPBSDict(jconfig):
      """Convert a JobConfig object into a dictionary corresponding to a PBS
      job batch script."""
      job_dict =  {
            "#PBS -A": " " + jconfig.accnt,
            "#PBS -N": " " + jconfig.label,
            "#PBS -q": " " + jconfig.queue,
            "#PBS -l": " " + "walltime={0:02d}:{1:02d}:{2:02d},nodes={3:d}".format(
               jconfig.time.hour, jconfig.time.minute, jconfig.time.second,
               jconfig.nodes),
            "# this": " script runs with {0:d} threads,"
            "{1:d} MPI tasks/node, and {2:d} nodes on titan".format(jconfig.omp,
               jconfig.mpn, jconfig.nodes),
            "export OMP_NUM_THREADS=" : "{0:d}".format(jconfig.omp),
            "aprun " : "-n {0:d} -N {1:d} -d $OMP_NUM_THREADS -ss ./{2:s} {3:s} ${{restartString}}".format(
               jconfig.nodes * jconfig.mpn, jconfig.mpn, jconfig.exe, jconfig.inputs)
            }
      return job_dict



class Supercomputer(object):
   """An abstract class representing a supercomputer generically.
   Provides an API for interacting with a supercomputer.  Details
   are implemented for any given machine."""
   __metaclass__ = abc.ABCMeta
   
   ##Constructor##
   @abc.abstractmethod
   def __init__(self, myconfig):
      """self          --> implicitly passed reference to this instance of Supercomputer
         myconfig      --> Configuration of this machine"""
      self.myconfig = myconfig

   ##Public Methods##
   @abc.abstractmethod
   def getProjsAndScratch(self):
      """Return a tuple of (projects base directory, scratch directory) for the host supercomputer."""
      raise NotImplementedError("You're calling an abstract method directly.")
         
   @abc.abstractmethod
   def generateJobScript(self, outpath, jobconfig, template_path=None):
      """Generate a PBS script to submit to the queue
      
      outpath       --> path to save the new job script in
      jobconfig     --> A supercomputer.JobConfig object for the new job
      template_path --> Full path to the template file for the job script"""
      raise NotImplementedError("You're calling an abstract method directly.")

   @abc.abstractmethod
   def isQueued(self, label):
      """Is this job in the queue?
      
      label         --> name of job you are checking for"""
      raise NotImplementedError("You're calling an abstract method directly.")

   @abc.abstractmethod
   def getQStatus(self, label):
      """Get a tuple of (queue status string, the number of times this job is in the queue)
      
      label         --> name of job you are checking for"""
      raise NotImplementedError("You're calling an abstract method directly.")

   @abc.abstractmethod
   def printUsageOverview(self):
      """Print an overview of the current user's usage including current queue status, 
      allocation's hours used, HPSS usage."""
      raise NotImplementedError("You're calling an abstract method directly.")

   @staticmethod
   def getCurrentSC(scomp_config=None):
      """Return a Supercomputer object for the current host machine."""
      import os
      import getpass
      
      host = os.environ.get('HOST')
      if host is None:
         raise OSError("No environment variable HOST.")
      elif host.count('titan') == 1:
         try:
            scmodule = __import__("titan")
            return scmodule.TitanSC(scomp_config)
         except ImportError:
            print('Cannot find module for host {}'.format(host))
            return None
      elif host.count('h2ologin') == 1:
         try:
            scmodule = __import__("bluewaters")
            return scmodule.BWSC()
         except ImportError:
            print('Cannot find module for host {}'.format(host))
            return None
      else:
         raise NotImplementedError("Unimplemented host: ", host)

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



#################
### Functions ###
#################

#################
### Execution ###
#################
#This is only for testing.  supercomputer.py is intended to be used as a module.
if __name__== "__main__":
   if(len(sys.argv) <= 1):
      #TODO: Add any desired testing
      pass

