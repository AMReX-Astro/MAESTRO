#!/usr/bin/env python

#Title:         Module to facilitate Sub-Chandra processing and analysis (subchandra.py)
#Author:        Adam Jacobs
#Creation Date: June 26, 2015

#Usage: load as a module

#Description: This module contains code utilized by the Sub-Chandra IPython notebooks and
#             Sub-Chandra processing/analysis scripts.
#             For now it's just a huge module that I hope to break up as subsets become clear.
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
class SCSimulationSuite(object):
   """A class representing a suite of sub-Chandra simulations and the different
   operations one might want to perform on such such a suite."""
   ##Global Data##

   ##Constructor##
   def __init__(self, stage_dir, scratch_dir, label):
      """self        --> implicitly passed reference to this instance of SCSimulation
         stage_dir   --> staging directory containing all simulations
         scratch_dir --> scratch directory where the work is done but data's purged
         label       --> this suite's label (e.g. SubChandraII)"""
      self._stage_dir = stage_dir.rstrip('/') #Get rid of any trailing '/'
      self._scratch_dir = scratch_dir.rstrip('/')
      self._label = label

   ##Public methods##
   def listRuns(self):
      """List all active runs in this suite."""
      from supercomputer import TermColors, Supercomputer
      from os.path import isfile, isdir, join
      from glob import glob

      active_runs = self._getActiveRuns()
      curcomp = Supercomputer.getCurrentSC()

      heading = '{0:29s}|{1:14s}|{2:14s}'.format('Label', 'In scratch?', 'In queue?')
      list_format = '{0:29s}|{1:14s}|{2:14s}'
      yep    = TermColors.START_GREEN + '{0:14s}'.format("Yes!")    + TermColors.RESET
      nope   = TermColors.START_RED   + '{0:14s}'.format("No!")     + TermColors.RESET
      purged = TermColors.START_BLUE  + '{0:14s}'.format("Purged!") + TermColors.RESET
      
      print(heading)
      for r in active_runs:
         #Check scratch: is it there, not, or there but purged?
         rundir = join(self._scratch_dir, r)
         if isdir(rundir):
            found_all_expected_files = (
                  len(glob( join(rundir, 'main.*') ))         > 0 and
                  len(glob( join(rundir, 'inputs*') ))        > 0 and
                  len(glob( join(rundir, 'helm_table.dat') )) > 0
                  )
            if found_all_expected_files:
               sc_str = yep
            else:
               sc_str = purged
         else:
            sc_str = nope
         
         #Check queue 
         #  ASSUMPTION: run directory is same as queue label
         if curcomp.isQueued(r):
            q_str = yep
         else:
            q_str = nope

         outstr = list_format.format(r, sc_str, q_str) 
         print(outstr)

   def printStatus(self, temp_tol=5.e8, mach_tol=0.3):
      """Go through all runs in the staging directory, printing out:
              run status (queued, running, or inactive)
              T_peak
              M_peak
              Simulation time"""
      print('run label:')
      print('run status | T_peak | M_peak | time | note \n')

      #Loop over all runs in staging directory assuming <run label>/[output, run, plots]
      #structure
      runs = self._getActiveRuns()
      for rr in runs:
         #Get properties
         rstat = self._getRStat(rr)
         tpeak = self._getTPeak(rr, temp_tol)
         mpeak = self._getMPeak(rr, mach_tol)
         time  = self._getTime(rr)
         note  = self._getNote(rr)

         #Print
         print(rr + ":\n" + rstat + "|" + tpeak + "|" + mpeak + "|" + time + "|" + note + '\n')

      return


   ##Private methods##
   def _getRStat(self, run_label):
      """Determine run status of run_label"""
      from supercomputer import Supercomputer
      
      curcomp = Supercomputer.getCurrentSC()
      (job_status, qcount) = curcomp.getQStatus(run_label)

      return job_status + '*{0:02d}'.format(qcount)

   def _getTPeak(self, run_label, temp_tol):
      """Return the largest temperature found in diagnostic .out files with a preference for
         files found in the work directory"""
      import os
      import numpy as np
      from supercomputer import Supercomputer, TermColors
     
      stg_dir, wrk_dir = self._stage_dir, self._scratch_dir

      #Check work directory
      if(os.path.isfile(wrk_dir + '/' + run_label + '/subchandra_temp_diag.out')):
         diag_file = wrk_dir + '/' + run_label + '/subchandra_temp_diag.out'
         temps = np.loadtxt(diag_file, usecols=(1,), unpack=True)
         maxt = max(temps)
         if(maxt > temp_tol): 
            return TermColors.START_RED + '{:6.2f}'.format(maxt/1.e6) + TermColors.RESET
         else:
            return '{:6.2f}'.format(maxt/1.e6) 
      #If no luck, check staging directory
      elif(os.path.isfile(stg_dir + '/' + run_label + '/output/subchandra_temp_diag.out')):
         diag_file = stg_dir + '/' + run_label + '/output/subchandra_temp_diag.out'
         temps = np.loadtxt(diag_file, usecols=(1,), unpack=True)
         maxt = max(temps)
         if(maxt > temp_tol): 
            return TermColors.START_RED + '{:6.2f}'.format(maxt/1.e6) + TermColors.RESET
         else:
            return '{:6.2f}'.format(maxt/1.e6) 
      else:
         return 'no files'
   
   def _getMPeak(self, run_label, mach_tol):
      """Return the largest Mach number found in diagnostic .out files with a preference for
         files found in the work directory"""
      import os
      import numpy as np
      from supercomputer import TermColors
      
      stg_dir, wrk_dir = self._stage_dir, self._scratch_dir
      
      #Check work directory
      if(os.path.isfile(wrk_dir + '/' + run_label + '/subchandra_vel_diag.out')):
         diag_file = wrk_dir + '/' + run_label + '/subchandra_vel_diag.out'
         machnums = np.loadtxt(diag_file, usecols=(3,), unpack=True)
         maxm = max(machnums)
         if(maxm > mach_tol): 
            return TermColors.START_RED + '{:5.3f}'.format(maxm) + TermColors.RESET
         else:
            return '{:5.3f}'.format(maxm)
      #If no luck, check staging directory
      elif(os.path.isfile(stg_dir + '/' + run_label + '/output/subchandra_vel_diag.out')):
         diag_file = stg_dir + '/' + run_label + '/output/subchandra_vel_diag.out'
         machnums = np.loadtxt(diag_file, usecols=(3,), unpack=True)
         maxm = max(machnums)
         if(maxm > mach_tol): 
            return TermColors.START_RED + '{:5.3f}'.format(maxm) + TermColors.RESET
         else:
            return '{:5.3f}'.format(maxm)
      else:
         return 'no files'
   
   def _getTime(self, run_label):
      """Return the latest time found in diagnostic .out files with a preference for
         files found in the work directory"""
      import os
      
      stg_dir, wrk_dir = self._stage_dir, self._scratch_dir

      #Check work directory
      diag_file = None
      if(os.path.isfile(wrk_dir + '/' + run_label + '/subchandra_temp_diag.out')):
         diag_file = wrk_dir + '/' + run_label + '/subchandra_temp_diag.out'
      #If no luck, check staging directory
      elif(os.path.isfile(stg_dir + '/' + run_label + '/output/subchandra_temp_diag.out')):
         diag_file = stg_dir + '/' + run_label + '/output/subchandra_temp_diag.out'
      else:
         return 'no files'
      
      with open(diag_file, 'r') as f:
         f.seek(-2,2)              #Jump to near the end of the file (second-to-last byte)
         while f.read(1) != '\n':  #Read to EOL
            f.seek(-2, 1)          #Jump back a bit
         last_line = f.readline()  #Read current line
      tokens = last_line.split()
      return '{:06.2f}'.format(float(tokens[0]))
   
   def _getNote(self, run_label):
      """Return the any short note found in the staging directory's out directory
         as note.txt."""
      import os
      
      stg_dir, wrk_dir = self._stage_dir, self._scratch_dir

      #Check for a note, return it if it exists
      if(os.path.isfile(stg_dir + '/' + run_label + '/output/note.txt')):
         note_file = stg_dir + '/' + run_label + '/output/note.txt'
         with open(note_file, 'r') as nf:
            for i, l in enumerate(nf):
               assert (i == 0), 'Error! note.txt should have only one line!'
               txt = l
         return txt.strip('\n')
      
      #Otherwise, return a blank
      return ' '

   def _getActiveRuns(self):
      from os import listdir
      from os.path import join
      from supercomputer import Supercomputer

      curcomp = Supercomputer.getCurrentSC()
      ret = []

      for d in listdir(self._stage_dir):
         if d == 'inactive':
            continue
         ret.append(d)
      
      return ret


class SCSimConfig(object):
   """struct-like class for containing a sub-Chandra simulation's configuration."""
   #Constants
   OUTLET   = 12
   SYMMETRY = 13

   #Sim metadata
   project_label = None
   run_label = None
   Maestro_home = None
   label_extras = []
   exe = 'main.Linux.Cray.mpi.omp.exe'
   inputs = ''

   #Initial model core parameters
   Mcore = None
   Mhe   = None
   tcore = None
   tbase = None

   #Initial model supplementary parameters
   delta = 5.e6
   rmin  = 0.e0
   rmax  = 1.1e9

   #Inputs parameters, files
   init_hse    = "<full path>/sub_chandra.M_WD*hse*"
   init_extras = "<full path>/sub_chandra.M_WD*extras*"
   an_cut         = None
   base_cutoff    = 1.e4
   sponge_cen_den = None
   sponge_st_fac  = 2.0

   #Grid properties
   coarse_res = None
   levels     = None
   drdx       = 5     #This is the standard value for the resolution jump 
                      #from fine dx to base state's dr
   octant     = None

   #Required computing resources
   node_count = None
   mpn        = None #MPI tasks / node
   omp_thds   = None

class SCSimulation(object):
   """A class representing a particular sub-Chandra simulation and the different
   operations such a simulation might perform."""
   ##Global Data##
   IGTEMP = 3.5e8 

   ##Constructor##
   def __init__(self, label, proj_label):
      """self        --> implicitly passed reference to this instance of SCSimulation
         label       --> this simulation's label (e.g. 12050-107-175-3levs)
         proj_label  --> the project/suite this simulation belongs to (e.g. SubChandraII)"""
      from supercomputer import Supercomputer
      from os.path import join, split, dirname

      curcomp = Supercomputer.getCurrentSC()
      pbase, sbase = curcomp.getProjsAndScratch()
      self._stage_dir   = join(pbase, proj_label, 'Runs', label)
      self._scratch_dir = join(sbase, 'sub_chandra', label)
      self._label = label
      self._initmod = SCInitialModel(self._stage_dir, self._label)
      self._initInputs() #Inits self._inputs
      self._out = SCOutput(self._stage_dir, self._scratch_dir, self._label, self)

   ##Public methods##
   def getLabel(self):
      """Get this simulation's label"""
      return self._label

   def getIMDict(self):
      """Get a dictionary of this simulation's initial model data, 
      including keys 'radius', 'density', 'temperature', 'pressure', 'species'"""
      return self._initmod.getDataDict()
   
   def getIMCuts(self):
      """Get a struct-like class of this simulation's initial model cutoffs"""
      return self._initmod.getCuts()
   
   def getIMLims(self):
      """Get a struct-like class of this simulation's initial model limits"""
      return self._initmod.getLims()

   def getIMConvBounds(self):
      """Get the radius boundaries of convection in the initial model"""
      return self._initmod.getConvBounds()
  
   def getStageDir(self):
      """Get the staging directory for this run"""
      return self._stage_dir

   def getDiagTemp(self):
      """Get the temperature diagnostic data"""
      return self._out._diags.getTempData()

   def getDiagMach(self):
      """Get the Mach # diagnostic data"""
      return self._out._diags.getMachData()

   def getDiagEnuc(self):
      """Get the e_nuc diagnostic data"""
      return self._out._diags.getEnucData()

   def getSliceArray(self):
      """Get an array of SCSlices for this simulation"""
      return self._out._slices

   @staticmethod
   def listRuns(project):
      """List all active runs in the given project."""
      from supercomputer import TermColors, Supercomputer
      from os.path import isfile, isdir, join
      from glob import glob

      active_runs = SCSimulation._getActiveRuns(project)
      curcomp = Supercomputer.getCurrentSC()

      heading = '{0:29s}|{1:14s}|{2:14s}'.format('Label', 'In scratch?', 'In queue?')
      list_format = '{0:29s}|{1:14s}|{2:14s}'
      yep    = TermColors.START_GREEN + '{0:14s}'.format("Yes!")    + TermColors.RESET
      nope   = TermColors.START_RED   + '{0:14s}'.format("No!")     + TermColors.RESET
      purged = TermColors.START_BLUE  + '{0:14s}'.format("Purged!") + TermColors.RESET
      
      print(heading)
      for r in active_runs:
         #Check scratch: is it there, not, or there but purged?
         rundir = join(curcomp.myconfig.scratch_base, 'sub_chandra', r)
         if isdir(rundir):
            found_all_expected_files = found_required_scratch_files(rundir)
            
            if found_all_expected_files:
               sc_str = yep
            else:
               sc_str = purged
         else:
            sc_str = nope
         
         #Check queue 
         #  ASSUMPTION: run directory is same as queue label
         if curcomp.isQueued(r):
            q_str = yep
         else:
            q_str = nope

         outstr = list_format.format(r, sc_str, q_str) 
         print(outstr)
   
   @staticmethod
   def generateSimulation(simconfig):
      """A static factory method for building a new SCSimulation"""
      from supercomputer import Supercomputer
      from shutil  import copy
      from os.path import join, basename, isfile

      #Create all the files and structures that a simulation consists of
      curcomp = Supercomputer.getCurrentSC()
      SCSimulation._generateDirStructure(simconfig, curcomp) 
      SCSimulation._generateInitModel(simconfig, curcomp) 
      SCSimulation._generateInputs(simconfig, curcomp) 
      SCSimulation._generateScripts(simconfig, curcomp) 
      
      #Copy needed files to scratch
      exe_src = join(simconfig.Maestro_home, 'Exec', 'SCIENCE', 'sub_chandra', 'build', simconfig.exe)
      exe_dst = join(curcomp.myconfig.scratch_base, 'sub_chandra', simconfig.run_label, simconfig.exe)
      copy(exe_src, exe_dst)
      
      hse_src = simconfig.init_hse
      hse_dst = join(curcomp.myconfig.scratch_base, 'sub_chandra', simconfig.run_label, basename(hse_src))
      copy(hse_src, hse_dst)
      
      inputs_src = join(curcomp.myconfig.projs_base, simconfig.project_label, 'Runs', 
            simconfig.run_label, 'run', simconfig.inputs)
      inputs_dst = join(curcomp.myconfig.scratch_base, 'sub_chandra', simconfig.run_label, simconfig.inputs)
      copy(inputs_src, inputs_dst)

      proc_src = join(curcomp.myconfig.projs_base, simconfig.project_label, 'Runs', 
            simconfig.run_label, 'run', curcomp.myconfig.proc_script)
      proc_dst = join(curcomp.myconfig.scratch_base, 'sub_chandra', simconfig.run_label, 
            curcomp.myconfig.proc_script)
      copy(proc_src, proc_dst)

      run_src = join(curcomp.myconfig.projs_base, simconfig.project_label, 'Runs', 
            simconfig.run_label, 'run', curcomp.myconfig.run_script)
      run_dst = join(curcomp.myconfig.scratch_base, 'sub_chandra', simconfig.run_label, 
            curcomp.myconfig.run_script)
      copy(run_src, run_dst)

      helmtab = join(curcomp.myconfig.scratch_base, 'sub_chandra', simconfig.run_label, 'helm_table.dat')
      if not isfile(helmtab):
         copy(join(simconfig.Maestro_home, 'Microphysics', 'EOS', 'helmeos', 'helm_table.dat'), 
            helmtab)

   def printTimestampOverview(self):
      """Print a text summary of timestamps with output data."""
      self._out.printTimestampOverview()

   def getFineResolution(self):
      """Get the physical length in cm corresponding to the 
      resolution at the finest level of refinement."""
      xmax  = float(self._inputs['prob_hi_x'].replace('d','e'))
      mlevs = int(self._inputs['max_levs'])
      nx    = int(self._inputs['n_cellx'])
      dxf = xmax/float(nx*2**(mlevs-1))
      return dxf

   #TODO: Move to supercomputer.py
   def getpltdata(self, machine, max=None):
      """Retrieve this run's pltfile data from HPSS, store in this run's scratch directory.
      'machine' specifies which system we're on (e.g. titan), max is used to set maximum
      number of files to extract."""
      from subprocess import PIPE, STDOUT, Popen
      from os.path import join, basename, isdir, isfile
      from os import listdir
      from re import match

      if machine.lower().strip() == 'titan':
         HPSS_BASE = '/proj/ast106/sub_chandra'
         PLT_RE = r'.*_plt[0-9]{5,6}[.]tar$' #<any characters>_plt<5-6 numbers>.tar<end of string>

         #Get list of files in this simulation's HPSS directory
         avail = []
         hsi_ls = Popen(['hsi', 'ls -1 {0}'.format(join(HPSS_BASE, self._label))], stdout=PIPE, stderr=STDOUT)
         out = hsi_ls.communicate()[0]
         tokens = out.split('\n')

         #Find index where directory list starts (i.e. skip the header)
         dir_idx = len(tokens) + 1 #force error if not initialized
         for i, s in enumerate(tokens):
            if s.startswith('Username'):
               dir_idx=i+1
               break

         #Build list of pltfiles
         pltfiles = []
         for s in tokens[dir_idx:]:
            if match(PLT_RE, s.strip()):
               pltfiles.append(s)

         #Build list of pltfiles already on scratch
         scdir = join(self._scratch_dir, self._label, 'plotfiles')
         exists = []
         for d in listdir(scdir):
            hfile = join(scdir, d, 'Header')
            pltdir = join(scdir, d)
            if isdir(pltdir) and isfile(hfile):
               exists.append(d)

         #Remove any pltfiles already residing on scratch
         #from list of files to extract form HPSS
         pltfiles = [pf for pf in pltfiles if exists.count(basename(pf)[:-4]) == 0] #The [:-4] removes .tar

         #Retrieve data in background
         # htar -xf <full path>.tar >> file.out 2>&1 &
         # will write "HTAR SUCCESS" to file.out
         self._titanHPSSRetrieve(pltfiles, max)

      else:
         #TODO: implement other systems
         pass

   def getTHistData(self, step):
      """Return the temperature histogram object for the given step."""
      return self._out.getTHistData(step)

   def printOutcome(self, ig_temp=IGTEMP):
      """Analyze output and print outcome details."""
      import numpy as np
      #First get the diagnostic data
      Tpeak, timepeak, rpeak, peakslice = self._out.getPeakState()
      tconv = self._out.getTurnover() 

      #Determine the outcome
      outcome = 'quasi-equilibrium'
      if Tpeak > ig_temp:
         outcome = 'ignition'
         #TODO: Add criteria for nova

      #Check peakslice
      dt = peakslice.timestamp[1] - timepeak
      if peakslice is None:
         print("Couldn't find slice data!")
         return
      if abs(dt) > 3.0:
         print("WARNING: Couldn't find a slice near (< 3 s) ignition time! dt: ", dt)

      #Determine peak radially averaged temperature, 
      #use its location as radius of "base" of convective envelope
      avgTpeak = max(peakslice.datalist[SCSlice.ITEMP])
      brho_i = np.where(peakslice.datalist[SCSlice.ITEMP] == avgTpeak)[0][0]
      avgBaseRho = peakslice.datalist[SCSlice.IRHO][brho_i]
      rBase = peakslice.datalist[SCSlice.IRADIUS][brho_i]

      print('Outcome |\n timepeak | <tconv> | (t_peak-50)/tconv | Tpeak | r_peak |\n  dt | dr | <Tpeak> | <rho_peak> | <r_base> ')
      print(outcome, '|\n', timepeak, tconv, (timepeak-50.)/tconv, Tpeak, rpeak, '\n', dt, rpeak - rBase, avgTpeak, avgBaseRho, rBase)
      
   def validateRun(self, verbose=False):
      """Do various checks to make sure this run is behaving well.
      Set verbose to True to get all possible warnings"""
      from glob import glob
      from os.path import join, isfile

      #Output check
      self._checkOutfiles(verbose)

      #Necessary files check
      if found_required_scratch_files(self._scratch_dir):
         print('SUCCESS! Found needed scratch files')
      else:
         print('WARNING! Did not find all needed scratch files')

      #Diagnostics check
      self._checkDiags(verbose)


   ##Private methods##
   def _titanHPSSRetrieve(self, pfl, max=None):
      """Assuming Titan's system configuration, retrieves files in the plotfiles list pfl
      from HPSS archival storage."""
      from subprocess import Popen, PIPE, STDOUT
      from os.path import join, isfile

      HTAR_EXE = '/opt/public/bin/htar'
      HTAR_ARGS = '-xf'
      PF_DIR = 'plotfiles'

      #For each plotfile spawn a subprocess to extract it using htar
      #The OLCF/Titan helpdesk recommends no more than two simultaneous,
      #active htar instances to make sure you aren't depriving other users
      #of fair access to shared resources.  So only have 2 going at a time.
      proc_list = []
      if not max:
         max = len(proc_list) - 1 #process all files
      else:
         max = max - 1 #Convert into 0-based

      for i, f in enumerate(pfl):
         curp = Popen([HTAR_EXE, HTAR_ARGS + f.strip()], stdout=PIPE, stderr=STDOUT, 
               cwd = join(self._scratch_dir, self._label, PF_DIR))
         print('Spawned process for ' + f + ': ' + str(curp.pid))
         proc_list.append(curp)
         if i == max:
            break
         if i % 2 == 0:
            #Wait for each to complete, printing their output
            print('\nWaiting for some processes to complete, may take a while...')
            for p in proc_list:
               out = p.communicate()[0]
               print(str(p.pid) + ' completed with')
               print('  return code: ', p.returncode)
               print('  output: ')
               print(out)
            proc_list = []

      #Wait for any remaining processes
      if proc_list:
         print('\nWaiting for some processes to complete, may take a while...')
         for p in proc_list:
            out = p.communicate()[0]
            print(str(p.pid) + ' completed with')
            print('  return code: ', p.returncode)
            print('  output: ')
            print(out)

   def _initInputs(self):
      """Initializes a map of parameters to their values for this run's inputs
      file.  It's assumed needed class member variables have been initialized."""
      from glob import glob

      #Find the file
      #  ASSUMPTION: There's only one inputs file with the prefix 'inputs' in
      #  <stage dir>/<run label>/run/
      inputs_fname = glob(self._stage_dir + '/run/inputs*')[0]
      inputs_file = open(inputs_fname)

      #Create inputs map
      self._inputs = {}

      #Go through inputs file, fill in the inputs map
      for line in inputs_file:
         tokens = line.partition('=')
         if tokens[1]: #Only do anything if a '=' was found
            key = tokens[0].strip()
            strval = tokens[2].strip()
            self._inputs[key] = strval

   def _checkOutfiles(self, verbose):
      """Check the outfiles (<run label>.o<jobid>)"""
      from glob    import glob
      from os.path import join, basename

      MIN_LINES = 250  #All output files should have more line than this,
                       #typical count for sub-Chandra is about 27000 lines for a 4 hour run

      #Check output files
      #Does the base-state physical resolution match initial model physical resolution?
      outfiles = glob(join(self._stage_dir, 'output', '{}.o*'.format(self._label)))
      outfiles += glob(join(self._scratch_dir, '{}.o*'.format(self._label)))
      outfiles.sort()
      if len(outfiles) == 0:
         if verbose:
            print('NOTE: No output files yet')
      lowcount = 0
      total = 0
      for fname in outfiles:
         if fname.endswith('~'):
            #Skip over vim junk
            continue
         with open(fname, 'r') as f:
            total += 1
            #Check linecount
            num_lines = sum(1 for line in f)
            if num_lines < MIN_LINES:
               lowcount += 1
               if verbose: 
                  print('WARNING! Low line count for {}: {}'.format(basename(fname), num_lines))
           
            #Check outfile data
            f.seek(0, 0) #Reset to beginning of file
            base_res = None
            im_res = None
            res_warned = False
            for line in f:
               #Look for base-state, initial model resolution, make sure they match
               if line.count('dr of MAESTRO base state') == 1:
                  base_res = float(line.partition('=')[2].strip())
               if line.count('dr of input file data') == 1:
                  im_res = float(line.partition('=')[2].strip())
               if base_res is not None and im_res is not None:
                  if not res_warned:
                     if base_res != im_res:
                        print('WARNING! Base state and initial model resolution do not match')
                        print('   Base res: {}'.format(base_res))
                        print('   IM res:   {}'.format(im_res))
                        print('   file:     {}'.format(basename(fname)))
                        res_warned = True
                     else:
                        print('SUCCESS! Base state/initial model resolutions match')
                        print('   file:     {}'.format(basename(fname)))
                        res_warned = True
      if lowcount > 0:
         if not verbose:
            print('WARNING! {}/{} outfiles had low line count'.format(lowcount, total))

   def _checkDiags(self, verbose):
      """Check the status of the diagnostics files"""
      from os.path import join, isfile, basename

      #Get lists of the diagnostics files
      sc_diag = []
      st_diag = []
      sc_diag_lc = [0, 0, 0]
      st_diag_lc = [0, 0, 0]
      sc_diag.append(join(self._scratch_dir,         'subchandra_temp_diag.out'))
      sc_diag.append(join(self._scratch_dir,         'subchandra_vel_diag.out'))
      sc_diag.append(join(self._scratch_dir,         'subchandra_enuc_diag.out'))
      st_diag.append(join(self._stage_dir, 'output', 'subchandra_temp_diag.out'))
      st_diag.append(join(self._stage_dir, 'output', 'subchandra_vel_diag.out'))
      st_diag.append(join(self._stage_dir, 'output', 'subchandra_enuc_diag.out'))

      #Count stage files
      st_count = 0
      for i, f in enumerate(st_diag):
         if isfile(f):
            st_count += 1
            st_diag_lc[i] = sum(1 for line in open(f, 'r'))
         else:
            if verbose:
               print('WARNING! Missing diag file in stage: {}'.format(f))

      #Count scratch files
      sc_count = 0
      for i, f in enumerate(sc_diag):
         if isfile(f):
            sc_count += 1
            sc_diag_lc[i] = sum(1 for line in open(f, 'r'))
         else:
            if verbose:
               print('WARNING! Missing diag file in scratch: {}'.format(f))
     
      #If both are 0, we probably haven't started running yet
      if st_count == 0 and sc_count == 0:
         if verbose:
            print('NOTE: Looks like run has not started, no diag files')

      #They don't match, not good
      if st_count != sc_count:
         print('WARNING! Stage and scratch do not have the same diag file count')

      #Do the line counts match?
      for i in range(1,3):
         if sc_diag_lc[0] != sc_diag_lc[i]:
            print("WARNING! line counts don't match")
            print("   <scratch>/{}:{}".format(basename(sc_diag[0]), sc_diag_lc[0]))
            print("   <scratch>/{}:{}".format(basename(sc_diag[i]), sc_diag_lc[i]))
            return
         if st_diag_lc[0] != st_diag_lc[i]:
            print("WARNING! line counts don't match")
            print("   <stage>/{}:{}".format(basename(st_diag[0]), st_diag_lc[0]))
            print("   <stage>/{}:{}".format(basename(st_diag[i]), st_diag_lc[i]))
            return

      #Time to archive?
      if sc_diag_lc[0] != st_diag_lc[0]:
         if verbose:
            print('NOTE: You need to archive scratch diagnostics to stage')

      #All's good if all diag files found in only scratch or in both scratch and stage
      if st_count == 3 and (sc_count==3 or sc_count==0):
         print('SUCCESS! All expected diagnostic files found')
         
 
   @staticmethod
   def _generateDirStructure(simconfig, curcomp):
      """Create a directory convention, including a run directory in backed-up home and 
      a scratch directory in the purged space of the filesystem.
      
      Run directory will be <home>/Projects/simconfig.project_label/Runs/<run_label>/ containing run, output, and plots.
      Scratch directory will be <scratch>/sub_chandra/<run_label>"""
      from os import makedirs, error
      from os.path import join
      import os
      from supercomputer import Supercomputer
      import math

      #Construct <run_label> string
      run_label = "{0:02d}{1:03d}-{2:03d}-{3:03d}-{4:1d}lev".format(int(simconfig.Mcore*10), 
            int(simconfig.Mhe*1000), int(simconfig.tcore/1.e6), int(simconfig.tbase/1.e6),
            simconfig.levels)
      if simconfig.octant == False:
         run_label += "-full"
      for note in simconfig.label_extras:
         run_label += "-" + note
      simconfig.run_label = run_label

      #Create Run and scratch directory
      (projs_base, scratch_base) = curcomp.getProjsAndScratch()
      proj_base = join(projs_base, simconfig.project_label)
      run_base = join(proj_base, 'Runs', run_label)
      scratch_base = join(scratch_base, 'sub_chandra')
      scratch = join(scratch_base, run_label)
      try:
         makedirs( join(run_base, 'run') )
         makedirs( join(run_base, 'output') )
         makedirs( join(run_base, 'plots') )
      except os.error:
         print('WARNING, already exists: {}'.format(run_base))
      try:
         makedirs(scratch)
      except os.error:
         print('WARNING, already exists: {}'.format(scratch))

   @staticmethod
   def _generateInitModel(simconfig, curcomp):
      """Generate 1D initial model."""
      from shutil  import move
      from os.path import join
      from glob import glob

      #Generate an initial parameters file,
      #this one will have an overly large radius, which we'll adjust
      proj_base = join(curcomp.myconfig.projs_base, simconfig.project_label)
      pfilename = SCSimulation._generateParams(simconfig, proj_base)
      
      #Build the initial model
      SCSimulation._buildIM(simconfig, pfilename)
      
      #Repeat the above, adjusting the radius down
      (simconfig.rmax, simconfig.an_cut, simconfig.sponge_cen_den) = SCSimulation._computeRmaxAndCutoffs(simconfig, 0.5)
      SCSimulation._generateParams(simconfig, proj_base)
      SCSimulation._buildIM(simconfig, pfilename)
      
      #Move the parameters and model file to run directory
      move(pfilename, join(proj_base, 'Runs', simconfig.run_label, 'run', pfilename))
      (hse, extras) = glob('sub_chandra.M_WD*')
      simconfig.init_hse    = join(proj_base, 'Runs', simconfig.run_label, 'run', hse)
      simconfig.init_extras = join(proj_base, 'Runs', simconfig.run_label, 'run', extras)
      move(hse, simconfig.init_hse)
      move(extras, simconfig.init_extras)

   @staticmethod
   def _generateParams(simconfig, proj_base):
      """Generate parameters input file for initial model builder."""
      from os.path import join
      from shutil  import move

      #Build param dict
      param_dict = {}
      if simconfig.octant:
         param_dict['nx']     = str(simconfig.coarse_res * 2**(simconfig.levels-1) * simconfig.drdx)
      else:
         #For full star, only half of the coarse res is used for the radial base state
         param_dict['nx']     = str(simconfig.coarse_res/2 * 2**(simconfig.levels-1) * simconfig.drdx)
      param_dict['M_tot']     = str(simconfig.Mcore)
      param_dict['M_He']      = str(simconfig.Mhe)
      param_dict['delta']     = str(simconfig.delta)
      param_dict['xmin']      = str(simconfig.rmin)
      param_dict['xmax']      = str(simconfig.rmax)
      param_dict['temp_core'] = str(simconfig.tcore)
      param_dict['temp_base'] = str(simconfig.tbase)

      #Build from template file
      tfilename = join(proj_base, 'Templates', '_params.template')
      pfilename = '_params.{}'.format(simconfig.run_label) 
      SCSimulation._writeKeyValFile(tfilename, pfilename, param_dict)

      return pfilename
        
   @staticmethod
   def _writeKeyValFile(template_file, outfile, kvdict):
      """Read from a template file, use the given dictionary to write key, value pairs
      in a parameter file containing pairs of the form "key = value"."""

      with open(outfile, 'w') as ofile, open(template_file, 'r') as tfile:
         for line in tfile:
            if line.find('=') != -1:
               tokens = line.partition('=')
               key = tokens[0]
               val = tokens[2].strip()
               if key.strip() in kvdict:
                  val = kvdict[key.strip()]
               ofile.write(key.rstrip() + ' = ' + val.strip() + '\n')
            else:
               ofile.write(line)

      return outfile

   @staticmethod
   def _buildIM(simconfig, pfilename):
      """Build the initial model data."""
      from subprocess import call, Popen, PIPE, STDOUT
      from os.path import join, isfile
      from os import remove
      from glob import glob
      from shlex import split

      #Make sure helmtable is linked
      if not isfile('helm_table.dat'):
         call(['ln', '-s', join(simconfig.Maestro_home, 'Microphysics', 'EOS', 'helmeos', 'helm_table.dat')])

      #Build the executable command
      init1d_exe = join(simconfig.Maestro_home, 'Util', 'initial_models', 'sub_chandra', 
            'init_1d.Linux.gfortran.debug.exe') + ' ' + pfilename

      #Execute, removing any old IM data files
      old_files = glob('sub_chandra.M_WD*')
      for f in old_files:
         remove(f)
      i1d_proc = Popen(split(init1d_exe), stdout=PIPE, stderr=PIPE)
      (i1d_out, i1d_err) = i1d_proc.communicate()
      if i1d_err:
         print('init1d error: ', i1d_err)
         print('init1d out, last 5 lines:')
         for line in i1d_out.split('\n')[-5:]:
            print(line)
   
   @staticmethod
   def _computeRmaxAndCutoffs(simconfig, rfac):
      """Read in the initial model data, return adjusted radius to be r_peak + r_peak*rfac,
      where r_peak is the radius of peak temperature. Also return the anelastic and
      sponge central density cutoffs based on top of convective zone."""
      from numpy import loadtxt
      from glob import glob

      #Read in data
      #ASSUMPTION: initial model data for exactly one model 
      #  exists in the current directory
      imfile = glob('sub_chandra.M_WD*.hse*')[0]
      rad, rho, temp = loadtxt(imfile, usecols=(0,1,2), unpack=True)
      rmax = 0
      ancut = 0.0
      te_old = 0.0
      for r, d, t in zip(rad, rho, temp):
         dt = t - te_old
         if dt < 0.0:
            rmax = r
         if rmax and dt == 0.0:
            ancut = d_old/simconfig.sponge_st_fac
            scd   = ancut
            break
         te_old = t
         d_old = d

      rmax = rmax + rmax*rfac
      return (round(rmax, -7), round(ancut, -3), round(scd, -3))

   @staticmethod
   def _generateInputs(simconfig, curcomp):
      """Generate inputs file for the Maestro executable."""
      from os.path import join, basename
      from shutil  import move
      from numpy   import log10

      #Build inputs dict
      inputs_dict = {}

      inputs_dict['model_file'] = '"' + basename(simconfig.init_hse) + '"'
      simconfig.job_name =  '"{0:d}^3 base grid, T_core = 10^{1:d}, '.format(simconfig.coarse_res, int(log10(simconfig.tcore)))
      simconfig.job_name += 'T_base = {0:3d} MK -- M_WD={1:3.1f}, M_He={2:4.2f}"'.format(int(simconfig.tbase/1.e6), 
            simconfig.Mcore, simconfig.Mhe)
      inputs_dict['job_name']   = simconfig.job_name
      inputs_dict['max_levs']   = str(simconfig.levels)
      inputs_dict['n_cellx']    = str(simconfig.coarse_res)
      inputs_dict['n_celly']    = str(simconfig.coarse_res)
      inputs_dict['n_cellz']    = str(simconfig.coarse_res)
      inputs_dict['anelastic_cutoff']      = str(simconfig.an_cut)
      inputs_dict['base_cutoff_density']   = str(simconfig.base_cutoff)
      inputs_dict['sponge_center_density'] = str(simconfig.sponge_cen_den)
      inputs_dict['sponge_start_factor']   = str(simconfig.sponge_st_fac)
      if simconfig.octant == True:
         inputs_dict['octant'] = '.true.'
         inputs_dict['prob_hi_x'] = str(simconfig.rmax)
         inputs_dict['prob_hi_y'] = str(simconfig.rmax)
         inputs_dict['prob_hi_z'] = str(simconfig.rmax)
         inputs_dict['bcx_lo']    = str(simconfig.SYMMETRY)
         inputs_dict['bcx_hi']    = str(simconfig.OUTLET)
         inputs_dict['bcy_lo']    = str(simconfig.SYMMETRY)
         inputs_dict['bcy_hi']    = str(simconfig.OUTLET)
         inputs_dict['bcz_lo']    = str(simconfig.SYMMETRY)
         inputs_dict['bcz_hi']    = str(simconfig.OUTLET)
      elif simconfig.octant == False:
         inputs_dict['octant'] = '.false.'
         inputs_dict['prob_hi_x'] = str(simconfig.rmax*2)
         inputs_dict['prob_hi_y'] = str(simconfig.rmax*2)
         inputs_dict['prob_hi_z'] = str(simconfig.rmax*2)
         inputs_dict['bcx_lo']    = str(simconfig.OUTLET)
         inputs_dict['bcx_hi']    = str(simconfig.OUTLET)
         inputs_dict['bcy_lo']    = str(simconfig.OUTLET)
         inputs_dict['bcy_hi']    = str(simconfig.OUTLET)
         inputs_dict['bcz_lo']    = str(simconfig.OUTLET)
         inputs_dict['bcz_hi']    = str(simconfig.OUTLET)
      else:
         raise ValueError('Invalid value for octant in SCSimConfig')
      inputs_dict['plot_base_name']  = '"' + simconfig.run_label + '_plt"'
      inputs_dict['check_base_name'] = '"' + simconfig.run_label + '_chk"'

      #Build from template file
      proj_base = join(curcomp.myconfig.projs_base, simconfig.project_label)
      tfilename = join(proj_base, 'Templates', 'inputs3d.template')
      ifilename = 'inputs3d.{}'.format(simconfig.run_label) 
      SCSimulation._writeKeyValFile(tfilename, ifilename, inputs_dict)

      #Move to run directory 
      ipathname = join(proj_base, 'Runs', simconfig.run_label, 'run', ifilename)
      move(ifilename, ipathname)
      simconfig.inputs = ifilename 

   @staticmethod
   def _generateScripts(simconfig, curcomp):
      """Generate the scripts for running the job in scratch."""
      from supercomputer import JobConfig
      from datetime import time
      from os.path import join
      from shutil  import copy
      #Configure job
      jconfig = JobConfig()
      jconfig.nodes  = simconfig.node_count
      jconfig.mpn    = simconfig.mpn
      jconfig.omp    = simconfig.omp_thds
      jconfig.label  = simconfig.run_label
      jconfig.time   = time(4,0,0)         
      jconfig.exe    = simconfig.exe
      jconfig.inputs = simconfig.inputs

      #Make strings
      outpath = join(curcomp.myconfig.projs_base, simconfig.project_label, 'Runs',
            simconfig.run_label, 'run')
      tpath = join(curcomp.myconfig.projs_base, simconfig.project_label, 'Templates',
            'titan.run.template')

      #Generate job script
      curcomp.generateJobScript(outpath, jconfig, tpath)

      #Copy process script over
      proc_temp = join(curcomp.myconfig.projs_base, simconfig.project_label, 'Templates',
            'process.titan.template')
      proc_dst = join(curcomp.myconfig.projs_base, simconfig.project_label, 'Runs',
            simconfig.run_label, 'run', 'process.titan')
      copy(proc_temp, proc_dst)

   @staticmethod
   def _getActiveRuns(project):
      from os import listdir
      from os.path import join
      from supercomputer import Supercomputer

      curcomp = Supercomputer.getCurrentSC()
      proj_root = join(curcomp.myconfig.projs_base, project, 'Runs')
      ret = []

      for d in listdir(proj_root):
         if d == 'inactive':
            continue
         ret.append(d)
      
      return ret


class SCInitialModel(object):
   """A class representing a Maestro initial model for the Sub-Chandra
   problem."""
   #TODO: Write Maestro initial model super class
   ##Data##
   #Static variables 
   _STAT_VAR = 1

   #Global constants
   #IRADIUS = 0

   #Public variables 
   nspec = 0

   #Consructor
   def __init__(self, stage_dir, label):
      """SCInitialModel constructor.
         self       --> reference to this instance of SCInitialModel (implicitly passed)
         stage_dir  --> directory this run will be staged in
         label      --> run label (e.g. 12050-107-175-3lev)"""
      from glob import glob

      #Initialize simple private variables
      self._label = label
      self._stage_dir = stage_dir.rstrip('/')

      #Initialize parameters and model data 
      self._parameters = {'nx': None,
          'M_tot': None,
          'M_He': None,
          'delta': None,
          'xmin': None,
          'xmax': None,
          'temp_core': None,        
          'temp_base': None,
          'mixed_co_wd': None,
          'low_density_cutoff': None,
          'temp_fluff': None,
          'smallt': None}
      self._model_data = {'radius': None,
          'density': None,
          'temperature': None,
          'pressure': None,
          'species': None}

      #Build filenames
      #ASSUMPTIONS:
      #  1) All run information can be found in stage_dir/label/run/
      #  2) params files are prefixed with '_params' and there's only one in
      #     stage_dir/label/run/
      #  2) data files are prefixed with 'sub_chandra.M_WD' and there's only one set
      #     of 'hse' and possibly 'extras' in stage_dir/label/run/
      params_file = glob(self._stage_dir + '/run/_params*')[0]
      data_file =   glob(self._stage_dir + '/run/sub_chandra.M_WD*hse*')[0]

      #Read file data
      self._read_params(params_file)
      self._read_data(data_file, extras=True)

      #Initialize plotting limits and convection zone boundaries
      self.initLimits()
      self.initConv()

   ##Public methods##
   def __str__(self):
      """String representation of an instance of SCInitialModel.  Used when passed to print()."""
      ret = "Label: " + self._label + "\n"
      ret += "nspec: " + str(self.nspec) + "\n"
      ret += "Parameters: \n"
      for key in self._parameters:
         ret += "  " + key + " --> " + self._parameters[key] + "\n"
      ret += "Data: \n"
      for key in self._model_data:
         ret += "  " + key + " --> " + str(self._model_data[key]) + "\n"
      ret += "\n\n"
      return ret

   def initLimits(self):
      """Initializes plotting limits based on parameters in 
      stage_dir/label/run/inputs*."""
      from glob import glob

      #Initialize Limits and Cutoffs struct-like objects
      self._mylim = _Limits()
      self._mycut = _Cutoffs()

      #Build inputs filename
      #ASSUMPTIONS: 
      # 1) Inputs file is in self._stage_dir/self._label/run/
      # 2) Inputs file is prefixed with 'inputs3d' and there's only one
      inputs_file = glob(self._stage_dir + '/run/inputs3d*')[0]
  
      #Determine axis limits
      #  Temperature
      #  For now fix it at 1e7 K to 5e8 K
      self._mylim.Tlims=(1.e7, 5.e8)
  
      #  Density
      #  For now, fix it at 1 to 1e9
      self._mylim.dlims = (1.0, 1.e9)
  
      #  Radius
      #  Go from 0 to xmax in parameters file
      rmax = float(self._parameters['xmax'].replace('d', 'e'))
      self._mylim.rlims = (0., rmax)
  
      #  Soundspeed
      self._mylim.clims = (1.e7, 1.e10)

      #  Entropy
      self._mylim.slims = (1.e7, 1.e10)
      #self._mylim.slims = (0, 1)
  
      #Determine zoom bounds
      rtop, dtop, Ttop = self._get_top()
      self._mylim.rzoom = (rtop - rtop*0.05, rtop + rtop*0.15)
      self._mylim.dzoom = (dtop - dtop*0.25, dtop + dtop*1.5)
      self._mylim.Tzoom = (Ttop - Ttop*0.25, Ttop + Ttop*1.5)
      self._mylim.zbounds = (self._mylim.rzoom, self._mylim.dzoom, self._mylim.Tzoom)
  
      #Read in cutoffs
      self._mycut.an_cut = float(get_param('anelastic_cutoff', inputs_file).replace('d', 'e'))
      self._mycut.sp_cen_den = float(get_param('sponge_center_density', inputs_file).replace('d', 'e'))
      self._mycut.sp_st_fac = float(get_param('sponge_start_factor', inputs_file).replace('d', 'e'))
      self._mycut.base_cut = float(get_param('base_cutoff_density', inputs_file).replace('d', 'e'))

      #Find cutoff radii
      # find the r coordinate of the anelastic cutoff
      for r, rho in zip(self._model_data['radius'], self._model_data['density']):
         #TODO add error checking
         if rho <= self._mycut.an_cut:
            self._mycut.r_an = r
            break

      # find the r coordinate of the start of the sponge
      sp_st_rho = self._mycut.sp_cen_den * self._mycut.sp_st_fac
      for r, rho in zip(self._model_data['radius'], self._model_data['density']):
         #TODO add error checking
         if rho <= sp_st_rho:
            self._mycut.r_sp = r
            break
  
      # find the r coordinate of the cutoff density
      for r, rho in zip(self._model_data['radius'], self._model_data['density']):
         #TODO add error checking
         if rho <= self._mycut.base_cut:
            self._mycut.r_bc = r
            break

      return

   def initConv(self):
      """Initializes the lower and upper radial bounds of the convection zone
      as determined by the entropy."""
      from glob import glob
      #The DS critical values are based on analyzing models at different extremes and
      #are found to accurately demark the convective zone where ds/dr <= 0.
      #TODO: I'm noticing that these critical values depend on resolution of the initial model,
      #      which makes sense as higher resolution will have smaller ds between cells.  I should
      #      normalize by dr.
      DSDR_TRIGGER = 5.0 #This triggers the search for the convective zone, 
      DS_TRIGGER = 1.0e6 #This triggers the search for the convective zone, 
                         #it must be after this point
      DSDR_THRESH = 0.4 #When DS falls below this, we're in the convectively unstable zone, 
      DS_THRESH = 7.5e4  #When DS falls below this, we're in the convectively unstable zone, 
                         #when it rises back above it we're out of the convective zone

      ###OLD 
      #S_TOL = 1.e-16
      #normalize to avoid differences between huge numbers
      #s_norm = self._model_data['entropy']/max(self._model_data['entropy'])
      ###END OLD 
  
      ent = self._model_data['entropy']

      #Find radial bounds in which entropy is flat (i.e. isentropic)
      #  ASSUMPTION: dr is fixed
      sold = ent[1]
      dr = self._model_data['radius'][2]-self._model_data['radius'][1]
      left = None
      right = None
      trigger = False
      #print('trigger: ', DS_TRIGGER)
      #print('dr trigger: ', DSDR_TRIGGER)
      #print('thresh:  ', DS_THRESH)
      #print('dr thresh:  ', DSDR_THRESH)
      for r, s in zip(self._model_data['radius'][2:], ent[2:]):
         ds = (s-sold)
         dsdr = ds/dr
         #print('ds: ', ds)
         #print('ds/dr: ', dsdr)
         if dsdr > DSDR_TRIGGER and r > 1.e8:
            trigger = True
            #print('trigger!')
         #TODO add error checking
         if trigger:
            if not left:
               if dsdr < DSDR_THRESH:
                  left = r
                  #print('found left!')
            else:
               if not right:
                  if dsdr >= DSDR_THRESH:
                     right = r
                     #print('found right!')
                     break
         sold = s

      #Make sure we got them both
      assert left  != None, "Found no lower convective bound!  Check initConv()"
      assert right != None, "Found no upper convective bound!  Check initConv()"

      self._conv_bounds = (left, right)

   def getDataDict(self):
      """Get a dictionary of this initial model's data, 
      including keys 'radius', 'density', 'temperature', 'pressure', 'species'"""
      return self._model_data

   def getCuts(self):
      """Get a struct-like class of this initial model's cutoffs"""
      return self._mycut
   
   def getLims(self):
      """Get a struct-like class of this initial model's limits"""
      return self._mylim

   def getConvBounds(self):
      """Get a tuple of the lower and upper radius bounds of convection."""
      return self._conv_bounds


   ##Private methods##
   def _read_params(self, params_file):
      """Reads in parameters from params_file and stores them in _parameters dictionary."""
      pfile = open(params_file)
      for line in pfile:
         if(line.find('=') > -1):
            tokens = line.partition('=')
            cur_param = tokens[0].strip()
            cur_val = tokens[2]
            if(cur_param in self._parameters):
               self._parameters[cur_param] = cur_val
      pfile.close()

   def _read_data(self, data_file, extras=False):
      """Reads in data from data_file and stores them in _model_data.  If
      extras=True then soundspeend data is loaded from data_file with 'hse'
      replaced with 'extras'."""
      import numpy as np
      #ASSUMPTION!: Data file has a header of this form:
        # npts = <some integer>
        # num of variables = <some integer>
        # <variable label 1>
        # ...
        # <variable label n>

      #The number of non-species variables (radius, density, temperature, pressure)
      NONSPEC = 4

      #Parse header
      dfile = open(data_file)
      if extras:
         efile = open(data_file.replace('hse', 'extras'))
      npts = -1
      varc = -1
      varcol = {'radius': 0}
      evarcol = {}
      i = 1
      for line in dfile:
         #If no # prefix, we've finished reading the header
         if(not line.strip().startswith('#')):
            break
         if(line.find('=') > -1): 
            #Strip comment characters
            sline = line.lstrip('#')
            #Extract npts and variable count
            tokens = sline.partition('=')
            if(tokens[0].strip() == 'npts'):
               npts = int(tokens[2])
            if(tokens[0].strip() == 'num of variables'):
               varc = int(tokens[2])
         else:
            #Store column number of variable
            tokens = line.partition(' ')
            varcol[tokens[2].strip()] = i
            i = i +1
      self.nspec = len(varcol) - NONSPEC

      #If needed, parse extras header
      #ASSUMPTION: extras file has same npts as hse data file
      i=1 #i=0 is still radius
      if extras:
         for line in efile:
            #If no # prefix, we've finished reading the header
            if(not line.strip().startswith('#')):
               break
            if(line.find('=') > -1): 
               continue
            else:
               #Store column number of variable
               tokens = line.partition(' ')
               evarcol[tokens[2].strip()] = i
               i = i +1

      #Use header info to build non-species _model_data dict
      for key in self._model_data:
         if key != 'species':
            data_arr = np.loadtxt(data_file, usecols=(varcol[key],))
            self._model_data[key] = data_arr
            del varcol[key]

      #Build species part of _model_data dict
      species_dict = {}
      for key in varcol:
         data_arr = np.loadtxt(data_file, usecols=(varcol[key],))
         species_dict[key] = data_arr
      self._model_data['species'] = species_dict

      #Build extras part of _model_data dict
      if extras:
         for key in evarcol:
            data_arr = np.loadtxt(data_file.replace('hse','extras'), usecols=(evarcol[key],))
            self._model_data[key] = data_arr

   def _get_top(self):
      """Return tuple of (radius, density, temperature) values at the top of the convective zone
      (where temperature levels out)."""
      #Get data arrays
      r = self._model_data['radius']
      rho = self._model_data['density']
      temp = self._model_data['temperature']

      #Search from end of temp array until the temp changes, this is the top of
      #the convective zone.  return values at this point.
      te_prev = temp[len(temp)-1]
      r_top = -1.0
      rho_top = -1.0
      T_top = -1.0
      for re, rhoe, te in zip(r[::-1], rho[::-1], temp[::-1]):
         if(te_prev != te):
            r_top = re
            rho_top = rhoe
            T_top = te
            break
      return (r_top, rho_top, T_top)    


class SCDiagnostics(object):
   """A class representing the diagnostics printed out every timestep: peak temperature
   details, peak velocity/Mach # details, and peak nuclear burning energy details."""
   ##Shared class data##

   ##Constructor##
   def __init__(self, stage_dir, scratch_dir, label):
      """self        --> implicitly passed reference to this instance of SCSimulation
         stage_dir   --> staging directory containing all simulations
         scratch_dir --> scratch directory where the work is done but data's purged
         label       --> this simulation's label (e.g. 12050-107-175-3levs)"""
      from glob import glob
      from os.path import join, getmtime, isfile
      from warnings import warn
      import numpy as np

      #Constants
      TIME_COL = 0

      MAXT_COL = 1
      MAXT_X_COL = 2
      MAXT_Y_COL = 3
      MAXT_Z_COL = 4
      MAXT_R_COL = 8
      MAXT_VX_COL = 5
      MAXT_VY_COL = 6
      MAXT_VZ_COL = 7
      MAXT_VR_COL = 9

      MAXVMAG_COL = 1
      MAXMACH_SUB_COL = 2
      MAXMACH_FULL_COL = 3
      DT_COL = 4

      MAXENUC_COL = 1
      MAXENUC_X_COL = 2
      MAXENUC_Y_COL = 3
      MAXENUC_Z_COL = 4
      MAXENUC_R_COL = 8
      MAXENUC_VX_COL = 5
      MAXENUC_VY_COL = 6
      MAXENUC_VZ_COL = 7
      MAXENUC_VR_COL = 9

      #Store args
      self._stage_dir = stage_dir.rstrip('/') #Get rid of any trailing '/'
      self._scratch_dir = scratch_dir.rstrip('/')
      self._label = label

      #Find the most up-to-date diagnostics data files
      temp_stfname = join(self._stage_dir, 'output', 'subchandra_temp_diag.out')
      temp_scfname = join(self._scratch_dir, 'subchandra_temp_diag.out')
      if not isfile(temp_scfname):
         if not isfile(temp_stfname):
            warn('No temperature diagnostic file found for {0}!'.format(self._label))
            temp_fname = None
         else:
            temp_fname = temp_stfname
      elif getmtime(temp_scfname) >= getmtime(temp_stfname):
         temp_fname = temp_scfname
      else:
         temp_fname = temp_stfname

      vel_stfname  = join(self._stage_dir, 'output', 'subchandra_vel_diag.out')
      vel_scfname = join(self._scratch_dir, 'subchandra_vel_diag.out')
      if not isfile(vel_scfname):
         if not isfile(vel_stfname):
            warn('No velocity diagnostic file found for {0}!'.format(self._label))
            vel_fname = None
         else:
            vel_fname = vel_stfname
      elif getmtime(vel_scfname) >= getmtime(vel_stfname):
         vel_fname = vel_scfname
      else:
         vel_fname = vel_stfname

      enuc_stfname  = join(self._stage_dir, 'output', 'subchandra_enuc_diag.out')
      enuc_scfname = join(self._scratch_dir, 'subchandra_enuc_diag.out')
      if not isfile(enuc_scfname):
         if not isfile(enuc_stfname):
            warn('No enuc diagnostic file found for {0}!'.format(self._label))
            enuc_fname = None
         else:
            enuc_fname = enuc_stfname
      elif getmtime(enuc_scfname) >= getmtime(enuc_stfname):
         enuc_fname = enuc_scfname
      else:
         enuc_fname = enuc_stfname

      #Load most recent diagnostics data
      #  max temperature and time (ASSUMPTION: All diag files have same time data)
      if temp_fname != None:
         (self._time, self._maxT, self._maxT_x, self._maxT_y, self._maxT_z, self._maxT_r, 
          self._maxT_vx, self._maxT_vy, self._maxT_vz, self._maxT_vr) = np.loadtxt(temp_fname, 
          usecols=(TIME_COL, MAXT_COL, MAXT_X_COL, MAXT_Y_COL, MAXT_Z_COL, MAXT_R_COL, 
          MAXT_VX_COL, MAXT_VY_COL, MAXT_VZ_COL, MAXT_VR_COL), unpack=True)
      #  max velocity
      if vel_fname != None:
         (self._time, self._maxvmag, self._maxmachsub, self._maxmachfull, self._dt) = np.loadtxt(vel_fname, 
          usecols=(TIME_COL, MAXVMAG_COL, MAXMACH_SUB_COL, MAXMACH_FULL_COL, DT_COL), unpack=True)
      #  max enuc
      if enuc_fname != None:
         (self._time, self._maxenuc, self._maxenuc_x, self._maxenuc_y, self._maxenuc_z, self._maxenuc_r,
          self._maxenuc_vx, self._maxenuc_vy, self._maxenuc_vz, self._maxenuc_vr) = np.loadtxt(enuc_fname, 
          usecols=(TIME_COL, MAXENUC_COL, MAXENUC_X_COL, MAXENUC_Y_COL, MAXENUC_Z_COL, MAXENUC_R_COL, 
          MAXENUC_VX_COL, MAXENUC_VY_COL, MAXENUC_VZ_COL, MAXENUC_VR_COL), unpack=True)

   ##Public methods##
   def stageOutput(self):
      """Check the scratch directory for any updated output.  If found, copy to
      the stage directory."""
      #TODO: Implement this
      pass

   def getPeakState(self):
      """Get the temp and time at the time of peak global temperature."""
      import numpy as np

      Tpeak    = self._maxT.max()
      peakdex  = np.where(self._maxT == Tpeak)[0][0]
      #Tpeak = self._maxT[:peakdex-100].max()
      #peakdex = np.where(self._maxT == Tpeak)[0][0]
      timepeak = self._time[peakdex]
      rpeak    = self._maxT_r[peakdex]

      return Tpeak, timepeak, rpeak

   def getTempData(self):
      """Get tuple of data:
         (time, peak temp, peak temp's radius, radial velocity at peak temp)"""
      return self._time, self._maxT, self._maxT_r, self._maxT_vr
   
   def getMachData(self):
      """Get tuple of data:
         (time, peak Mach #)"""
      return self._time, self._maxmachfull

   def getEnucData(self):
      """Get tuple of data:
         (time, peak e_nuc)"""
      return self._time, self._maxenuc

   ##Private methods##


#TODO: Migrate self-plotting objects to plotting via SCPlotter
class SCPlotter(object):
   """A class for executing plots relevant to the Sub-Chandra problem."""
   ##Data##
   #Static variables 

   #Private variables 

   #Global constants

   #Public variables 

   #Constructor
   def __init__(self):
      """SCPlotter constructor.
         self       --> reference to this instance of SCPlotter (implicitly passed)"""
      pass

   ##Private methods##

   ##Public methods##
   def __str__(self):
      pass

   def printHotRhoTemp(self, sims):
      """Print 'step | time | <rho>_base | <temp>_base' for all sims"""
      import numpy as np
      import matplotlib.pyplot as plt

      print('step | time | <rho>_base | <temp>_base')
      for sim in sims:
         max_temp = 0.0
         max_i = -1
         max_rho = 0.0
         max_ts = None
         for s in sim.getSliceArray():
            avgTpeak = max(s.datalist[SCSlice.ITEMP])
            if avgTpeak > max_temp:
               max_temp = avgTpeak
               max_i = np.where(s.datalist[SCSlice.ITEMP] == avgTpeak)
               max_rho = s.datalist[SCSlice.IRHO][max_i][0]
               max_ts = s.timestamp
         print(max_ts[0], max_ts[1], max_rho, max_temp)   

   def plotHotRhoTemp(self, sims):
      """Print 'step | time | <rho>_base | <temp>_base' for all sims"""
      import numpy as np
      import matplotlib.pyplot as plt

      fig, ax = plt.subplots(nrows=1, ncols=1)
      y = []; x = []
      for sim in sims:
         max_temp = 0.0
         max_i = -1
         max_rho = 0.0
         max_ts = None
         for s in sim.getSliceArray():
            avgTpeak = max(s.datalist[SCSlice.ITEMP])
            if avgTpeak > max_temp:
               max_temp = avgTpeak
               max_i = np.where(s.datalist[SCSlice.ITEMP] == avgTpeak)
               max_rho = s.datalist[SCSlice.IRHO][max_i][0]
               max_ts = s.timestamp
         y.append(max_rho)
         x.append(max_temp)

      ax.scatter(x,y)

   def plotTempTimeSeries(self, sim):
      """Plot the core mass vs shell mass."""
      import numpy as np
      import matplotlib.pyplot as plt
      import matplotlib
      import scipy.integrate as spint

      #Convenience aliases for initial model data
      an_cut = sim._initmod._mycut.an_cut
      sp_cen_den = sim._initmod._mycut.sp_cen_den
      sp_st_fac = sim._initmod._mycut.sp_st_fac

      #Convenience aliases for diagnostic data
      t_T, Tpeak, Tpeak_r, Tpeak_vr  = sim.getDiagTemp()

      #Prepare variables for use in slice loop
      H = []
      iface = []
      rbot = []
      rtop = []
      t_slice = []

      tconv = []
      tconvb = []
      tnuc_x = []
      tnuc_xb = []
      tnuc_wk = []
      ratio = []

      avgTpeak = []
      avgBaseRho = []
      rhoCrit = []

      #Loop over slices in chronological order
      for s in sorted(sim.getSliceArray(), key=lambda sl: sl.timestamp[0]):
         #Build radius data from slices
         rbot.append(s.derived_datalist[SCSlice.IRCONV])
         rtop.append(s.derived_datalist[SCSlice.IRCONV] + s.derived_datalist[SCSlice.ILCONV])
         H.append(s.derived_datalist[SCSlice.IRCONV] + s.derived_datalist[SCSlice.IPSCALE])
         iface.append(s.derived_datalist[SCSlice.IINTFACE])
         t_slice.append(s.timestamp[1])

         #Estimate of convective turnover timescale and minimum nuclear timescale
         tconv.append(s.derived_datalist[SCSlice.ILCONV] / s.derived_datalist[SCSlice.IUCONV])
         tnuc_x.append(min(s.derived_datalist[SCSlice.ITNUC][0]))
         tnuc_wk.append(min(s.derived_datalist[SCSlice.ITNUC][1]))
         tnuc_xb.append(min(s.derived_datalist[SCSlice.ITNUC][2]))

         #Get the peak radially averaged temperature as an estimate of the background
         #conditions the hottest spot is being generated in.
         avgTpeak.append(max(s.datalist[SCSlice.ITEMP]))
         brho_i = np.where(s.datalist[SCSlice.ITEMP] == avgTpeak[len(avgTpeak)-1])
         avgBaseRho.append(s.datalist[SCSlice.IRHO][brho_i])
         t8 = avgTpeak[-1:][0] / 1.e8
         rctemp = (1.68e-4*np.exp(20.0/t8))**(1.0/2.3)
         rhoCrit.append(rctemp*1.e6)

         ###TEST: Calculate global convective timescale
         #Calculate cutoff radii
         sp_st_den = sp_cen_den*sp_st_fac
         r_anelastic = None
         r_sp_start = None
         for r, rho in zip(s.datalist[SCSlice.IRADIUS], s.datalist[SCSlice.IRHO]):
           #TODO add error checking
           if not r_anelastic and rho <= an_cut:
             r_anelastic = r
             if r_sp_start:
               break
           if not r_sp_start and rho <= sp_st_den:
             r_sp_start = r
             if r_anelastic:
               break
         cbot = s.derived_datalist[SCSlice.IRCONV]
         ctop = cbot + s.derived_datalist[SCSlice.ILCONV]
         magvel = s.datalist[SCSlice.IMV]
         mv_rad = s.datalist[SCSlice.IRADIUS]

         #Change to dimensionless variables, only care about the convective zone
         li = np.where(mv_rad == cbot)[0]
         ri = np.where(mv_rad == ctop)[0]
         r_norm = mv_rad[li:ri]/ctop
         magvel_norm = magvel[li:ri] / magvel.max()
         mvn_inv = 1.0 / magvel_norm

         #Calculate global convective timescale as integral of 1/v over the convective zone
         #Convert back to physical units
         tcg = (ctop/magvel.max())*spint.trapz(mvn_inv, r_norm)
         tconvb.append(tcg)

         ###END TEST

         ratio.append(tnuc_x[len(tnuc_x)-1]/tconvb[len(tconvb)-1])

      #Build plots
      fig, ax_list = plt.subplots(nrows=1, ncols=1)
      fig.set_size_inches(12.0, 10.0)

      #Temp
      ax_list.plot(t_T, Tpeak, color='red')
      ax_list.plot(t_slice, avgTpeak, color='red', marker='x', linestyle='None', label='avg peak temp')
      ax_list.set_ylabel(r'T$_{\mathrm{peak}}$ (K)', color='red')
      ax_list.set_title(sim.getLabel() + ' | Peak Temperature')
      #ax_list[0].set_ylim(1.74e8, 4.0e8)
      #ax_list[0].set_xlim(150, 155)
      tw = ax_list.twinx()
      #tw.set_xlim(150, 155)
      tw.plot(t_T, Tpeak_r, color='green')
      tw.plot(t_slice, iface, color='cyan', marker='o', linestyle='None', label='CO/He')
      tw.plot(t_slice, H, color='cyan', marker='^', linestyle='None', label='H')
      tw.plot(t_slice, rbot, color='cyan', marker='v', linestyle='None', label='rbot')
      tw.plot(t_slice, rtop, color='cyan', marker='v', linestyle='None', label='rtop')
      tw.set_ylabel(r'T$_{\mathrm{peak}}$ radius (cm)', color='green')
      #handles, labels = ax_list[0].get_legend_handles_labels()
      #fig.legend(handles, labels)
      tw.legend(loc=2)

   def plotMinMass(self, core, shell, fontsize='medium'):
      """Plot the core mass vs shell mass."""
      import numpy as np
      import matplotlib.pyplot as plt
      import matplotlib
      from matplotlib.patches import Polygon

      # Make Bildsten et al. 2007 data
      bswn_core = np.array([0.8, 1.0, 1.2])
      bswn_shell =     np.array([0.095, 0.0425, 0.012])
      bswn_shell_err = np.array([[0.07, 0.0275, 0.0075],
                                 [0.12, 0.055,  0.018]])
      vertices = [(0.8, 0.07),  (1.0, 0.0275), (1.2, 0.0075), 
                  (1.2, 0.018), (1.0, 0.055),  (0.8, 0.12)]
      bswn_patch = Polygon(vertices, facecolor='0.9', edgecolor='1.0', label='BSWN2007')

      # Make Woosley & Kasen 2011 data
      wk_core = np.array([0.8, 1.0, 1.1])
      wk_shell =     np.array([0.085,  0.05,   0.0375])
      wk_shell_err = np.array([[0.075, 0.0375, 0.025],
                               [0.095, 0.062,  0.05]])
      vertices = [(0.8, 0.075),  (1.0, 0.0375), (1.1, 0.025), 
                  (1.1, 0.05), (1.0, 0.062),  (0.8, 0.095)]
      wk_patch = Polygon(vertices, hatch='x', facecolor='1.0', edgecolor='0.0', label=r'W\&K2011')

      #make the plots
      fig, ax = plt.subplots(nrows=1, ncols=1)

      #plot the data
      ax.plot(core, shell, 'r', label='this work')
      ax.add_patch(bswn_patch)
      ax.add_patch(wk_patch)
      #ax.errorbar(bswn_core, bswn_shell, yerr=bswn_shell_err, fmt='--k', label='BSWN2007')
      #ax.errorbar(wk_core, wk_shell, yerr=wk_shell_err, fmt='-k', label='W&K2011')

      #tune plot settings
      ax.set_yscale('log')
      ax.set_ylabel(r'Helium shell mass (M$_\odot$)', fontsize=fontsize)
      ax.set_ylim(0.001, 0.2)

      ax.set_xlabel(r'Core mass (M$_\odot$)', fontsize=fontsize)
      ax.set_xlim(0.75, 1.35)

      ax.tick_params(labelsize=fontsize)
      ax.legend(loc=1)

      #save fig for pub
      fig.savefig("minmass.eps", bbox_inches='tight')

   def plotHotspots(self, scsim, step, rlim=None, templog=False, reset=False, paper=True, plot_top=False):
      """Plot radius and temperature histograms for this timestep's top hotspots 
      as well as temperature contours."""
      import matplotlib
      #matplotlib.use('Agg')
      import matplotlib.pyplot as plt
      import numpy as np
      from matplotlib import cm, colors
      from mpl_toolkits_ext.basemap import Basemap#, cm

      TCRIT = 2.25e7
      TMAX = 2.4e9
      PFONT_SIZE = 'large'
      SFONT_SIZE = 'x-large'
      RSCALE = 1.0e8
      TSCALE = 1.0e8
      RAD2DEG = 180./np.pi
      CM2M =  1./100.
      CM2KM = 1.0/1.e5
      OCTANT_OMEGA = 90.*90. #Square degrees of an octant surface
         
      matplotlib.rc('text', usetex=False)

      print(matplotlib.__version__)
      #First I need the interesting radii details.  
      #TODO: This is a stupid inefficient way to get them.  Need to rewrite/restructure.
      radii = (None, None, None) 
      for s in scsim._out._slices:
         if s.timestamp[0] == step:
            radii = (s.derived_datalist[SCSlice.IRCONV],
                  s.derived_datalist[SCSlice.IPSCALE],
                  s.derived_datalist[SCSlice.IINTFACE])

      dx = scsim.getFineResolution()
      print('rad, fine_dx: ')
      print(radii)
      print(dx)
      #Find the right set of hotspots
      hspots = None
      for hs in scsim._out._hotspots:
         if hs.timestamp[0] == step:
            hspots = hs
            break

      if (hspots is None):
         print("Couldn't find step {0:6d}".format(step))
         return

      #When tweaking plots it's nice if I can store the data, but then when I load new
      #data I need to reset.  Here I delete all data arrays so they'll be rebuilt.
      if(reset and hasattr(hspots, 'r')):
         del hspots.r
         del hspots.theta
         del hspots.phi
         del hspots.temp
         del hspots.logtemp
         del hspots.rho6
         del hspots.temp_lons
         del hspots.temp_lats

      #Get data, and save data so we only go through all the arrays once
      if not hspots._hsarray: #If no array yet, build it
         hspots._buildHSArr()
      if not hasattr(hspots, 'r'):
         hspots.r         = np.array([hs.loc[1][0] for hs in hspots._hsarray if hs.temp > TCRIT and hs.temp < TMAX])
      if not hasattr(hspots, 'theta'):
         hspots.theta     = np.array([hs.loc[1][1] for hs in hspots._hsarray if hs.temp > TCRIT and hs.temp < TMAX])
      if not hasattr(hspots, 'phi'):
         hspots.phi       = np.array([hs.loc[1][2] for hs in hspots._hsarray if hs.temp > TCRIT and hs.temp < TMAX])
      if not hasattr(hspots, 'temp'):
         hspots.temp      = np.array([hs.temp/TSCALE for hs in hspots._hsarray if hs.temp > TCRIT and hs.temp < TMAX])
         #hspots.temp      = np.array([hs.temp for hs in hspots._hsarray if hs.temp > TCRIT and hs.temp < TMAX])
      if not hasattr(hspots, 'logtemp'):
         hspots.logtemp   = np.log10(hspots.temp)
      if not hasattr(hspots, 'rho6'):
         hspots.rho6      = np.array([hs.rho/1.e6  for hs in hspots._hsarray if hs.temp > TCRIT and hs.temp < TMAX])
      if not hasattr(hspots, 'temp_lons'):
         hspots.temp_lons = np.array([deg for deg in np.degrees(hspots.phi)])
      if not hasattr(hspots, 'temp_lats'):
         hspots.temp_lats = np.array([-(deg-90.) for deg in np.degrees(hspots.theta)])

      #Local aliases for data
      r         = hspots.r
      avg_r     = np.average(r)  #Average radius in cm
      #theta    = hspots.theta
      #phi      = hspots.phi
      temp      = hspots.temp
      logtemp   = hspots.logtemp
      rho6      = hspots.rho6
      temp_lons = hspots.temp_lons
      temp_lats = hspots.temp_lats

      #Critical hotspots
      #ctemp      = np.array([hs.temp for hs in self._hsarray if hs.rho/1.e6 > (1.68e-4*np.exp(20.0/(hs.temp/1.e8)))**(1.0/2.3)])
      #ctheta     = np.array([hs.loc[1][1] for hs in self._hsarray if hs.rho/1.e6 > (1.68e-4*np.exp(20.0/(hs.temp/1.e8)))**(1.0/2.3)])
      #cphi       = np.array([hs.loc[1][2] for hs in self._hsarray if hs.rho/1.e6 > (1.68e-4*np.exp(20.0/(hs.temp/1.e8)))**(1.0/2.3)])
      #ctemp_lons = np.array([deg for deg in np.degrees(cphi)])
      #ctemp_lats = np.array([-(deg-90.) for deg in np.degrees(ctheta)])

      crit_temp = 2.3e8
      #ctemp      = np.array([hs.temp for hs in self._hsarray if hs.temp > crit_temp])
      #ctheta     = np.array([hs.loc[1][1] for hs in self._hsarray if hs.temp > crit_temp])
      #cphi       = np.array([hs.loc[1][2] for hs in self._hsarray if hs.temp > crit_temp])
      #ctemp_lons = np.array([deg for deg in np.degrees(cphi)])
      #ctemp_lats = np.array([-(deg-90.) for deg in np.degrees(ctheta)])

      #Get important radii of interest
      rbot = radii[0]
      H = radii[1]
      iface = radii[2]

      #Get min, max temp
      min_temp = temp.min()
      max_temp = temp.max()
      hs_count = len(temp)

      #Calculate temperature bins for color map
      #shsarr = sorted(hspots._hsarray, key=lambda hs: hs.temp)
      stemp = sorted(temp)
      tlevs = [temp.min()]
      tllevs = [logtemp.min()]
      count = 0
      for t in stemp:
         count += 1
         if count > (1./9.)*hs_count:
            tlevs.append(t)
            tllevs.append(np.log10(t))
            count = 0
         #For hottest temps break into top 9% and top 1%
         #if len(tlevs) == 9:
         #  if count > 0.09*hs_count:
         #    tlevs.append(hs.temp)
         #    tllevs.append(np.log10(hs.temp))
         #    count = 0
         #else:
         #  if count > 0.1*hs_count:
         #    tlevs.append(hs.temp)
         #    tllevs.append(np.log10(hs.temp))
         #    count = 0
      else:
         tlevs.append(t)
         tllevs.append(np.log10(t))

      #Build colormap
      #TODO: Get nice paper-worthy plot for 'ignition' section
      #TODO: Figure out discrepancy between pltfile cells and valid cells
      #      Update: currently working with OLCF on this, seems to only happen with optimized Cray compiler
      #rr = np.linspace(1.0, 0.7, 11)
      #gb = np.linspace(1.0, 0.0, 11)
      #temp_cmap = colors.ListedColormap([
      #  (rr[0],  gb[0],  gb[0]),   #Coldest
      #  (rr[1],  gb[1],  gb[1]),
      #  (rr[2],  gb[2],  gb[2]),
      #  (rr[3],  gb[3],  gb[3]),
      #  (rr[4],  gb[4],  gb[4]),
      #  (rr[5],  gb[5],  gb[5]),
      #  (rr[6],  gb[6],  gb[6]),
      #  (rr[7],  gb[7],  gb[7]),
      #  (rr[8],  gb[8],  gb[8]),
      #  (rr[9],  gb[9],  gb[9]),
      #  (rr[10], gb[10], gb[10])]) #Hottest
      #  #(1,      1,      0)]) #Hottest
      # RGB colors from 9-class OrRd at colorbrewer2.org
      temp_cmap = colors.ListedColormap([
         (255./255., 247./255., 236./255.),   #Coldest
         (254./255., 232./255., 200./255.),
         (253./255., 212./255., 158./255.),
         (253./255., 187./255., 132./255.),
         (252./255., 141./255.,  89./255.),
         (239./255., 101./255.,  72./255.),
         (215./255.,  48./255.,  31./255.),
         (179./255.,   0./255.,   0./255.),
         (127./255.,   0./255.,   0./255.)]) #Hottest
      tc_bounds = tlevs
      tc_norm = colors.BoundaryNorm(tc_bounds, temp_cmap.N)
      tmap = cm.ScalarMappable(norm=tc_norm, cmap=temp_cmap)
      tmap._A = [] # mpl >= 1.2 want ScalarMAppables to have an _A array.

      #Calculate critical density range based on eqn (8) of Woosley & Kasen 2011
      t8 = max_temp #/1.e8
      rho_critl = (1.68e-4*np.exp(20.0/t8))**(1.0/2.3)
      t8 = min_temp #/1.e8
      rho_critr = (1.68e-4*np.exp(20.0/t8))**(1.0/2.3)

      print('Min, Max temp of {0} hottest cells: {1}, {2}'.format(hs_count, min_temp, max_temp))
      print('Critical density range: [{0}, {1}]'.format(rho_critl, rho_critr))
      print('Scale height: {0}'.format(H))
      sys.stdout.flush()

      if paper:
         #Build plots
         print('mark A')
         sys.stdout.flush()
         plt.clf()
         fig = plt.figure()
         #subplot2grid call signature: (grid_row, grid_cols), (subplot_row, subplot_col), colspan=1, rowspan=1 
         ax_rad  = plt.subplot2grid((3,3), (0,0))  
         ax_temp = plt.subplot2grid((3,3), (0,1))  
         ax_rho  = plt.subplot2grid((3,3), (0,2))  
         ax_proj = plt.subplot2grid((3,3), (1,0), colspan=2, rowspan=2)  

         #Plot temperature histogram
         ax_temp.hist(temp, bins=1000)
         ax_temp.set_xlabel("temperature (K)")

         #Build projection map for temperature theta, phi locations
         map = Basemap(projection='nsper', lon_0=45, lat_0=45, 
               #llcrnrlon=-180, llcrnrlat=-90, urcrnrlon=180, urcrnrlat=90, 
               resolution=None, ax=ax_proj, rsphere=avg_r*CM2M)
         #map.drawmeridians(np.arange(0, 90, 15), color="0.65", latmax=90)
         map.drawparallels([0, 80], color="0.65", latmax=90) #, labels=[1,0,0,1])
         map.drawmapscale(-15,45,45,45,H*CM2KM,labelstyle=False,
               format='%6.2f', fontsize=11) #Draw scale height

         #It's stupid that I have to do this, but below I "erase" extraneous latitude lines
         #by writing over them with thick white lines.
         #for lat in range(0,90,15):
         #  for long in range(15,180,15):
         #    #Erase extraneous latitude lines on the left
         #    left_long=-long
         #    map.drawgreatcircle(left_long+15, lat, left_long, lat, linewidth=5, color="w")
         #    #Erase extraneous latitude lines on the right
         #    right_long=long+90
         #    map.drawgreatcircle(right_long-15, lat, right_long, lat, linewidth=5, color="w")
         ##Same with extraneous longitude lines at the bottom
         #map.drawgreatcircle(0,    0,  0, -25, linewidth=5, color="w")
         #map.drawgreatcircle(15,   0, 15, -30, linewidth=5, color="w")
         #map.drawgreatcircle(30,   0, 30, -35, linewidth=5, color="w")
         #map.drawgreatcircle(45,   0, 45, -30, linewidth=5, color="w")
         #map.drawgreatcircle(60,   0, 60, -35, linewidth=5, color="w")
         #map.drawgreatcircle(75,   0, 75, -30, linewidth=5, color="w")

         print('mark B')
         sys.stdout.flush()
         # draw the boundary of our domain -- we want great circles here
         # note that we draw in 15 degree increments.  Otherwise the lat/long grid
         # doesn't line up with the boundary 
         #Left boundary
         for lat in range(0,90,15):
            map.drawgreatcircle(0, lat, 0, lat+15, linewidth=1, color="k", zorder=max_temp+10)
         #Right boundary
         for lat in range(0,90,15):
            map.drawgreatcircle(90, lat, 90, lat+15, linewidth=1, color="k", zorder=max_temp+10)
         #Bottom boundary
         for lon in range(0,90,15):
            map.drawgreatcircle(lon, 0, lon+15, 0, linewidth=1, color="k", zorder=max_temp+10)

         print('mark C')
         sys.stdout.flush()
         if templog:
            clevs = np.linspace(logtemp.min(), logtemp.max(), 11) 
            #cs = map.contourf(temp_lons, temp_lats, logtemp, clevs, latlon=True, tri=True, cmap=cm.Reds)
            cs = map.contourf(temp_lons, temp_lats, logtemp, tllevs, latlon=True, tri=True, cmap=cm.jet)
         else:
            #clevs = np.linspace(temp.min(), temp.max(), 11) 
            clevs = np.linspace(2.25e8, crit_temp, 11) 
            #cs = map.contourf(temp_lons, temp_lats, temp, tlevs, latlon=True, tri=True, cmap=temp_cmap, norm=tc_norm)
            #map.contourf(ctemp_lons, ctemp_lats, ctemp, clevs, latlon=True, tri=True, cmap=cm.Greens)
         #cbar = map.colorbar(cs, location='right', pad='5%', ticks=tlevs)
         #cbar.set_label('Kelvin')
         temp_cols =  tmap.to_rgba(temp)
         deg_fine = dx/avg_r * RAD2DEG 
         #deg_fine = deg_fine*1.e2
         print('mark Cx')
         sys.stdout.flush()
         tx, ty = map(temp_lons, temp_lats)
         num_bins = int(np.round(90./deg_fine))
         #map.hexbin(tx,ty,C=temp,gridsize=num_bins,cmap=temp_cmap,norm=tc_norm,
         #     reduce_C_function=np.max)
         #    reduce_C_function=np.max,rasterized=True)
         i = 0
         for tln, tlt, t, c in zip(temp_lons, temp_lats, temp, temp_cols):
            i = i + 1
            if i % 10 == 0:
               map.tissot(tln, tlt, deg_fine, 4, zorder=t, color=c)
            
            #map.tissot(tln, tlt, deg_fine, 4, zorder=t, color=c, rasterized=True)
         #map.tissot(45,45,10,4)

         cbar = map.colorbar(tmap, location='right', pad='5%', ticks=tlevs, ax=ax_proj,
                             format='%.3f')
         cbar.set_label(r'temperature ($\times 10^8$ Kelvin)', fontsize='x-large')

         #Plot radius histogram with temperature color-coding
         # color-coding is achieved by plotting several bars instead of using
         # ax.hist(), and we use the projection map's colorbar
         # WD/He iface = radius where <X_He> = 0.9
         #ax_rad.hist(r, bins=1000)

         dr = dx #use a dr of dx, which is roughly the radial resolution in these simulations
         #dr = 1.e5 #use a dr of 1 km, which is roughly the radial resolution in these simulations
         radii, counts, cols = hspots._binData(r, temp, dr, cbar)

         print('mark D')
         sys.stdout.flush()
         for r, c, col in zip(radii, counts, cols):
            ax_rad.bar(r/RSCALE, c, width=dr/RSCALE, color=col, edgecolor=col, align='center') 
         #ax_rad.bar(radii, counts, width=dr, color=(1.0, 1.0, 1.0), align='center') 
         ax_rad.set_xlabel(r"radius ($\times 10^8$ cm)", fontsize='x-large')
         if rlim:
            ax_rad.set_xlim(rlim[0]/RSCALE, rlim[1]/RSCALE)
         # plot the radii (CO/He interface, start of convective region, top of convective region)
         ax_rad.plot([iface/RSCALE, iface/RSCALE], [0, ax_rad.get_ylim()[1]], color="k", linestyle='--')
         ax_rad.plot([rbot/RSCALE, rbot/RSCALE], [0, ax_rad.get_ylim()[1]], color="k", linestyle='--')
         if plot_top:
            ax_rad.plot([(rbot + H)/RSCALE, (rbot + H)/RSCALE], [0, ax_rad.get_ylim()[1]], color="k", linestyle='--')

         #Annotate the interface line
         ifrac = (iface/RSCALE - ax_rad.get_xlim()[0]) / (ax_rad.get_xlim()[1] - ax_rad.get_xlim()[0])
         ax_rad.annotate(r'$\langle X_\mathrm{He}\rangle_r = 0.9$', 
               xy=(ifrac,0.9), 
               xytext=(-80,-30),
               xycoords='axes fraction', textcoords='offset points',
               fontsize = PFONT_SIZE, fontweight='bold',
               arrowprops=dict(facecolor='black', arrowstyle='simple',
                  connectionstyle='arc3,rad=-0.2'))

         #Annotate the convective base line
         cbfrac = (rbot/RSCALE - ax_rad.get_xlim()[0]) / (ax_rad.get_xlim()[1] - ax_rad.get_xlim()[0])
         ax_rad.annotate('Convective\nbase', 
               xy=(cbfrac,0.8), 
               xytext=(30,-30),
               xycoords='axes fraction', textcoords='offset points',
               arrowprops=dict(facecolor='black', arrowstyle='simple',
                  connectionstyle='arc3,rad=-0.2'))

         #Annotate the pressure scale height
         ax_proj.annotate(r'Scale Height [km]'.format(H*CM2KM),
               xy=(0.06,0.675), 
               xytext=(0,0),
               xycoords='axes fraction', textcoords='offset points')

         print('mark E')
         sys.stdout.flush()
         #Plot density histogram with color-coding
         drho6 = 0.001
         rho6_bins, rho_counts, rho_colors = hspots._binData(rho6, temp, drho6, cbar)
         for r6, c, col in zip(rho6_bins, rho_counts, rho_colors):
            ax_rho.bar(r6, c, width=drho6, color=col, edgecolor=col, align='center') 
         #ax_rho.hist(rho6, bins=1000)
         ax_rho.set_xlabel(r"density ($\times 10^6$ g cm$^{-3}$)", fontsize='x-large')
         #ax_rho.fill_between([rho_critl, rho_critr], 0, 500, facecolor='0.9', edgecolor='1.0')
         ax_rho.plot([rho_critr, rho_critr], [0, ax_rho.get_ylim()[1]], color="k", linestyle='--')

         #Annotate the critical density line
         rhofrac = (rho_critr - ax_rho.get_xlim()[0]) / (ax_rho.get_xlim()[1] - ax_rho.get_xlim()[0])
         ax_rho.annotate(r'$\rho_{\mathrm{cr},\mathrm{WK}}$', 
               xy=(rhofrac,0.8), size=12.5,
               xytext=(30,-30),
               xycoords='axes fraction', textcoords='offset points',
               arrowprops=dict(facecolor='black', arrowstyle='simple',
                  connectionstyle='arc3,rad=-0.2'))

         print('mark F')
         sys.stdout.flush()
         #Set plot properties 
         fig.set_size_inches(15.0, 12.5)
         #This fixes a problem with mpl's pstoeps converter when using ghostscript as distiller
         #matplotlib.rc('ps', usedistiller='xpdf')
         #fig.savefig("test_ig.png", bbox_inches='tight')
         #fig.savefig("test_ig.pdf")
         print('mark G')
         sys.stdout.flush()

      else:
         #TODO: Need to implement non-inline plotting
         pass

   def plotTHist(self, scsim, step):
      """Plot temperature histograms from the given timestep and simulation."""
      import matplotlib
      #matplotlib.use('Agg')
      import matplotlib.pyplot as plt
      import numpy as np

      #Get the data
      thd = scsim.getTHistData(step)
      ells, temps, counts = thd.ells, thd.temps, thd.counts
      TMax = thd.TMax

      #make the plots
      fig, ax = plt.subplots(nrows=len(ells), ncols=2)

      #Plot first lengthscale
      for i in range(len(ells)):
         #Data for cell's stored temperature
         ax[i][0].bar(temps[i], counts[0][i], width=1.0e5)
         ax[i][0].set_title('l={0}'.format(ells[i]))
         text_bbox = dict(boxstyle="round", fc="white")
         ax[i][0].text(0.3, 0.8, r'$T_\mathrm{{max}} = {0}$'.format(TMax[i]), ha='center', va='center', transform=ax[i][0].transAxes,
                    bbox=text_bbox, size=15)
         
         #Data for cell's temperature derived from rho, h, X
         ax[i][1].bar(temps[i], counts[1][i], width=1.0e5)
         ax[i][1].set_title('l={0}'.format(ells[i]))
         text_bbox = dict(boxstyle="round", fc="white")
         ax[i][1].text(0.3, 0.8, r'$T_\mathrm{{max}} = {0}$'.format(TMax[i]), ha='center', va='center', transform=ax[i][1].transAxes,
                    bbox=text_bbox, size=15)

      #ax[1].bar(temps[1], counts[1], width=1.0e5)
      #ax[1].set_title('l={0}'.format(ells[1]))
      #ax[2].bar(temps[2], counts[2], width=1.0e5)
      #ax[2].set_title('l={0}'.format(ells[2]))
      #ax[3].bar(temps[3], counts[3], width=1.0e5)
      #ax[3].set_title('l={0}'.format(ells[3]))

      fig.set_size_inches(10.0, 15.0)

   def plotInitModel(self, sim, writefile=False, fontsize='x-large'):
      """Plots this initial model using matplotlib."""
      import matplotlib.pyplot as plt
      import numpy as np
      from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
      from mpl_toolkits.axes_grid1.inset_locator import mark_inset
 
      #NOTE: fontsize can be [size in points | 'xx-small' | 'x-small' |
      #  'small' | 'medium' | 'large' | 'x-large' | 'xx-large'
      
      #This function is just an alias for the following call 
      self.plotSliceOverview(sim, None, writefile, fontsize)
  
   def plotSliceOverview(self, sim, time, writefile=False, fontsize='x-large'):
      """Plot overview of the available slice data nearest to time.  Set time to None
      to plot the initial model data."""
      import matplotlib.pyplot as plt
      import numpy as np
      from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
      from mpl_toolkits.axes_grid1.inset_locator import mark_inset
  
      #NOTE: fontsize can be [size in points | 'xx-small' | 'x-small' |
      #  'small' | 'medium' | 'large' | 'x-large' | 'xx-large'
  
      # make the plots
      fig, ax_list = plt.subplots(nrows=3, ncols=1)

      #Find the right slice, if requested
      slicedata = None
      if time is not None:
         for s in sorted(sim._out._slices, key=lambda sl: sl.timestamp[0]):
            if s.timestamp[1] > time:
               slicedata = s
               break

      #Plot the three subplots
      self._plotRhoT(ax_list[0], sim, slicedata, fontsize)
      self._plotX(   ax_list[1], sim, slicedata, fontsize)
      self._plotS(   ax_list[2], sim, slicedata, fontsize)

      #Adjust figure settings
      fig.set_size_inches(6.0,12.0)
      fig.tight_layout()
  
      #If requested, save a figure 
      if writefile:
         if time is None:
            fig.savefig(sim._label + "_initial_model.eps", bbox_inches='tight')
         else:
            fig.savefig(sim._label + "_slice_overview.eps", bbox_inches='tight')

      #Write header for interactive viewing (not for the file)
      if time is None:
         fig.suptitle(sim._label + ' Initial Model', fontsize=18, y=1.00)
      else:
         fig.suptitle(sim._label + ' Slice', fontsize=18, y=1.00)

   def displayParameters(self, sim):
      """Use IPython's display framework to write out a table of key parameters"""
      from glob import glob
      from IPython.display import HTML, display
      from os.path import dirname

      #Assume <run label>/[output | run | plots] naming scheme

      #Write opening
      htmlstr = r'<table border="1">' + "\n"

      #I'm repurposing functions that originally needed a list, so put the 
      #label in a list
      #TODO: Need to rework this so I'm not using single item lists
      labellist = [sim.getLabel(),]
      stage_dir = dirname(sim.getStageDir())

      #Add each row
      htmlstr += self._labelRow(labellist)
      htmlstr += self._massRow(labellist, stage_dir)
      htmlstr += self._paramRow(labellist, stage_dir, r'T$_{\mathrm{CO}}$ [K]', 'temp_core')
      htmlstr += self._paramRow(labellist, stage_dir, r'T$_{\mathrm{base}}$ [K]', 'temp_base')
      htmlstr += self._rhobaseRow(labellist, stage_dir, r'$\rho_{\mathrm{base}}$ [$\times 10^5$ g cm$^{-3}$]')
      htmlstr += self._paramRow(labellist, stage_dir, r'x$_{\mathrm{max}}$ [cm]', 'xmax')
      htmlstr += self._inputsRow(labellist, stage_dir, r'anelastic_cutoff [g cm$^{-3}$]', 'anelastic_cutoff')
      htmlstr += self._inputsRow(labellist, stage_dir, r'base_cutoff_density [g cm$^{-3}$]', 'base_cutoff_density')
      htmlstr += self._inputsRow(labellist, stage_dir, r'species_pred_type', 'species_pred_type')
      htmlstr += self._inputsRow(labellist, stage_dir, r'octant', 'octant')
      htmlstr += self._inputsRow(labellist, stage_dir, r'Levs', 'max_levs')
      htmlstr += self._dxRow(labellist, stage_dir, r'$\Delta x_\mathrm{fine}$ [km]')

      #Write closing
      htmlstr += r'</table>'
      
      display(HTML(htmlstr))

   def plotTimeSeries(self, sim):  #an_cut, sp_cen_den, sp_st_fac):
      """Plot data of interest vs time for this simulation's output."""
      import numpy as np
      import matplotlib.pyplot as plt
      import matplotlib
      import scipy.integrate as spint

      matplotlib.rc('text', usetex=False)

      #Convenience aliases for diagnostic data
      t_T, Tpeak, Tpeak_r, Tpeak_vr  = sim.getDiagTemp()
      t_M, Mpeak = sim.getDiagMach()
      t_e, epeak = sim.getDiagEnuc()

      #Prepare variables for use in slice loop
      H = []
      iface = []
      rbot = []
      rtop = []
      t_slice = []

      tconv = []
      tconvb = []
      tnuc_x = []
      tnuc_xb = []
      tnuc_wk = []
      ratio = []

      avgTpeak = []
      avgBaseRho = []
      rhoCrit = []

      #Get some data
      slices = sim.getSliceArray()
      cuts = sim.getIMCuts()
      an_cut = cuts.an_cut
      sp_cen_den = cuts.sp_cen_den
      sp_st_fac = cuts.sp_st_fac

      #Loop over slices in chronological order
      for s in sorted(slices, key=lambda sl: sl.timestamp[0]):
         #Build radius data from slices
         rbot.append(s.derived_datalist[SCSlice.IRCONV])
         rtop.append(s.derived_datalist[SCSlice.IRCONV] + s.derived_datalist[SCSlice.ILCONV])
         H.append(s.derived_datalist[SCSlice.IRCONV] + s.derived_datalist[SCSlice.IPSCALE])
         iface.append(s.derived_datalist[SCSlice.IINTFACE])
         t_slice.append(s.timestamp[1])

         #Estimate of convective turnover timescale and minimum nuclear timescale
         #tconv.append(s.derived_datalist[SCSlice.IPSCALE] / s.derived_datalist[SCSlice.IUCONV])
         tconv.append(s.derived_datalist[SCSlice.ILCONV] / s.derived_datalist[SCSlice.IUCONV])
         #tconv.append(s.derived_datalist[SCSlice.IUCONV])
         tnuc_x.append(min(s.derived_datalist[SCSlice.ITNUC][0]))
         tnuc_wk.append(min(s.derived_datalist[SCSlice.ITNUC][1]))
         tnuc_xb.append(min(s.derived_datalist[SCSlice.ITNUC][2]))
         #ratio.append(tnuc_x[len(tconv)-1]/tconv[len(tconv)-1])
         #ratio.append(tnuc_wk[len(tconv)-1]/tconv[len(tconv)-1])

         #Get the peak radially averaged temperature as an estimate of the background
         #conditions the hottest spot is being generated in.
         avgTpeak.append(max(s.datalist[SCSlice.ITEMP]))
         brho_i = np.where(s.datalist[SCSlice.ITEMP] == avgTpeak[len(avgTpeak)-1])
         avgBaseRho.append(s.datalist[SCSlice.IRHO][brho_i])
         t8 = avgTpeak[-1:][0] / 1.e8
         rctemp = (1.68e-4*np.exp(20.0/t8))**(1.0/2.3)
         rhoCrit.append(rctemp*1.e6)

         ###TEST: Calculate global convective timescale
         #Calculate cutoff radii
         sp_st_den = sp_cen_den*sp_st_fac
         r_anelastic = None
         r_sp_start = None
         for r, rho in zip(s.datalist[SCSlice.IRADIUS], s.datalist[SCSlice.IRHO]):
            #TODO add error checking
            if not r_anelastic and rho <= an_cut:
               r_anelastic = r
               if r_sp_start:
                  break
            if not r_sp_start and rho <= sp_st_den:
               r_sp_start = r
               if r_anelastic:
                  break
         cbot = s.derived_datalist[SCSlice.IRCONV]
         ctop = cbot + s.derived_datalist[SCSlice.ILCONV]
         magvel = s.datalist[SCSlice.IMV]
         mv_rad = s.datalist[SCSlice.IRADIUS]

         #Change to dimensionless variables, only care about the convective zone
         li = np.where(mv_rad == cbot)[0]
         ri = np.where(mv_rad == ctop)[0]
         r_norm = mv_rad[li:ri]/ctop
         magvel_norm = magvel[li:ri] / magvel.max()
         mvn_inv = 1.0 / magvel_norm

         #Calculate global convective timescale as integral of 1/v over the convective zone
         #Convert back to physical units
         tcg = (ctop/magvel.max())*spint.trapz(mvn_inv, r_norm)
         tconvb.append(tcg)

         ###END TEST

         ratio.append(tnuc_x[len(tnuc_x)-1]/tconvb[len(tconvb)-1])

      if 'inline' in matplotlib.get_backend():
         #Build plots
         fig, ax_list = plt.subplots(nrows=5, ncols=1)

         #Temp
         ax_list[0].plot(t_T, Tpeak, color='red')
         ax_list[0].plot(t_slice, avgTpeak, color='red', marker='x', linestyle='None', label='avg peak temp')
         ax_list[0].set_ylabel(r'T$_{\mathrm{peak}}$ (K)', color='red')
         ax_list[0].set_title(sim.getLabel() + ' | Peak Temperature')
         #ax_list[0].set_ylim(1.74e8, 4.0e8)
         #ax_list[0].set_xlim(150, 155)
         tw = ax_list[0].twinx()
         #tw.set_xlim(150, 155)
         #tw.set_ylim(4.26e8, 4.38e8)
         tw.plot(t_T, Tpeak_r, color='green')
         tw.plot(t_slice, iface, color='cyan', marker='o', linestyle='None', label='CO/He')
         tw.plot(t_slice, H, color='cyan', marker='^', linestyle='None', label='H')
         tw.plot(t_slice, rbot, color='cyan', marker='v', linestyle='None', label='rbot')
         tw.plot(t_slice, rtop, color='cyan', marker='v', linestyle='None', label='rtop')
         tw.set_ylabel(r'T$_{\mathrm{peak}}$ radius (cm)', color='green')
         #handles, labels = ax_list[0].get_legend_handles_labels()
         #fig.legend(handles, labels)
         tw.legend(loc=2)

         #Mach
         ax_list[1].plot(t_M, Mpeak, color='blue')
         ax_list[1].set_title(sim.getLabel() + ' | Peak Mach #')
         #ax_list[1].set_xlim(0, 70)
         #ax_list[1].set_ylim(0, .045)

         #enuc
         ax_list[2].plot(t_e, epeak, color='blue')
         ax_list[2].set_title(sim.getLabel() + ' | Peak H_nuc')
         #ax_list[2].set_xlim(0, 70)
         #ax_list[2].set_ylim(0.5e13, 4.e13)

         #Timescales
         ax_list[3].set_yscale('log')

         #11050-107-175-3lev
         #ax_list[3].set_ylim(1.e3, 1.e4)
         #ax_list[3].set_xlim(145.0, 155.0)

         #12025-107-175-3lev
         #ax_list[3].set_ylim(1.e3, 1.e4)
         #ax_list[3].set_xlim(200.0, 225.0)

         #12050-107-165-3lev
         #ax_list[3].set_ylim(1.e3, 1.e4)
         #ax_list[3].set_xlim(60.0, 70.0)

         ax_list[3].plot(t_slice, tconv, color='blue', label=r'$\tau_{conv}$')
         ax_list[3].plot(t_slice, tconvb, color='blue', linestyle='--', label=r'$\tau_{conv,b}$')
         ax_list[3].set_ylabel(r'$\tau_{conv}$ (s)', color='blue')
         #ax_list[3].plot(t_slice, tnuc, color='green', label=r'$\tau_{nuc}$')
         ax_list[3].set_title(sim.getLabel() + ' | Timescales')
         ax_list[3].set_xlabel(r'time [s]')
         ax_list[3].grid(which='both')
         tw2 = ax_list[3].twinx()
         tw2.plot(t_slice, ratio, color='cyan', label=r'$\tau_{nuc,x}/\tau_{conv}$')
         tw2.set_ylabel(r'$\tau_{nuc,x}/\tau_{conv,b}$ ratio', color='cyan')

         #tw2.set_ylim(200.0, 350.0)
         #tw2.set_xlim(200.0, 225.0)
         #tw2.set_xlim(60.0, 70.0)
         #tw2.set_xlim(145.0, 155.0)
         ax_list[3].plot(t_slice, tnuc_x, color='green', label=r'$\tau_{nuc,x}$')
         ax_list[3].plot(t_slice, tnuc_xb, color='green', label=r'$\tau_{nuc,xb}$')
         ax_list[3].plot(t_slice, tnuc_wk, color='green', linestyle='--', label=r'$\tau_{nuc,wk}$')
         ax_list[3].set_ylabel(r'$\tau_{nuc}$ (s)', color='green')
         ax_list[3].legend()

         #Plot base density (density at radius of peak temp)
         ax_list[4].plot(t_slice, avgBaseRho)
         ax_list[4].plot(t_slice, rhoCrit, 'b--')
         #ax_list[4].set_yscale('log')
         #ax_list[4].set_ylim(9.e5, 1.5e6)
         ax_list[4].set_ylabel(r'$\rho$ [g cm$^{-3}$]')
         ax_list[4].set_title(sim.getLabel() + ' | Base Density')

         #Set plot properties 
         fig.set_size_inches(10.0, 40.0)
         #fig.tight_layout()
         #fig.savefig("convt.png", bbox_inches='tight')
      else:
         #TODO: Need to implement non-inline plotting
         pass


   def plotEntropyDetail(self, sim, time, r_lims=(3.65e8, 4.25e8), ent_lims=(1.e8, 1.e9), ds_lims=(-50000.0, 50000.0)):
      """Plot entropy with ability to focus in on details"""
      import matplotlib.pyplot as plt
      import numpy as np
      from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
      from mpl_toolkits.axes_grid1.inset_locator import mark_inset
  
      #NOTE: fontsize can be [size in points | 'xx-small' | 'x-small' |
      #  'small' | 'medium' | 'large' | 'x-large' | 'xx-large'
      fontsize = 'x-large'
  
      # make the plots
      fig, ax = plt.subplots(nrows=1, ncols=1)

      #Find the right slice, if requested
      slicedata = None
      if time is not None:
         for s in sorted(sim._out._slices, key=lambda sl: sl.timestamp[0]):
            if s.timestamp[1] > time:
               slicedata = s
               break

      #Plot the entropy subplot
      self._plotSDetail(ax, sim, slicedata, fontsize, r_lims, ent_lims, ds_lims)

      #Adjust figure settings
      fig.set_size_inches(8.0,6.0)
      fig.tight_layout()
  
      #Write header for interactive viewing
      fig.suptitle(sim._label + ' Slice', fontsize=18, y=1.00)


   ##Private methods##
   def _plotRhoT(self, sp, sim, slicedata, fontsize):
      """Plot a density,temperature vs. radius subplot sp. If slicedata is given, use it,
      if None use initial model data."""
      import pylab
      import numpy as np

      # grab data
      if slicedata is None:
         model_data = sim.getIMDict()
         r    = model_data['radius']
         rho  = model_data['density']
         temp = model_data['temperature']
      else:
         r    = slicedata.datalist[SCSlice.IRADIUS]
         rho  = slicedata.datalist[SCSlice.IRHO]
         temp = slicedata.datalist[SCSlice.ITEMP]
  
      # plot density
      sp.plot(r, rho, color="blue")
  
      # plot the sponge start, anelastic cutoff, and base cutoff density
      cuts = sim.getIMCuts()
      lims = sim.getIMLims()
      sp.plot([cuts.r_an, cuts.r_an ], [lims.dlims[0], 10.0*lims.dlims[1]], color="k", linestyle='--')
      sp.plot([cuts.r_sp, cuts.r_sp ], [lims.dlims[0], 10.0*lims.dlims[1]], color="k", linestyle='--')
      sp.plot([cuts.r_bc, cuts.r_bc ], [lims.dlims[0], 10.0*lims.dlims[1]], color="k", linestyle='--')

      # shade convection zone
      conv_bounds = sim.getIMConvBounds()
      lidx = np.where(r >= conv_bounds[0])[0][0]
      ridx = np.where(r <= conv_bounds[1])[0][-1]
      r_conv = r[lidx:ridx]
      sp.fill_between(r_conv, lims.dlims[0], 10.0*lims.dlims[1], facecolor='0.9', edgecolor='1.0')

      # set plot properties
      # NOTE: the offset text is the 1.e8 that appears under the axis labels when
      # doing scientific notation
      sp.set_ylabel(r"density (g cm$^{-3}$)", color="blue", fontsize=fontsize)
      sp.set_yscale('log')
      sp.set_xlim(lims.rlims[0], lims.rlims[1])
      sp.set_ylim(lims.dlims[0], lims.dlims[1])
      sp.tick_params(labelbottom='off', labelsize=fontsize)
      sp.xaxis.offsetText.set_visible(False)
  
      # plot temperature on a twin axis
      sp2 = sp.twinx()
      sp2.plot(r, temp, color="red")
  
      # set plot properties
      sp2.yaxis.tick_right()
      sp2.set_yscale('log')
      sp2.axis(labelcolor="red")
      sp2.tick_params(labelsize=fontsize)
      sp2.set_ylabel(r"temperature (K)", color="red", fontsize=fontsize)
      sp2.set_xlim(lims.rlims[0], lims.rlims[1])
      sp2.set_ylim(lims.Tlims[0], lims.Tlims[1])
  
      # plot zoomed inset
      fig = sp.get_figure()
      axins = fig.add_axes([0.2, 0.75, 0.30, 0.10])
      axins.plot(r, rho, color="blue")
      axins_twin = axins.twinx() 
      axins_twin.plot(r, temp, color="red")
      axins.plot([cuts.r_an, cuts.r_an ], [lims.dlims[0], 10.0*lims.dlims[1]], color="k", linestyle='--')
      axins.plot([cuts.r_sp, cuts.r_sp ], [lims.dlims[0], 10.0*lims.dlims[1]], color="k", linestyle='--')
      axins.plot([cuts.r_bc, cuts.r_bc ], [lims.dlims[0], 10.0*lims.dlims[1]], color="k", linestyle='--')
  
      # set inset properties
      axins.set_yscale('log')
      axins_twin.set_yscale('log')
      axins.set_xlim(lims.zbounds[0][0], lims.zbounds[0][1])
      axins.set_ylim(lims.zbounds[1][0], lims.zbounds[1][1])
      axins.xaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))
      axins.set_yticklabels([' '], visible=False)
      axins_twin.set_ylim(lims.zbounds[2][0], lims.zbounds[2][1])
      axins_twin.set_yticklabels([' '], visible=False)

      return

   def _plotX(self, sp, sim, slicedata, fontsize):
      """Plot a mass fraction vs. radius subplot sp. If slicedata is given, use it,
      if None use initial model data."""
      import pylab
      import numpy as np

      # grab data
      if slicedata is None:
         model_data = sim.getIMDict()
         r    = model_data['radius']
         xhe  = model_data['species']['helium-4']
         xc   = model_data['species']['carbon-12']
      else:
         r    = slicedata.datalist[SCSlice.IRADIUS]
         xhe  = slicedata.datalist[SCSlice.ISPEC]['X(He4)'][0]
         xc   = slicedata.datalist[SCSlice.ISPEC]['X(C12)'][0]
      
      # Plot mass fractions and a legend
      sp.plot(r, xhe, color="red",  label=r"$^{4}\mathrm{He}$")
      sp.plot(r, xc,  color="blue", label=r"$^{12}\mathrm{C}$")
      sp.legend(loc=2, fontsize=fontsize)
  
      # shade convection zone
      conv_bounds = sim.getIMConvBounds()
      lidx = np.where(r >= conv_bounds[0])[0][0]
      ridx = np.where(r <= conv_bounds[1])[0][-1]
      r_conv = r[lidx:ridx]
      sp.fill_between(r_conv, 0, 1.05, facecolor='0.9', edgecolor='1.0')

      # draw in the sponge star, anelastic cutoff, and base cutoff density
      cuts = sim.getIMCuts()
      lims = sim.getIMLims()
      sp.plot([cuts.r_an, cuts.r_an], [0, 1.05], color="k", linestyle='--')
      sp.plot([cuts.r_sp, cuts.r_sp], [0, 1.05], color="k", linestyle='--')
      sp.plot([cuts.r_bc, cuts.r_bc], [0, 1.05], color="k", linestyle='--')

      # set plot properties
      sp.set_ylabel("mass fraction", fontsize=fontsize)
      sp.xaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))
      sp.tick_params(labelbottom='off', labelsize=fontsize)
      sp.xaxis.offsetText.set_visible(False)
      sp.set_xlim(lims.rlims[0], lims.rlims[1])
      sp.set_ylim(0.0,1.05)

   def _plotCs(self, sp, sim, slicedata, fontsize):
      """Plot a soundspeed vs. radius subplot sp. If slicedata is given, use it,
      if None use initial model data."""

      # grab data
      if slicedata is None:
         # Make sure we have soundspeed data
         model_data = sim.getIMDict()
         if not 'cs' in model_data.keys():
            raise ValueError("No soundspeed data loaded! Can't plot.")
         r    = model_data['radius']
         cs   = model_data['cs']
      else:
         raise NotImplementedError("Plotting soundspeed from slicedata not implemented.")
 
      sp.set_yscale('log')
      sp.plot(r, cs, color="red")
  
      # draw in the sponge start, anelastic cutoff, and base cutoff density
      cuts = sim.getIMCuts()
      lims = sim.getIMLims()
      sp.plot([cuts.r_an, cuts.r_an], [lims.clims[0], lims.clims[1]], color="0.65", linestyle='--')
      sp.plot([cuts.r_sp, cuts.r_sp], [lims.clims[0], lims.clims[1]], color="0.65", linestyle='--')
      sp.plot([cuts.r_bc, cuts.r_bc], [lims.clims[0], lims.clims[1]], color="0.65", linestyle='--')

      sp.set_xlabel("radius (cm)", fontsize=fontsize)
      sp.set_ylabel("sound speed (cm/s)", fontsize=fontsize)

      sp.xaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))
      sp.tick_params(labelsize=fontsize)
      sp.xaxis.offsetText.set_size(fontsize)
      #sp.yaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))

      sp.set_xlim(lims.rlims[0], lims.rlims[1])
      sp.set_ylim(lims.clims[0], lims.clims[1])

   def _plotS(self, sp, sim, slicedata, fontsize):
      """Plot entropy vs. radius subplot sp. If slicedata is given, use it,
      if None use initial model data."""
      import pylab
      import numpy as np

      # grab data
      if slicedata is None:
         # Make sure we have entropy data
         model_data = sim.getIMDict()
         if not 'entropy' in model_data.keys():
            raise ValueError("No entropy data loaded! Can't plot.")
         r = model_data['radius']
         S = model_data['entropy']
      else:
         r = slicedata.datalist[SCSlice.IRADIUS]
         S = slicedata.datalist[SCSlice.IS]

      #Plot entropy
      sp.plot(r, S, color="red")
  
      # draw in the sponge start, anelastic cutoff, and base cutoff density
      cuts = sim.getIMCuts()
      lims = sim.getIMLims()
      sp.plot([cuts.r_sp, cuts.r_sp], [lims.slims[0], lims.slims[1]], color="k", linestyle='--')
      sp.plot([cuts.r_an, cuts.r_an], [lims.slims[0], lims.slims[1]], color="k", linestyle='--')
      sp.plot([cuts.r_bc, cuts.r_bc], [lims.slims[0], lims.slims[1]], color="k", linestyle='--')

      # shade convection zone
      conv_bounds = sim.getIMConvBounds()
      lidx = np.where(r >= conv_bounds[0])[0][0]
      ridx = np.where(r <= conv_bounds[1])[0][-1]
      r_conv = r[lidx:ridx]
      sp.fill_between(r_conv, lims.slims[0], lims.slims[1], facecolor='0.9', edgecolor='1.0')

      # set plot properties
      sp.set_yscale('log')
      sp.set_xlabel("radius (cm)", fontsize=fontsize)
      sp.set_ylabel(r"specific entropy (erg g$^{-1}$ K$^{-1}$)", fontsize=fontsize)
      sp.xaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))
      sp.tick_params(labelsize=fontsize)
      sp.xaxis.offsetText.set_size(fontsize)
      sp.set_xlim(lims.rlims[0], lims.rlims[1])
      sp.set_ylim(lims.slims[0], lims.slims[1])

   def _plotSDetail(self, sp, sim, slicedata, fontsize, r_lims, ent_lims, ds_lims):
      """Plot entropy vs. radius"""
      import pylab
      import numpy as np

      # grab data
      r = slicedata.datalist[SCSlice.IRADIUS]
      S = slicedata.datalist[SCSlice.IS]
      P = slicedata.datalist[SCSlice.IP]
      S_SMOOTH = [(S[i+1] + S[i] + S[i-1])/3. for i in range(1,len(S)-1)]
      S_SMOOTH.insert(0,S_SMOOTH[0])
      S_SMOOTH.append(S_SMOOTH[-1])
      DS = [S_SMOOTH[i] - S_SMOOTH[i-1] for i in range(1,len(S_SMOOTH))]
      DS.insert(0,DS[0])

      #Plot entropy
      #sp.plot(r, S, 'rx') #color="red")
      sp.plot(r, P, 'rx') #color="red")
      #s_lo = ent_lims[0] #2.930e8
      #s_hi = ent_lims[1] #2.940e8
      #sp.set_xlim(r_lims[0], r_lims[1])
      #sp.set_ylim(s_lo, s_hi)
      #ent_level = slicedata.derived_datalist[SCSlice.IRCONV]
      #sp.plot([ent_level, ent_level], [s_lo, s_hi], color="k", linestyle='--')

      #On same axis, plot ds
      #tw = sp.twinx()
      #tw.plot(r, DS, 'bx')#color="blue")
      #tw.set_ylim(ds_lims[0], ds_lims[1])
      #tw.set_xlim(r_lims[0], r_lims[1])
  
      # set plot properties
      #sp.set_yscale('log')
      sp.set_xlabel("radius (cm)", fontsize=fontsize)
      sp.set_ylabel(r"specific entropy (erg g$^{-1}$ K$^{-1}$)", fontsize=fontsize)
      sp.xaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))
      sp.tick_params(labelsize=fontsize)
      sp.xaxis.offsetText.set_size(fontsize)

   def _labelRow(self, dirlist):
      """Return string corresponding to html table row describing labels of dirs in dirlist."""
      ret = r'  <tr>' + "\n"
      ret += r'    <th>Model Label</th>' + "\n"
      for curdir in dirlist:
         ret += r'    <th>' + curdir + r'</th>' + "\n"
      ret += r'  </tr>' + "\n"
      return ret
   
   def _massRow(self, dirlist, stg_dir):
      """Return string corresponding to html table row describing mass configuration of runs in dirlist."""
      from glob import glob
  
      ret = r'  <tr>' + "\n"
      ret += r'    <th>(Core, Shell) Mass [M$_\odot$]</th>' + "\n"
      for curdir in dirlist:
         pfile = glob(stg_dir + '/' + curdir + '/run/_params.*') #Grab first parameter file found, should only be one
         m_core = get_param('M_tot', pfile[0])
         m_shell = get_param('M_He', pfile[0])
         mstr = '(' + str(m_core) + ', ' + str(m_shell) + ')'
         ret += r'    <td>' + mstr + r'</td>' + "\n"
      ret += r'  </tr>' + "\n"
      return ret
   
   def _paramRow(self, dirlist, stg_dir, label, param):
      """Return string corresponding to html table row with the given label and parameter for all runs in dirlist."""
      from glob import glob
      ret = r'  <tr>' + "\n"
      ret += r'    <th>' + label + r'</th>' + "\n"
      for curdir in dirlist:
         pfile = glob(stg_dir + '/' + curdir + '/run/_params.*') #Grab first parameter file found, should only be one
         pstr = get_param(param, pfile[0])
         ret += r'    <td>' + pstr + r'</td>' + "\n"
      ret += r'  </tr>' + "\n"
      return ret
   
   def _inputsRow(self, dirlist, stg_dir, label, param):
      """Return string corresponding to html table row with the given label and input parameter for all runs in dirlist."""
      from glob import glob
      ret = r'  <tr>' + "\n"
      ret += r'    <th>' + label + r'</th>' + "\n"
      for curdir in dirlist:
         pfile = glob(stg_dir + '/' + curdir + '/run/inputs*') #Grab first inputs file found, should only be one
         pstr = get_param(param, pfile[0])
         ret += r'    <td>' + pstr + r'</td>' + "\n"
      ret += r'  </tr>' + "\n"
      return ret
   
   def _dxRow(self, dirlist, stg_dir, label):
      """Return string corresponding to html table row with the given label and the fine dx resolution for all runs in dirlist."""
      from glob import glob
      ret = r'  <tr>' + "\n"
      ret += r'    <th>' + label + r'</th>' + "\n"
      for curdir in dirlist:
         pfile = glob(stg_dir + '/' + curdir + '/run/inputs*') #Grab first inputs file found, should only be one
         # get max_levs, n_cellx, and xmax
         lev = int(get_param('max_levs', pfile[0]))
         nx = int(get_param('n_cellx', pfile[0]))
         xmax = float(get_param('prob_hi_x', pfile[0]).replace('d','e'))
         dxf = xmax/float(nx*2**(lev-1))
         dxf = dxf / 1.e5 #Convert to km
         ret += r'    <td>' + str(dxf) + r'</td>' + "\n"
      ret += r'  </tr>' + "\n"
      return ret
   
   def _rhobaseRow(self, dirlist, stg_dir, label):
      """Return string corresponding to html table row with the given label and 
      the density at the base of the convective layer (where T peaks)."""
      from glob import glob
      import numpy as np
      ret = r'  <tr>' + "\n"
      ret += r'    <th>' + label + r'</th>' + "\n"
      for curdir in dirlist:
         imfile = glob(stg_dir + '/' + curdir + '/run/sub_chandra*hse*') #Grab first hse initial model file, should only be one
         # find peak T, store corresponding density
         rho, T = np.loadtxt(imfile[0], usecols=(1, 2), unpack=True)
         idx = T.argmax()
         rbase = rho[idx]
         rbase = rbase / 1.e5
   
         ret += r'    <td>' + str(rbase) + r'</td>' + "\n"
      ret += r'  </tr>' + "\n"
      return ret


class SCSlice(object):
   """A class representing a 1D radial slice of angle-averaged data from a
   sub-Chandra run's plotfile."""
   ##Class/pseudo-static data, constructor##
   #Indices for the datalist 
   #TODO:I don't like explicitly manipulating so many quantities, but not sure what a better alternative is
   IRADIUS   = 0
   IRHO      = 1  #Density
   IRHO_RMS  = 2
   ITEMP     = 3  #Temperature
   ITEMP_RMS = 4  
   IP        = 5  #Pressure
   IP_RMS    = 6
   IMV       = 7  #Magnitude of Velocity vector U
   IHNUC     = 8  #\sum_k{q_k * omegadot_k}
   IS        = 9  #Entropy
   IS_RMS    = 10
   ISPEC     = 11 #Species, mass fractions and their time derivatives

   #Indices for derived_datalist
   IINTFACE = 0 #CO/He Interface radius (radius where X_He = 0.9)
   IRCONV   = 1 #Radius of the bottom of the convective region
   ILCONV   = 2 #Lengthscale of the convective region.  I define the lengthscale
                #as the region in which ds/dr is numerically <= 0 (s is specific entropy)
   IUCONV   = 3 #Average velocity magnitude in the convective region
   IPSCALE  = 4 #Length of a pressure scale height.  From the radius of
                #T_peak to the radius where the pressure at that point has
                #fallen by 1/e
                #TODO: Compare this to common equations used to determine H
   ITNUC    = 5 #Tuple of radial arrays of different nuclear timescales
                #  (X/dXdt, WK 2011 eqn 3)

   #Size of the non-species part of datalist
   _NS_SIZE = 11 

   #Constructor
   def __init__(self, slicefile, parent):
      """self      --> implicitly passed reference to this instance of SCSlice
         slicefile --> a .slice file generated by fsubchandra.f90
         parent    --> reference to this slice's parent SCOutput"""
      from os.path import basename

      (self.timestamp, self.nspec, self.glbs, 
       self.datalist) = self._parseSlicefile(slicefile)
      self.derived_datalist = self._deriveData()
      self.parent = parent
      #ASSUMPTION: Slicefiles are of form '<pltfile name>.slice'
      self.my_pltfile = basename(slicefile[:-6])

   ##Public methods##

   ##Private methods##
   def _deriveData(self):
      """Derive interesting data from slice data: 
      pressure scale height, interface radius, convection radii, enuc timescales."""
      import numpy as np
      #The DS critical values are based on analyzing models at different extremes and
      #are found to accurately demark the convective zone where ds/dr <= 0.
      DS_TRIGGER = 1.e6  #This triggers the search for the convective zone, 
                         #it must be after this point
      DS_THRESH = 2.0e4  #When DS falls below this, we're in the convectively unstable zone, 
                         #when it rises back above it we're out of the convective zone
      ENT_TOL = 2.0e5    #Tolerance for distance from entropy level, 
                         #beyond which we consider to not be convective
      ENT_TOL = 0.001
      TNUC_MAX = 5.e15

      retlist = [None, None, None, None, None, None]
      tpeak = self.datalist[SCSlice.ITEMP].max()
      ppeak = None
      l1 = None
      l2 = None
      vavg = 0.0
      vsum = 0.0
      vcnt = 0
      #tnuc = (tnuc based on X/dXdt, tnuc based on W&K 2011 eqn 3)
      tnuc = (np.zeros(len(self.datalist[SCSlice.IRADIUS])), np.zeros(len(self.datalist[SCSlice.IRADIUS])), np.zeros(len(self.datalist[SCSlice.IRADIUS])))

      #These are the specific binding energies of [He4, C12, O16] from
      #AstroDev/networks/triple_alpha_plus_cago/network.f90 
      #in units of erg/g
      q = [6.8253797e18, 7.4103097e18, 7.6959581e18]
      q_tot = sum(q)
      #E_CRIT = 1.7e-7  # erg
      #TEST_VOL = 4./3.*np.pi*(1e7)**3 # cm^3

      trigger = False
      trigger_count = 0
      min_R = 2.0e8
      ent_level = 0.0
      ent_count = 0
      calc_ent_level = True
      for i in range(len(self.datalist[SCSlice.IRADIUS])):
         #Calculate estimate of interface radius based on when the helium mass fraction hits 0.9
         if not retlist[SCSlice.IINTFACE] and self.datalist[SCSlice.ISPEC]['X(He4)'][0][i] > 0.9:
            retlist[SCSlice.IINTFACE] = self.datalist[SCSlice.IRADIUS][i]

         #Calculate the entropy level.  This is the average value of the specific entropy
         #in the region where its radial profile flattens out, which is one criteria for
         #determining the convective region
         if calc_ent_level and i > 2 and i < len(self.datalist[SCSlice.IRADIUS]):
            #Calculate smoothed entropy difference, ds
            s_smooth1 = (self.datalist[SCSlice.IS][i-1] + self.datalist[SCSlice.IS][i] + self.datalist[SCSlice.IS][i+1])/3.0
            s_smooth2 = (self.datalist[SCSlice.IS][i] + self.datalist[SCSlice.IS][i+1] + self.datalist[SCSlice.IS][i+2])/3.0
            ds = s_smooth2 - s_smooth1

            #Trigger the search for the entropy level
            if ds > DS_TRIGGER and self.datalist[SCSlice.IRADIUS][i] > min_R:
               trigger_count += 1
               if trigger_count > 3:
                  #To avoid a false trigger due to noise, 
                  #make sure we've been above the trigger consistently for the
                  #last few elements in the array.  If yes, then trigger away.
                  trigger = True
            else:
               trigger_count = 0
            
            #If the entropy trigger was set off and we're below the threshold
            #we've found the bottom of the convective envelope
            if trigger and ds < DS_THRESH:
               ent_level += self.datalist[SCSlice.IS][i]
               ent_count += 1
           
            #We're outside of the convecting region, stop calculation of ent_level
            if trigger and ent_level > 0.0 and ds >= DS_THRESH:
               calc_ent_level = False

         #Calculate nuclear timescale with limits on how large it can get
         for species in self.datalist[SCSlice.ISPEC]:
            x = self.datalist[SCSlice.ISPEC][species][0][i]
            dxdt = self.datalist[SCSlice.ISPEC][species][1][i]
            if dxdt < 0.0: #If negative, then fuel
               dxdt = -dxdt
               #Normalize with x
               if tnuc[0][i] == 0.0:
                  tnuc[0][i] = x/dxdt
               elif x/dxdt < tnuc[0][i]:
                  tnuc[0][i] = x/dxdt
               if tnuc[0][i] > TNUC_MAX: 
                  tnuc[0][i] = TNUC_MAX
               #Raw dxdt
               if tnuc[2][i] == 0.0:
                  tnuc[2][i] = 1.0/dxdt
               elif 1.0/dxdt < tnuc[2][i]:
                  tnuc[2][i] = 1.0/dxdt
               if tnuc[2][i] > TNUC_MAX: 
                  tnuc[2][i] = TNUC_MAX
         if tnuc[0][i] == 0.0:
            tnuc[0][i] = TNUC_MAX
         if tnuc[2][i] == 0.0:
            tnuc[2][i] = TNUC_MAX

         rho6 = self.datalist[SCSlice.IRHO][i]/1.e6
         T8 = self.datalist[SCSlice.ITEMP][i]/1.e8
         tnuc[1][i] = 3.4e-4 * np.exp(20./T8) * pow(rho6, -2.3)

      ent_level = ent_level / ent_count
      #print('ent_level: {}'.format(ent_level))
      top_i = None
      for i in range(len(self.datalist[SCSlice.IRADIUS])):
         rr  = self.datalist[SCSlice.IRADIUS][i]
         ent_norm = self.datalist[SCSlice.IS][i]/ent_level

         if np.abs(ent_norm - 1.0) < ENT_TOL:
            #We're in the convecting region, calculate various quantities of interest

            #First, get the base of this region
            if not l1:
               l1 = self.datalist[SCSlice.IRADIUS][i]
               ppeak = self.datalist[SCSlice.IP][i]
               retlist[SCSlice.IRCONV] = l1
            #Calculate average velocity magnitude in convective region
            vavg += self.datalist[SCSlice.IMV][i] 
            vcnt += 1
            vsum += 1.0 / self.datalist[SCSlice.IMV][i]
            
            #Get the index of the top of the convecting region
            #The last value assigned to top_i will be this index
            top_i = i

         #Calculate pressure scale height
         if ppeak and not retlist[SCSlice.IPSCALE]: #ppeak has been found and we haven't found the scale height yet
            if self.datalist[SCSlice.IP][i] <= ppeak/np.e: #Found the scale height
               retlist[SCSlice.IPSCALE] = rr - l1

      l2 = self.datalist[SCSlice.IRADIUS][top_i]
      retlist[SCSlice.ILCONV] = l2 - l1
      retlist[SCSlice.IUCONV] = vavg / vcnt
      #dr = self.datalist[SCSlice.IRADIUS][2] - self.datalist[SCSlice.IRADIUS][1]
      #retlist[SCSlice.IUCONV] = dr*vsum
      retlist[SCSlice.ITNUC] = tnuc

      return retlist

   def _parseSlicefile(self, slicefile):
      """Return (timestamp tuple, nspec, globals, data list) based on slicefile."""
      import numpy as np
      import re
      SPEC_COL = 7

      sf = open(slicefile)
      tret = (None, None) #(time step, time in seconds)
      nspec = None
      specmap = {}
      gret = _Globals()
      dlret = [None for i in range(SCSlice._NS_SIZE+1)] #Use list comprehension
  
      #Parse header / global info
      #TODO: I make strict assumptions about the structure of the header,
      #might be nice to develop a more portable/generic algorithm
      for line in sf:
         if not line.strip():
            continue #If blank line, just loop to next line
         if not line.strip().startswith('#') :
            break #Reached end of header

         if line.count('time') == 1:
            tokens = line.split()
            time = float( tokens[len(tokens)-1] )
         elif line.count('peak temperature') == 1:
            tokens = line.split()
            gret.Tpeak = float( tokens[len(tokens)-1] )
         elif line.count('peak temp loc') == 1:
            tokens = line.split()
            tcart = tuple( float(val) for val in tokens[-3:]  )
         elif line.count('peak temp radius') == 1:
            tokens = line.split()
            trad = float( tokens[len(tokens)-1] )
         elif line.count('velocity @ peak T loc (vx') == 1:
            tokens = line.split()
            tvcart = tuple( float(val) for val in tokens[-3:]  )
         elif line.count('radial velocity @ peak T') == 1:
            tokens = line.split()
            tvrad = float( tokens[len(tokens)-1] )
         elif line.count('peak enucdot =') == 1:
            tokens = line.split()
            gret.epeak = float( tokens[len(tokens)-1] ) 
         elif line.count('peak enucdot loc') == 1:
            tokens = line.split()
            ecart = tuple( float(val) for val in tokens[-3:]  )
         elif line.count('peak enucdot radius') == 1:
            tokens = line.split()
            erad = float( tokens[len(tokens)-1] )
         elif line.count('pressure') == 2:
            #This line labels the data, use it to calculate nspec and get species names
            #Split with regex representing two or more spaces -- some labels have one space in them 
            tokens = re.split(r'\s{2,}', line) 
            nspec = len(tokens) - SCSlice._NS_SIZE - 1 #-1 for the comment character
            nspec = nspec/2 #Divide by 2 because we have both X and dX/dt for each species
            speclist = [spec for spec in tokens[SPEC_COL+1:SPEC_COL+nspec+1]] #+1 for comment char

      sf.close()

      gret.tploc = tcart + (trad,)
      gret.tpvel = tvcart + (tvrad,)
      gret.eploc = ecart + (erad,)

      #Record time info
      li = slicefile.index('_plt') + 4
      ri = slicefile.index('.slice') 
      step = int(slicefile[li:ri])
      tret = (step, time)

      #Read in non-species data.  I know, this is hideous.
      rs = 7+2*nspec #RMS start
      (dlret[SCSlice.IRADIUS], dlret[SCSlice.IRHO], dlret[SCSlice.IRHO_RMS], 
       dlret[SCSlice.ITEMP], dlret[SCSlice.ITEMP_RMS], dlret[SCSlice.IP], 
       dlret[SCSlice.IP_RMS], dlret[SCSlice.IMV], dlret[SCSlice.IHNUC],
       dlret[SCSlice.IS], dlret[SCSlice.IS_RMS]) = np.loadtxt(slicefile, 
                 usecols=(0,1,rs,2,rs+1,3,rs+2,4,5,6,rs+3), unpack=True)

      #Read species data
      #Species maps 'X(<species>)' --> (<array of X(r)>, <array of dX/dt(r)>)
      for (i, key) in enumerate(speclist):
         specmap[key] = (np.loadtxt(slicefile, usecols=(SPEC_COL+i,)), np.loadtxt(slicefile, usecols=(SPEC_COL+nspec+i,)))
      dlret[SCSlice.ISPEC] = specmap 

      return (tret, nspec, gret, dlret)

   def _write3DRVIn(self, pltfile, min_rv, conv_bot, conv_top, H, rmax, r_anelastic, r_sp_start, savefile):
      """Write the input needed for the VisIt plotting script."""
      from numpy import sqrt
      #TODO: If the unused args prove not helpful, clean up

      PLOTRV_IN = '/ccs/home/ajacobs/Projects/SubChandra/lensSubmit/plotRV3D.in'

      with open(PLOTRV_IN, 'w') as prv_in:
         #Write header
         prv_in.write('#Inputs for the plotRV3D VisIt script\n\n')

         #Write pltfile's location
         prv_in.write('#Database/file to be read\n')
         prv_in.write('pltfile = {0}\n\n'.format(pltfile))

         #Convective region
         #prv_in.write('#The radius of the bottom and top of the convective layer\n')
         #prv_in.write('#(where entropy is flat, ds/dr <=0) in cm\n')
         #prv_in.write('#Radius of conv_bot + 1 pressure scale height\n')
         #prv_in.write('conv_bot = {0}\n'.format(conv_bot))
         #prv_in.write('conv_top = {0}\n'.format(conv_top))
         #prv_in.write('H = {0}\n\n'.format(H))

         #Maximum x
         prv_in.write('#Maximum x (same as max y, z)\n')
         xmax = rmax / sqrt(3.0)
         prv_in.write('xmax = {0}\n\n'.format(xmax))

         #Minimum radial velocity to plot
         prv_in.write('#Minimum radial velocity to plot\n')
         prv_in.write('min_rv = {0}\n\n'.format(min_rv))

         #Radii of Maestro params
         #prv_in.write('#Radii of the anelastic cutoff density and start of the sponge\n')
         #prv_in.write('r_anelastic = {0}\n'.format(r_anelastic))
         #prv_in.write('r_sp_start = {0}\n\n'.format(r_sp_start))

         #File to save to
         prv_in.write('#File to save the plot to\n')
         prv_in.write('savefile = {0}\n'.format(savefile))

   def _writeRVIn(self, infilename, pltfile, conv_bot, conv_top, H, rmax, r_anelastic, r_sp_start, savefile):
      """Write the input needed for the VisIt plotting script."""
      from numpy import sqrt

      with open(infilename, 'w') as prv_in:
         #Write header
         prv_in.write('#Inputs for the plotRV VisIt script\n\n')

         #Write pltfile's location
         prv_in.write('#Database/file to be read\n')
         prv_in.write('pltfile = {0}\n\n'.format(pltfile))

         #Convective region
         prv_in.write('#The radius of the bottom and top of the convective layer\n')
         prv_in.write('#(where entropy is flat, ds/dr <=0) in cm\n')
         prv_in.write('#Radius of conv_bot + 1 pressure scale height\n')
         prv_in.write('conv_bot = {0}\n'.format(conv_bot))
         prv_in.write('conv_top = {0}\n'.format(conv_top))
         prv_in.write('H = {0}\n\n'.format(H))

         #Maximum radius
         prv_in.write('#Maximum x (same as max y, z)\n')
         xmax = rmax / sqrt(3.0)
         prv_in.write('xmax = {0}\n\n'.format(xmax))

         #Radii of Maestro params
         prv_in.write('#Radii of the anelastic cutoff density and start of the sponge\n')
         prv_in.write('r_anelastic = {0}\n'.format(r_anelastic))
         prv_in.write('r_sp_start = {0}\n\n'.format(r_sp_start))

         #File to save to
         prv_in.write('#File to save the plot to\n')
         prv_in.write('savefile = {0}\n'.format(savefile))




#################
### Functions ###
#################

#TODO: Would be better to have objects with data dictionaries for this sort of thing.
#      Should rework code this way
def get_param(param, param_file):
   """Return the parameter value (as a string) found in a simple parameter 
      file with <param> = <val> assignments.  None is returned if not found."""
   pfile = open(param_file)
   ret = None
   for line in pfile:
      if(line.find('=') > -1):
         tokens = line.partition('=')
         cur_param = tokens[0].strip()
         cur_val = tokens[2]
         if(cur_param == param):
            ret = cur_val
   pfile.close()
   return ret

def found_required_scratch_files(rundir):
   """Check that all files required for a sub-Chandra simulation are present in rundir"""
   from glob import glob
   from os.path import join, isdir, isfile

   #Make sure either valid chk file or initial model data is available
   if not isfile(join(rundir, 'sub_chandra.M_WD*hse*')):
      #No initial model, better have a chkfile
      chkfiles = []
      for f in glob(join(rundir, '*_chk*')):
         if isdir(f):
            chkfiles.append(f)
      chkfiles.sort()
      last_chk_file = chkfiles[-1]
      if not isfile(join(last_chk_file, 'Header')):
         return False

   #Make sure we have the executable, inputs file, and EoS data
   return (len(glob( join(rundir, 'main.*.exe') ))     > 0 and
           len(glob( join(rundir, 'inputs*') ))        > 0 and
           len(glob( join(rundir, 'helm_table.dat') )) > 0
          )


#################
### Execution ###
#################
#This is only for testing.  subchandra.py is intended to be used as a module.
if __name__== "__main__":
   if(len(sys.argv) <= 1):
      #TODO: Add any desired testing
      pass


######################
### OLD!!!! DELETE ###
######################

class _Limits(object):
   """struct-like class for containing plot limits."""
   #Axis limits
   Tlims = (None, None)
   rlims = (None, None)
   dlims = (None, None)
   clims = (None, None)
   slims = (None, None)

   #Zoomed inset limits
   rzoom = (None, None)
   dzoom = (None, None)
   Tzoom = (None, None)
   zbounds = (rzoom, dzoom, Tzoom)

class _Cutoffs(object):
   """struct-like class for containing cutoffs."""
   an_cut = None
   sp_cen_den = None
   sp_st_fac = None
   base_cut = None

class SCOutput(object):
   """A class representing the output of a particular sub-Chandra simulation."""
   ##Shared class data##

   ##Constructor##
   def __init__(self, stage_dir, scratch_dir, label, parent):
      """self        --> implicitly passed reference to this instance of SCSimulation
         stage_dir   --> staging directory for this sub-Chandra simulation
         scratch_dir --> scratch directory where the work is done but data's purged
         label       --> this simulation's label (e.g. 12050-107-175-3levs)
         parent      --> the SCSimulation parent of this SCOutput"""
      from glob import glob
      from os.path import join
      self._stage_dir = stage_dir.rstrip('/') #Get rid of any trailing '/'
      self._scratch_dir = scratch_dir.rstrip('/')
      self._label = label
      self.parent = parent

      #Load all available slicefiles
      filelist = glob(join(self._stage_dir, 'output', '*.slice'))
      self._slices = [SCSlice(sfile, self) for sfile in filelist]

      #Load all available hotspots files
      filelist = glob(join(self._stage_dir, 'output', '*.hotspots'))
      if parent is None:
         self._hotspots = [SCHotspots(hsfile) for hsfile in filelist]
      else:
         self._hotspots = [SCHotspots(hsfile, parent._inputs) for hsfile in filelist]

      #Load available temperature histogram files
      filelist = glob(join(self._stage_dir, 'output', '*.temphist'))
      self._thists = [SCTempHist(thfile, self) for thfile in filelist]

      #Load diagnostics data
      self._diags = SCDiagnostics(stage_dir, scratch_dir, label)

   ##Public methods##
   def getTHistData(self, step):
      """Return the temperature histogram object for the given step."""
      for th in self._thists:
         if th.timestamp[0] == step:
            return th
     
   def getPeakState(self, t_cut=None):
      """Get the temperature and time at the time of peak global temperature,
      as well as the slice with a timestamp nearest the peak temp time."""
      Tpeak, timepeak, rpeak = self._diags.getPeakState()

      #Find the slice with a timestamp nearest timepeak
      dt = 1.e9
      peakslice = None
      for s in self._slices:
         t = s.timestamp[1]
         if t_cut and t > t_cut: #Optional time cutoff
            break
         if (abs(t - timepeak) < dt):
            dt = abs(t - timepeak)
            peakslice = s

      return Tpeak, timepeak, rpeak, peakslice

   def getTurnover(self):
      """Get the average convective turnover time, skipping the first 50 seconds during which
      convection is being established."""
      import numpy as np

      tconv = []
      for s in self._slices:
         t = s.timestamp[1]
         if t < 50.:
            continue
         tconv.append(s.derived_datalist[SCSlice.ILCONV] / s.derived_datalist[SCSlice.IUCONV])
      tcavg = sum(tconv)/len(tconv)

      return tcavg


   def stageOutput(self):
      """Check the scratch directory for any updated output.  If found, copy to
      the stage directory."""
      #TODO: Implement this
      pass

   def analyzePlts(self, hpcnt=-1.0, full_star=False):
      """Run the Fortran analysis routine 'fsubchandra.f90' on all of this 
      simulation's plotfiles in the scratch directory. Copy the generated
      output to the staging directory's 'output' directory.
  
      hpcnt is optional.  It's the percentage of cells to track 
      when calculating the hottest cells for doing hotspot statistics. If <= 0.0 
      then no hotspot statistics are calculated."""
      from os import listdir
      from os.path import join, basename, dirname, isfile
      from re import match
      from subprocess import call
      from glob import glob
      from shutil import copy2
      FSUB_CMD = '/u/sciteam/ajacobs/Codebase/AmrPostprocessing/F_Src/MAESTRO_sub_chandra/fsubchandra.Linux.gfortran.exe'

      #Get a list of valid plotfiles from this simulation's 
      #base directory and plot directory
      pltfiles = self._getPltFiles()
  
      #Loop over plotfiles, analyzing each
      for plt in pltfiles:
         print('checking ', plt)
         outdir = join(self._stage_dir, 'output')
         pltfname = basename(plt)
         #Only do the analysis if it hasn't been done yet, 
         #and only calculate hotspot data if asked
         #TODO: For now I'm just assuming if the output exists it's current.
         #      Would be better to check file timestamps
         if(hpcnt > 0.0):
            found_slice = isfile(join(outdir, pltfname + '.slice'))
            found_hs = isfile(join(outdir, pltfname + '.hotspots'))
            if not (found_slice and found_hs):
               args = ['-h', str(hpcnt), plt + '/']
               if(full_star):
                 args.insert(0, '-f')
               print('Running fsubchandra, may take several minutes...')
               #Do a flush to make sure the print makes it to stdout before
               #call() starts blocking
               sys.stdout.flush()  
               call_list = []
               call_list.append(FSUB_CMD)
               for a in args:
                 call_list.append(a)
               call(call_list)

               #Copy the output files to the staging directory
               ofiles = [join(dirname(plt), pltfname + '.slice'), join(dirname(plt), pltfname + '.hotspots')]
               for out in ofiles:
                  dst = join(outdir, basename(out))
                  copy2(out, dst)
         else:
            found_slice = isfile(join(outdir, pltfname + '.slice'))
            if not found_slice:
               args = [plt + '/',]
               if(full_star):
                 args.insert(0, '-f')
               print('Running fsubchandra, may take several minutes...')
               #Do a flush to make sure the print makes it to stdout before
               #call() starts blocking
               sys.stdout.flush()  
               call_list = []
               call_list.append(FSUB_CMD)
               for a in args:
                 call_list.append(a)
               call(call_list)

               #Copy the output files to the staging directory
               slice_out = join(dirname(plt), pltfname + '.slice') 
               dst = join(outdir, basename(slice_out))
               copy2(slice_out, dst)

   def printTimestampOverview(self):
      """Print a summary of available timestamps with data."""
      print('{0:>4s}{1:>16s}{2:>16s}{3:>16s}{4:>16s}{5:>16s}'.format('Step', 'Time (s)', 'CO/He r (cm)', 'l_conv (cm)', 'U_conv (cm/s)', 'H (cm)'))
      print('{0:-<4s}{1:-<16s}{2:-<16s}{3:-<16s}{4:-<16s}{5:-<16s}'.format('-', '-', '-', '-', '-', '-'))
      for s in sorted(self._slices, key=lambda sl: sl.timestamp[0]):
         print('{0:>4d}{1:>16.4f}{2:>16.4g}{3:>16.4g}{4:>16.4g}{5:>16.4g}'.format(
               s.timestamp[0], s.timestamp[1], s.derived_datalist[SCSlice.IINTFACE],
               s.derived_datalist[SCSlice.ILCONV], s.derived_datalist[SCSlice.IUCONV],
               s.derived_datalist[SCSlice.IPSCALE]))

   def plotDiags(self):
      """Plot this simulation's diagnostics output: peak temperature, peak Mach #,
      etc... vs time."""
      import numpy as np
      import matplotlib.pyplot as plt
      import matplotlib

      #Get data
      t_T, Tpeak, Tpeak_r, Tpeak_vr  = self._diags._time, self._diags._maxT, self._diags._maxT_r, self._diags._maxT_vr
      t_M, Mpeak = self._diags._time, self._diags._maxmachfull

      if 'inline' in matplotlib.get_backend():
         #Build plots
         fig, ax_list = plt.subplots(nrows=2, ncols=1)

         #Temp
         ax_list[0].plot(t_T, Tpeak, color='red')
         ax_list[0].set_ylabel(r'T$_{\mathrm{peak}}$ (K)', color='red')
         ax_list[0].set_title(self._label + ' | Peak Temperature')
         tw = ax_list[0].twinx()
         tw.plot(t_T, Tpeak_r, color='green')
         tw.set_ylabel(r'T$_{\mathrm{peak}}$ radius (cm)', color='green')

         #Mach
         ax_list[1].plot(t_M, Mpeak, color='blue')
         ax_list[1].set_title(self._label + ' | Peak Mach #')
         ax_list[1].set_xlabel(r'time [s]')

         #Set plot properties 
         fig.set_size_inches(5.0, 8.0)
         fig.tight_layout()
      else:
         #TODO: Need to implement non-inline plotting
         pass

   def plotIgTimeFig(self, an_cut, sp_cen_den, sp_st_fac):
      """Plot a figure demonstrating ignition over time for publication, posters, presentation, etc..."""
      #TODO: I'm hacking this for a presentation, probably going to be ugly and need reworking
      import numpy as np
      import matplotlib.pyplot as plt
      import matplotlib
      from matplotlib.ticker import MaxNLocator
      import scipy.integrate as spint

      #Convenience aliases for diagnostic data
      t_T, Tpeak, Tpeak_r, Tpeak_vr  = self._diags._time, self._diags._maxT, self._diags._maxT_r, self._diags._maxT_vr
      t_M, Mpeak = self._diags._time, self._diags._maxmachfull
      t_e, epeak = self._diags._time, self._diags._maxenuc

      #Prepare variables for use in slice loop
      H = []
      iface = []
      rbot = []
      rtop = []
      t_slice = []

      tconv = []
      tconvb = []
      tnuc_x = []
      tnuc_xb = []
      tnuc_wk = []
      ratio = []

      avgTpeak = []
      avgBaseRho = []
      rhoCrit = []

      #Loop over slices in chronological order
      for s in sorted(self._slices, key=lambda sl: sl.timestamp[0]):
         #Build radius data from slices
         rbot.append(s.derived_datalist[SCSlice.IRCONV])
         rtop.append(s.derived_datalist[SCSlice.IRCONV] + s.derived_datalist[SCSlice.ILCONV])
         H.append(s.derived_datalist[SCSlice.IRCONV] + s.derived_datalist[SCSlice.IPSCALE])
         iface.append(s.derived_datalist[SCSlice.IINTFACE])
         t_slice.append(s.timestamp[1])

         #Estimate of convective turnover timescale and minimum nuclear timescale
         #tconv.append(s.derived_datalist[SCSlice.IPSCALE] / s.derived_datalist[SCSlice.IUCONV])
         tconv.append(s.derived_datalist[SCSlice.ILCONV] / s.derived_datalist[SCSlice.IUCONV])
         #tconv.append(s.derived_datalist[SCSlice.IUCONV])
         tnuc_x.append(min(s.derived_datalist[SCSlice.ITNUC][0]))
         tnuc_wk.append(min(s.derived_datalist[SCSlice.ITNUC][1]))
         tnuc_xb.append(min(s.derived_datalist[SCSlice.ITNUC][2]))
         #ratio.append(tnuc_x[len(tconv)-1]/tconv[len(tconv)-1])
         #ratio.append(tnuc_wk[len(tconv)-1]/tconv[len(tconv)-1])

         #Get the peak radially averaged temperature as an estimate of the background
         #conditions the hottest spot is being generated in.
         avgTpeak.append(max(s.datalist[SCSlice.ITEMP]))
         brho_i = np.where(s.datalist[SCSlice.ITEMP] == avgTpeak[len(avgTpeak)-1])
         avgBaseRho.append(s.datalist[SCSlice.IRHO][brho_i]/1.e5)
         t8 = avgTpeak[-1:][0] / 1.e8
         rctemp = (1.68e-4*np.exp(20.0/t8))**(1.0/2.3)
         rhoCrit.append(rctemp*1.e6/1.e5)

         ###TEST: Calculate global convective timescale
         #Calculate cutoff radii
         sp_st_den = sp_cen_den*sp_st_fac
         r_anelastic = None
         r_sp_start = None
         for r, rho in zip(s.datalist[SCSlice.IRADIUS], s.datalist[SCSlice.IRHO]):
            #TODO add error checking
            if not r_anelastic and rho <= an_cut:
               r_anelastic = r
               if r_sp_start:
                  break
            if not r_sp_start and rho <= sp_st_den:
               r_sp_start = r
               if r_anelastic:
                  break
         cbot = s.derived_datalist[SCSlice.IRCONV]
         ctop = cbot + s.derived_datalist[SCSlice.ILCONV]
         magvel = s.datalist[SCSlice.IMV]
         mv_rad = s.datalist[SCSlice.IRADIUS]

         #Change to dimensionless variables, only care about the convective zone
         li = np.where(mv_rad == cbot)[0]
         ri = np.where(mv_rad == ctop)[0]
         r_norm = mv_rad[li:ri]/ctop
         magvel_norm = magvel[li:ri] / magvel.max()
         mvn_inv = 1.0 / magvel_norm

         #Calculate global convective timescale as integral of 1/v over the convective zone
         #Convert back to physical units
         tcg = (ctop/magvel.max())*spint.trapz(mvn_inv, r_norm)
         tconvb.append(tcg)

         ###END TEST

         ratio.append(tnuc_x[len(tnuc_x)-1]/tconvb[len(tconvb)-1])

      #Fix any bad sorting of time array
      t_Ts, Tpeak = zip(*sorted(zip(t_T, Tpeak)))
      t_T, Tpeak_r = zip(*sorted(zip(t_T, Tpeak_r)))


      #NOTE: fontsize can be [size in points | 'xx-small' | 'x-small' |
      #  'small' | 'medium' | 'large' | 'x-large' | 'xx-large'

      #Build plots
      fig, ax_list = plt.subplots(nrows=2, ncols=1)

      #Temp
      ax_list[0].plot(t_T, Tpeak, color='red')
      ax_list[0].plot(t_slice, avgTpeak, color='red', marker='x', linestyle='None', label='avg peak temp')
      ax_list[0].set_ylabel(r'T$_{\mathrm{peak}}$ [$\times 10^8$ K]', color='red', fontsize='xx-large')
      ax_list[0].set_title('11030 Temperature and Radii', fontsize='xx-large')
      ax_list[0].set_ylim(1.85e8, 3.5e8)
      ax_list[0].tick_params(labelsize='xx-large')
      ax_list[0].yaxis.offsetText.set_visible(False)
      #ax_list[0].set_xlim(0, 40)
      tw = ax_list[0].twinx()
      tw.tick_params(labelsize='xx-large')
      tw.yaxis.offsetText.set_visible(False)
      tw.yaxis.set_major_locator(MaxNLocator(prune='lower'))
      #tw.set_xlim(0, 40)
      #tw.set_ylim(4.26e8, 4.38e8)
      tw.plot(t_T, Tpeak_r, color='green')
      tw.plot(t_slice, iface, color='black', label='CO/He Int.')
      tw.plot(t_slice, H, color='cyan', label='H')
      tw.plot(t_slice, rbot, color='cyan', label=r'$r_\mathrm{conv}$')
      #tw.plot(t_slice, rtop, color='cyan', marker='v', linestyle='None', label='rtop')
      tw.set_ylabel(r'T$_{\mathrm{peak}}$ radius [$\times 10^8$ cm]', color='green', fontsize='xx-large')
      #handles, labels = ax_list[0].get_legend_handles_labels()
      #fig.legend(handles, labels)
      tw.legend(loc=2, fontsize='x-large')

      #Plot base density (density at radius of peak temp)
      ax_list[1].plot(t_slice, avgBaseRho, label=r'$\rho_\mathrm{base}$')
      ax_list[1].plot(t_slice, rhoCrit, 'b--', 
            #label=r'$\rho_{WK} = \left( \frac{0.0607}{4 T_8^2} exp(20/T_8) \right)^{\frac{1}{2.3}} $')
            label=r'$\rho_{WK}$')

      #Customize axis labels (TODO: figure out how this works at some point)
      #label_text = [r'%i' % int(n/10**4) for n in plt.yticks()[0]]
      #ax_list[1].set_yticklabels(label_text)

      ax_list[1].tick_params(labelsize='xx-large')
      #ax_list[4].set_yscale('log')
      #ax_list[4].set_ylim(9.e5, 1.5e6)
      ax_list[1].set_ylabel(r'$\rho$ [$\times 10^5$ g cm$^{-3}$]', fontsize='xx-large')
      ax_list[1].set_xlabel(r'time [s]', fontsize='xx-large')
      ax_list[1].set_title('Density', fontsize='xx-large')
      ax_list[1].legend(loc=2, fontsize='x-large')

      #Set plot properties 
      fig.set_size_inches(10.0, 10.0)
      fig.tight_layout()
      fig.savefig("TRD.png", bbox_inches='tight')

   def plotConvTimescales(self, an_cut, sp_cen_den, sp_st_fac):
      """Plot a figure demonstrating ignition over time for publication, posters, presentation, etc..."""
      #TODO: I'm hacking this for a presentation, probably going to be ugly and need reworking
      import numpy as np
      import matplotlib.pyplot as plt
      import matplotlib
      from matplotlib.ticker import MaxNLocator
      import scipy.integrate as spint

      #Convenience aliases for diagnostic data
      t_T, Tpeak, Tpeak_r, Tpeak_vr  = self._diags._time, self._diags._maxT, self._diags._maxT_r, self._diags._maxT_vr
      t_M, Mpeak = self._diags._time, self._diags._maxmachfull
      t_e, epeak = self._diags._time, self._diags._maxenuc

      #Prepare variables for use in slice loop
      H = []
      iface = []
      rbot = []
      rtop = []
      t_slice = []

      tconv = []
      tconvb = []
      tnuc_x = []
      tnuc_xb = []
      tnuc_wk = []
      ratio = []

      avgTpeak = []
      avgBaseRho = []
      rhoCrit = []

      #Loop over slices in chronological order
      for s in sorted(self._slices, key=lambda sl: sl.timestamp[0]):
         #Build radius data from slices
         rbot.append(s.derived_datalist[SCSlice.IRCONV])
         rtop.append(s.derived_datalist[SCSlice.IRCONV] + s.derived_datalist[SCSlice.ILCONV])
         H.append(s.derived_datalist[SCSlice.IRCONV] + s.derived_datalist[SCSlice.IPSCALE])
         iface.append(s.derived_datalist[SCSlice.IINTFACE])
         t_slice.append(s.timestamp[1])

         #Estimate of convective turnover timescale and minimum nuclear timescale
         #tconv.append(s.derived_datalist[SCSlice.IPSCALE] / s.derived_datalist[SCSlice.IUCONV])
         tconv.append(s.derived_datalist[SCSlice.ILCONV] / s.derived_datalist[SCSlice.IUCONV])
         #tconv.append(s.derived_datalist[SCSlice.IUCONV])
         tnuc_x.append(min(s.derived_datalist[SCSlice.ITNUC][0]))
         tnuc_wk.append(min(s.derived_datalist[SCSlice.ITNUC][1]))
         tnuc_xb.append(min(s.derived_datalist[SCSlice.ITNUC][2]))
         #ratio.append(tnuc_x[len(tconv)-1]/tconv[len(tconv)-1])
         #ratio.append(tnuc_wk[len(tconv)-1]/tconv[len(tconv)-1])

         #Get the peak radially averaged temperature as an estimate of the background
         #conditions the hottest spot is being generated in.
         avgTpeak.append(max(s.datalist[SCSlice.ITEMP]))
         brho_i = np.where(s.datalist[SCSlice.ITEMP] == avgTpeak[len(avgTpeak)-1])
         avgBaseRho.append(s.datalist[SCSlice.IRHO][brho_i])
         t8 = avgTpeak[-1:][0] / 1.e8
         rctemp = (1.68e-4*np.exp(20.0/t8))**(1.0/2.3)
         rhoCrit.append(rctemp*1.e6)

         ###TEST: Calculate global convective timescale
         #Calculate cutoff radii
         sp_st_den = sp_cen_den*sp_st_fac
         r_anelastic = None
         r_sp_start = None
         for r, rho in zip(s.datalist[SCSlice.IRADIUS], s.datalist[SCSlice.IRHO]):
            #TODO add error checking
            if not r_anelastic and rho <= an_cut:
               r_anelastic = r
               if r_sp_start:
                  break
            if not r_sp_start and rho <= sp_st_den:
               r_sp_start = r
               if r_anelastic:
                  break
         cbot = s.derived_datalist[SCSlice.IRCONV]
         ctop = cbot + s.derived_datalist[SCSlice.ILCONV]
         magvel = s.datalist[SCSlice.IMV]
         mv_rad = s.datalist[SCSlice.IRADIUS]

         #Change to dimensionless variables, only care about the convective zone
         li = np.where(mv_rad == cbot)[0]
         ri = np.where(mv_rad == ctop)[0]
         r_norm = mv_rad[li:ri]/ctop
         magvel_norm = magvel[li:ri] / magvel.max()
         mvn_inv = 1.0 / magvel_norm

         #Calculate global convective timescale as integral of 1/v over the convective zone
         #Convert back to physical units
         tcg = (ctop/magvel.max())*spint.trapz(mvn_inv, r_norm)
         tconvb.append(tcg)

         ###END TEST

         ratio.append(tnuc_x[len(tnuc_x)-1]/tconvb[len(tconvb)-1])

      #Fix any bad sorting of time array
      t_Ts, Tpeak = zip(*sorted(zip(t_T, Tpeak)))
      t_T, Tpeak_r = zip(*sorted(zip(t_T, Tpeak_r)))


      #NOTE: fontsize can be [size in points | 'xx-small' | 'x-small' |
      #  'small' | 'medium' | 'large' | 'x-large' | 'xx-large'

      #Build plots
      fig, ax_list = plt.subplots(nrows=2, ncols=1)

      #Temp
      ax_list[0].plot(t_T, Tpeak, color='red')
      ax_list[0].plot(t_slice, avgTpeak, color='red', marker='x', linestyle='None', label='avg peak temp')
      ax_list[0].set_ylabel(r'T$_{\mathrm{peak}}$ [$\times 10^8$ K]', color='red', fontsize='xx-large')
      ax_list[0].set_title('Temperature and Radii', fontsize='xx-large')
      ax_list[0].set_ylim(1.85e8, 1.95e8)
      ax_list[0].tick_params(labelsize='xx-large')
      ax_list[0].yaxis.offsetText.set_visible(False)
      #ax_list[0].set_xlim(0, 40)
      tw = ax_list[0].twinx()
      tw.tick_params(labelsize='xx-large')
      tw.yaxis.offsetText.set_visible(False)
      tw.yaxis.set_major_locator(MaxNLocator(prune='lower'))
      #tw.set_xlim(0, 40)
      #tw.set_ylim(4.26e8, 4.38e8)
      tw.plot(t_T, Tpeak_r, color='green')
      tw.plot(t_slice, iface, color='black', label='CO/He Int.')
      tw.plot(t_slice, H, color='cyan', label='H')
      tw.plot(t_slice, rbot, color='cyan', label=r'$r_\mathrm{conv}$')
      #tw.plot(t_slice, rtop, color='cyan', marker='v', linestyle='None', label='rtop')
      tw.set_ylabel(r'T$_{\mathrm{peak}}$ radius [$\times 10^8$ cm]', color='green', fontsize='xx-large')
      #handles, labels = ax_list[0].get_legend_handles_labels()
      #fig.legend(handles, labels)
      tw.legend(loc=2, fontsize='x-large')

      #Plot base density (density at radius of peak temp)
      ax_list[1].plot(t_slice, avgBaseRho, label=r'$\rho_\mathrm{base}$')
      ax_list[1].plot(t_slice, rhoCrit, 'b--', label=r'$\rho_{WK}$')

      #Customize axis labels (TODO: figure out how this works at some point)
      #label_text = [r'%i' % int(n/10**4) for n in plt.yticks()[0]]
      #ax_list[1].set_yticklabels(label_text)

      ax_list[1].tick_params(labelsize='xx-large')
      #ax_list[4].set_yscale('log')
      #ax_list[4].set_ylim(9.e5, 1.5e6)
      ax_list[1].set_ylabel(r'$\rho$ [g cm$^{-3}$]', fontsize='xx-large')
      ax_list[1].set_xlabel(r'time [s]', fontsize='xx-large')
      ax_list[1].set_title('Base Density', fontsize='xx-large')

      #Set plot properties 
      fig.set_size_inches(10.0, 10.0)
      fig.tight_layout()
      fig.savefig("TRD.png", bbox_inches='tight')


   def plotTimescales(self, step):
      """Plot averaged timescales as a function of radius for the given timestep."""
      for s in self._slices:
         if s.timestamp[0] == step:
            s.plotTimescales()

   def plotHotspots(self, step, rlim=None, templog=False, reset=False):
      """Plot radius and temperature histograms for the given timestep's top hotspots 
      as well as temperature contours."""
      #First I need the interesting radii details.  
      #TODO: This is a stupid inefficient way to get them.  Need to rewrite/restructure.
      radii = (None, None, None) 
      for s in self._slices:
         if s.timestamp[0] == step:
            radii = (s.derived_datalist[SCSlice.IRCONV],
                     s.derived_datalist[SCSlice.IPSCALE],
                     s.derived_datalist[SCSlice.IINTFACE])
      #Find the right hotspot, call its plotting function 
      for hs in self._hotspots:
         if hs.timestamp[0] == step:
            hs.plotHotspots(radii, rlim=rlim, templog=templog, reset=reset)

   def plotHotspotsProj(self, step):
      """Plot a projection of hotspots onto a sphere for all of this timestep's hotspots."""
      for hs in self._hotspots:
         if hs.timestamp[0] == step:
            hs.plotHotspotsProj()

   def plotEP(self, step):
      """Plot entropy and pressure vs radius for the given timestep."""
      for s in self._slices:
         if s.timestamp[0] == step:
            s.plotEP()

   def plotRVSlice(self, step, poll, anelastic_cutoff, sp_cen_den, sp_st_fac):
      """Display a slice of a 3D pseudocolor plot of radial velocity generated by VisIt
      for this timestep. If the plot image isn't available, then generate it with VisIt on lens."""
      for s in self._slices:
         if s.timestamp[0] == step:
            s.plotRVSlice(poll, anelastic_cutoff, sp_cen_den, sp_st_fac)

   def genRVSlices(self, anelastic_cutoff, sp_cen_den, sp_st_fac):
      """Generate RV Slices for all available pltfiles. Each slice can be displayed with 
      plotRVSlice()."""
      from subprocess import PIPE, STDOUT, Popen

      #Make the input files needed by the VisIt script
      self._genRVSliceInfiles(anelastic_cutoff, sp_cen_den, sp_st_fac)

      #Spawn subprocess to execute the VisIt script
      lenssub = Popen(['lensSubmit/lenssub.sh', '/ccs/home/ajacobs/Projects/SubChandra/lensSubmit/genRVPlots.py'], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
      #The old version of the IPython notebook available on Titan has no ability to
      #accept stdin, so remind the user they must enter their password for lens in the
      #terminal running the notebook.
      print('Enter your lens password into the terminal running this notebook.')

   def gen3DRVPlots(self, min_rv, xmax, anelastic_cutoff, sp_cen_den, sp_st_fac):
      """Generate 3D RV plots for all available pltfiles. Each plots can be displayed with 
      plotRV3D()."""
      from subprocess import PIPE, STDOUT, Popen

      #Make the input files needed by the VisIt script
      self._gen3DRVInfiles(min_rv, xmax, anelastic_cutoff, sp_cen_den, sp_st_fac)

      #Spawn subprocess to execute the VisIt script
      lenssub = Popen(['lensSubmit/lenssub.sh', '/ccs/home/ajacobs/Projects/SubChandra/lensSubmit/gen3DRVPlots.py'], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
      #The old version of the IPython notebook available on Titan has no ability to
      #accept stdin, so remind the user they must enter their password for lens in the
      #terminal running the notebook.
      print('Enter your lens password into the terminal running this notebook.')

   def plotRV3D(self, step, min_rv, poll, anelastic_cutoff, sp_cen_den, sp_st_fac):
      """Display a 3D volume rendering plot of positive radial velocity generated by VisIt
      for this timestep. If the plot image isn't available, then generate it with VisIt on lens."""
      for s in self._slices:
         if s.timestamp[0] == step:
            s.plotRV3D(min_rv, poll, anelastic_cutoff, sp_cen_den, sp_st_fac)

   ##Private methods##
   def _getPltFiles(self):
      """Get a list of the full path to all valid pltfile directories for this simulation."""
      from os import listdir
      from os.path import join, basename, isfile
      from re import match

      #All plotfiles assumed to be form of <some label>_plt##### 
      #where the name can end in 5 or 6 numbers
      PLTFILERE = '.+_plt[0-9]{5,6}$' 
      PLOTDIR = 'plotfiles'

      pltdir = join(self._scratch_dir, PLOTDIR)
      candidates = [join(base, c) for c in listdir(base)] #Use a list comprehension to add the full path to each entry
      candidates += [join(pltdir, c) for c in listdir(pltdir)]
      pltfiles = []
      for c in candidates:
         cc = basename(c)
         if match(PLTFILERE, cc): #Match the regex
            if (cc not in [basename(f) for f in pltfiles]): #Check for duplicate
               #We make sure this isn't an empty or incomplete pltfile 
               #by making sure the last file written, Header, is present
               if isfile(join(c, 'Header')):
                  pltfiles.append(c)

      return pltfiles

   def _genRVSliceInfiles(self, anelastic_cutoff, sp_cen_den, sp_st_fac):
      """Generate the infiles needed to run the script that plots multiple RV slices."""
      from os.path import join, basename
      from numpy import sqrt

      LENSSUB_DIR = '/ccs/home/ajacobs/Projects/SubChandra/lensSubmit'
      PLOTRVS_IN = join(LENSSUB_DIR, 'genRVPlots.in')
      PLTS_DIR = join(self._stage_dir, 'plots')

      #Create the top-level infile
      with open(PLOTRVS_IN, 'w') as prvs_in:
         #Get a reverse-sorted list of all valid pltfiles
         pltfiles = self._getPltFiles()
         pltfiles = sorted(pltfiles, key=basename)
         pltfiles.reverse()

         #All plts have same xmax, so grab one and use it to calculate xmax
         rmax = self._slices[0].datalist[SCSlice.IRADIUS].max()
         xmax = rmax / sqrt(3.0)

         prvs_in.write('#Maximum x (same as max y, z)\n')
         prvs_in.write('xmax = {0}\n\n'.format(xmax))

         #These 2D RV slices make use of information from the radial slices, so
         #we only attempt to plot 2D slices for pltfiles with a corresponding
         #1D slice.  For each such pair, generate an input file used by the VisIt
         #plotting script
         i = 0
         sub_infiles = []
         curplt = pltfiles.pop()
         for s in sorted(self._slices, key=lambda sl: sl.timestamp[0]):
            if s.my_pltfile == basename(curplt):
               #Build and store sub-infile name
               infilename = join(LENSSUB_DIR, 'genRVPlots-{0}.in'.format(i))
               sub_infiles.append(infilename)
               i += 1

               #Store important radii
               cbot = s.derived_datalist[SCSlice.IRCONV]
               ctop = cbot + s.derived_datalist[SCSlice.ILCONV]
               H = s.derived_datalist[SCSlice.IPSCALE]

               #Calculate cutoff radii
               sp_st_den = sp_cen_den*sp_st_fac
               r_anelastic = None
               r_sp_start = None
               for r, rho in zip(s.datalist[SCSlice.IRADIUS], s.datalist[SCSlice.IRHO]):
                  #TODO add error checking
                  if not r_anelastic and rho <= anelastic_cutoff:
                     r_anelastic = r
                     if r_sp_start:
                        break
                  if not r_sp_start and rho <= sp_st_den:
                     r_sp_start = r
                     if r_anelastic:
                        break

               #Build savefilename
               base = self._label + '_RV' + str(s.timestamp[0]) + '.png'
               savefile = join(PLTS_DIR, base)

               #TODO: Only write infile if png doesn't exist?
               #Write the sub-infile
               pltfull = join(curplt, 'Header')
               s._writeRVIn(infilename, pltfull, cbot, ctop, H, rmax, r_anelastic, r_sp_start, savefile)

               #Advance to the next pltfile
               curplt = pltfiles.pop()

         prvs_in.write('#Input files for each individual frame\n')
         prvs_in.write('plt_count = {0}\n'.format(i))
         for si in sub_infiles:
            prvs_in.write('{0}\n'.format(si))

   def _gen3DRVInfiles(self, min_rv, xmax, anelastic_cutoff, sp_cen_den, sp_st_fac):
      """Generate the infiles needed to run the script that plots multiple 3D RV plots."""
      from os.path import join, basename
      from numpy import sqrt

      LENSSUB_DIR = '/ccs/home/ajacobs/Projects/SubChandra/lensSubmit'
      PLOT3DRVS_IN = join(LENSSUB_DIR, 'gen3DRVPlots.in')
      PLTS_DIR = join(self._stage_dir, 'plots')

      #Create the top-level infile
      with open(PLOT3DRVS_IN, 'w') as prvs_in:
         #Get a list of all valid pltfiles
         pltfiles = self._getPltFiles()

         #Write xmax and min_rv
         prvs_in.write('#Maximum x (same as max y, z)\n')
         prvs_in.write('xmax = {0}\n\n'.format(xmax))
         prvs_in.write('#Minimum radial velocity to plot\n')
         prvs_in.write('min_rv = {0}\n\n'.format(min_rv))

         #Write sub-infiles
         sub_infiles = []
         i = 0
         for pf in pltfiles:
            #Build and store sub-infile name
            infilename = join(LENSSUB_DIR, 'gen3DRVPlots-{0}.in'.format(i))
            sub_infiles.append(infilename)
            i += 1

            #Store important radii
            #cbot = s.derived_datalist[SCSlice.IRCONV]
            #ctop = cbot + s.derived_datalist[SCSlice.ILCONV]
            #H = s.derived_datalist[SCSlice.IPSCALE]

            #Calculate cutoff radii
            #sp_st_den = sp_cen_den*sp_st_fac
            #r_anelastic = None
            #r_sp_start = None
            #for r, rho in zip(s.datalist[SCSlice.IRADIUS], s.datalist[SCSlice.IRHO]):
            #  #TODO add error checking
            #  if not r_anelastic and rho <= anelastic_cutoff:
            #    r_anelastic = r
            #    if r_sp_start:
            #      break
            #  if not r_sp_start and rho <= sp_st_den:
            #    r_sp_start = r
            #    if r_anelastic:
            #      break

            #Build savefilename
            bpf = basename(pf)
            ts_idx = bpf.find('_plt') + 4
            timestep = bpf[ts_idx:]
            base = self._label + '_3DRV' + timestep + '.png'
            savefile = join(PLTS_DIR, base)

            #TODO: Only write infile if png doesn't exist?
            #Write the sub-infile
            pltfull = join(pf, 'Header')
            self._write3DRVIn(infilename, pltfull, xmax, min_rv, savefile)

         prvs_in.write('#Input files for each individual frame\n')
         prvs_in.write('plt_count = {0}\n'.format(i))
         for si in sub_infiles:
            prvs_in.write('{0}\n'.format(si))

   def _write3DRVIn(self, infilename, pltfile, xmax, min_rv, savefile):
      """Write infile for a single frame of 3D RV plot."""
      from numpy import sqrt

      with open(infilename, 'w') as prv_in:
         #Write header
         prv_in.write('#Inputs for the plotRV3D VisIt script\n\n')

         #Write pltfile's location
         prv_in.write('#Database/file to be read\n')
         prv_in.write('pltfile = {0}\n\n'.format(pltfile))

         #Maximum x
         prv_in.write('#Maximum x (same as max y, z)\n')
         prv_in.write('xmax = {0}\n\n'.format(xmax))

         #Minimum radial velocity to plot
         prv_in.write('#Minimum radial velocity to plot\n')
         prv_in.write('min_rv = {0}\n\n'.format(min_rv))

         #File to save to
         prv_in.write('#File to save the plot to\n')
         prv_in.write('savefile = {0}\n'.format(savefile))

class SCTempHist(object):
   """Temperature histograms for all refinement lenghtscales in a pltfile."""

   #Constructor
   def __init__(self, thistfile, parent):
      """self      --> implicitly passed reference to this instance of SCSlice
         slicefile --> a .temphist file generated by fsubchandra.f90
         parent    --> reference to this histogram's parent SCOutput"""
      from os.path import basename

      (self.ells, self.temps, self.counts, self.timestamp) = self._parseTHfile(thistfile)
      self.sco_parent = parent
      self.TMax = self._calcMax(self.temps, self.counts[1])

   ## Public Methods ##

   ## Private Methods ##
   def _calcMax(self, bin_arr, counts_arr):
      """Calculate the largest non-zero bin for histogram data in bin and counts arrays."""
      import numpy as np

      vmax_arr = []
      for b, c in zip(bin_arr, counts_arr):
         vmax = -1.0
         for val, num in zip(b,c):
            if num > 0 and val > vmax:
               vmax = val
         vmax_arr.append(vmax)

      return np.array(vmax_arr)

   def _parseTHfile(self, thistfile):
      """Return (ells, temp array, counts array, and timestep) based on thistfile."""
      import numpy as np
  
      thf = open(thistfile)
      tret = (None, None) #(time step, time in seconds)
  
      #Parse header / global info
      #TODO: I make strict assumptions about the structure of the header,
      #might be nice to develop a more portable/generic algorithm
      for line in thf:
         if not line.strip():
            continue #If blank line, just loop to next line
         if not line.strip().startswith('#') :
            break #Reached end of header
  
         if line.count('time') == 1:
            tokens = line.split()
            time = float( tokens[len(tokens)-1] )
            break
  
  
      #Record time info
      li = thistfile.index('_plt') + 4
      ri = thistfile.index('.temphist') 
      step = int(thistfile[li:ri])
      tret = (step, time)
 
      #For each lengthscale, collect the counts for each temperature bin
      #The format of the file after the header is
      #    <lengthscale of finest level>
      #    <temp bin 1>  <counts>
      #    ...
      #    <lengthscale of next level down>
      #    <temp bin 1>  <counts>
      #    ...
      #    ...
      #    <lengthscale of coarsest level>
      #    <temp bin 1>  <counts>
      #    ...
      ells = []
      temps = []   #temps[i][j]  --> temperature bin j for lengthscale i
      counts = []  #counts[i][j] --> counts for temperature in bin j for lengthscale i
      counts_eos = []
      for line in thf:
         if line.strip().startswith('#'): #Skip comments
            continue
         tokens = line.split()
         if len(tokens) == 1:
            #Only one entry means it's a lengthscale
            ells.append(float(tokens[0]))
            temps.append([])
            counts.append([])
            counts_eos.append([])
         elif len(tokens) == 3:
            #Two entires means a temperature bin and its counts,
            #store them in the corresponding lengthscale's array
            temps[len(temps)-1].append(float(tokens[0]))
            counts[len(counts)-1].append(int(tokens[1]))
            counts_eos[len(counts_eos)-1].append(int(tokens[2]))
         else:
            print('Unexpected line format in file!')
            thf.close()
            sys.exit(2)
            
      thf.close()

      ells = np.array(ells)
      temps = np.array(temps)
      countsret = [np.array(counts), np.array(counts_eos)]

      return (ells, temps, countsret, tret)
 


class SCHotspots(object):
   """The hotspot data from a particular sub-Chandra simulation plotfile/timestamp."""
   ##Shared class data##
   TEMP_COL = 0
   RHO_COL = 1
   X_COL, Y_COL, Z_COL = 2, 3, 4
   LEV_COL  = 5

   ##Constructor##
   def __init__(self, hsfile, in_dict=None):
      """SCHotspots constructor.
 
      self    --> reference to this instance of SCHotspots (implicitly passed)
      hsfile  --> .hotspots file
      in_dict --> the inputs dictionary for this simulation"""
      self._hsfile = hsfile
      self.timestamp = self._parseHSFile()
      self._hsarray = None
      self.in_dict = in_dict

   ##Public methods##
   def plotHotspots(self, radii, rlim=None, templog=False, plot_top=False, reset=False):
      """Plot radius and temperature histograms for this timestep's top hotspots 
      as well as temperature contours."""
      import numpy as np
      import matplotlib.pyplot as plt
      import matplotlib
      from matplotlib import cm, colors
      from mpl_toolkits_ext.basemap import Basemap#, cm

      TCRIT = 2.25e7
      TMAX = 2.4e9

      #When tweaking plots it's nice if I can store the data, but then when I load new
      #data I need to reset.  Here I delete all data arrays so they'll be rebuilt.
      if(reset and hasattr(self, 'r')):
         del self.r
         del self.theta
         del self.phi
         del self.temp
         del self.logtemp
         del self.rho6
         del self.temp_lons
         del self.temp_lats

      #Get data, and save data so we only go through all the arrays once
      if not self._hsarray: #If no array yet, build it
         self._buildHSArr()
      if not hasattr(self, 'r'):
         self.r         = np.array([hs.loc[1][0] for hs in self._hsarray if hs.temp > TCRIT and hs.temp < TMAX])
      if not hasattr(self, 'theta'):
         self.theta     = np.array([hs.loc[1][1] for hs in self._hsarray if hs.temp > TCRIT and hs.temp < TMAX])
      if not hasattr(self, 'phi'):
         self.phi       = np.array([hs.loc[1][2] for hs in self._hsarray if hs.temp > TCRIT and hs.temp < TMAX])
      if not hasattr(self, 'temp'):
         self.temp      = np.array([hs.temp for hs in self._hsarray if hs.temp > TCRIT and hs.temp < TMAX])
      if not hasattr(self, 'logtemp'):
         self.logtemp   = np.log10(self.temp)
      if not hasattr(self, 'rho6'):
         self.rho6      = np.array([hs.rho/1.e6  for hs in self._hsarray if hs.temp > TCRIT and hs.temp < TMAX])
      if not hasattr(self, 'temp_lons'):
         self.temp_lons = np.array([deg for deg in np.degrees(self.phi)])
      if not hasattr(self, 'temp_lats'):
         self.temp_lats = np.array([-(deg-90.) for deg in np.degrees(self.theta)])

      #Local aliases for data
      r         = self.r
      #theta    = self.theta
      #phi      = self.phi
      temp      = self.temp
      logtemp   = self.logtemp
      rho6      = self.rho6
      temp_lons = self.temp_lons
      temp_lats = self.temp_lats

      #Critical hotspots
      #ctemp      = np.array([hs.temp for hs in self._hsarray if hs.rho/1.e6 > (1.68e-4*np.exp(20.0/(hs.temp/1.e8)))**(1.0/2.3)])
      #ctheta     = np.array([hs.loc[1][1] for hs in self._hsarray if hs.rho/1.e6 > (1.68e-4*np.exp(20.0/(hs.temp/1.e8)))**(1.0/2.3)])
      #cphi       = np.array([hs.loc[1][2] for hs in self._hsarray if hs.rho/1.e6 > (1.68e-4*np.exp(20.0/(hs.temp/1.e8)))**(1.0/2.3)])
      #ctemp_lons = np.array([deg for deg in np.degrees(cphi)])
      #ctemp_lats = np.array([-(deg-90.) for deg in np.degrees(ctheta)])

      crit_temp = 2.3e8
      #ctemp      = np.array([hs.temp for hs in self._hsarray if hs.temp > crit_temp])
      #ctheta     = np.array([hs.loc[1][1] for hs in self._hsarray if hs.temp > crit_temp])
      #cphi       = np.array([hs.loc[1][2] for hs in self._hsarray if hs.temp > crit_temp])
      #ctemp_lons = np.array([deg for deg in np.degrees(cphi)])
      #ctemp_lats = np.array([-(deg-90.) for deg in np.degrees(ctheta)])

      #Get important radii of interest
      rbot = radii[0]
      H = radii[1]
      iface = radii[2]

      #Get min, max temp
      min_temp = temp.min()
      max_temp = temp.max()
      hs_count = len(temp)

      #Calculate temperature bins for color map
      shsarr = sorted(self._hsarray, key=lambda hs: hs.temp)
      stemp = sorted(temp)
      tlevs = [temp.min()]
      tllevs = [logtemp.min()]
      count = 0
      for t in stemp:
         count += 1
         if count > (1./9.)*hs_count:
            tlevs.append(t)
            tllevs.append(np.log10(t))
            count = 0
         #For hottest temps break into top 9% and top 1%
         #if len(tlevs) == 9:
         #  if count > 0.09*hs_count:
         #    tlevs.append(hs.temp)
         #    tllevs.append(np.log10(hs.temp))
         #    count = 0
         #else:
         #  if count > 0.1*hs_count:
         #    tlevs.append(hs.temp)
         #    tllevs.append(np.log10(hs.temp))
         #    count = 0
      else:
         tlevs.append(t)
         tllevs.append(np.log10(t))

      #Build colormap
      #TODO: Get nice paper-worthy plot for 'ignition' section
      #TODO: Figure out discrepancy between pltfile cells and valid cells
      #      Update: currently working with OLCF on this, seems to only happen with optimized Cray compiler
      #rr = np.linspace(1.0, 0.7, 11)
      #gb = np.linspace(1.0, 0.0, 11)
      #temp_cmap = colors.ListedColormap([
      #  (rr[0],  gb[0],  gb[0]),   #Coldest
      #  (rr[1],  gb[1],  gb[1]),
      #  (rr[2],  gb[2],  gb[2]),
      #  (rr[3],  gb[3],  gb[3]),
      #  (rr[4],  gb[4],  gb[4]),
      #  (rr[5],  gb[5],  gb[5]),
      #  (rr[6],  gb[6],  gb[6]),
      #  (rr[7],  gb[7],  gb[7]),
      #  (rr[8],  gb[8],  gb[8]),
      #  (rr[9],  gb[9],  gb[9]),
      #  (rr[10], gb[10], gb[10])]) #Hottest
      #  #(1,      1,      0)]) #Hottest
      # RGB colors from 9-class OrRd at colorbrewer2.org
      temp_cmap = colors.ListedColormap([
         (255./255., 247./255., 236./255.),   #Coldest
         (254./255., 232./255., 200./255.),
         (253./255., 212./255., 158./255.),
         (253./255., 187./255., 132./255.),
         (252./255., 141./255.,  89./255.),
         (239./255., 101./255.,  72./255.),
         (215./255.,  48./255.,  31./255.),
         (179./255.,   0./255.,   0./255.),
         (127./255.,   0./255.,   0./255.)]) #Hottest
      tc_bounds = tlevs
      tc_norm = colors.BoundaryNorm(tc_bounds, temp_cmap.N)

      #Calculate critical density range based on eqn (8) of Woosley & Kasen 2011
      t8 = max_temp/1.e8
      rho_critl = (1.68e-4*np.exp(20.0/t8))**(1.0/2.3)
      t8 = min_temp/1.e8
      rho_critr = (1.68e-4*np.exp(20.0/t8))**(1.0/2.3)

      print('Min, Max temp of {0} hottest cells: {1}, {2}'.format(hs_count, min_temp, max_temp))
      print('Critical density range: [{0}, {1}]'.format(rho_critl, rho_critr))
      sys.stdout.flush()

      if 'inline' in matplotlib.get_backend():
         #Build plots
         fig = plt.figure()
         #subplot2grid call signature: (grid_row, grid_cols), (subplot_row, subplot_col), colspan=1, rowspan=1 
         ax_rad  = plt.subplot2grid((3,3), (0,0))  
         ax_temp = plt.subplot2grid((3,3), (0,1))  
         ax_rho  = plt.subplot2grid((3,3), (0,2))  
         ax_proj = plt.subplot2grid((3,3), (1,0), colspan=2, rowspan=2)  

         #Plot temperature histogram
         ax_temp.hist(temp, bins=1000)
         ax_temp.set_xlabel("temperature (K)")

         #Build projection map for temperature theta, phi locations
         map = Basemap(projection='nsper', lon_0=45, lat_0=45, 
               llcrnrlon=-180, llcrnrlat=-90, urcrnrlon=180, urcrnrlat=90, resolution=None, ax=ax_proj)
         #map.drawmeridians(np.arange(0, 90, 15), color="0.65", latmax=90)
         #map.drawparallels(np.arange(0, 90, 15), color="0.65", latmax=90) #, labels=[1,0,0,1])

         #It's stupid that I have to do this, but below I "erase" extraneous latitude lines
         #by writing over them with thick white lines.
         #for lat in range(0,90,15):
         #  for long in range(15,180,15):
         #    #Erase extraneous latitude lines on the left
         #    left_long=-long
         #    map.drawgreatcircle(left_long+15, lat, left_long, lat, linewidth=5, color="w")
         #    #Erase extraneous latitude lines on the right
         #    right_long=long+90
         #    map.drawgreatcircle(right_long-15, lat, right_long, lat, linewidth=5, color="w")
         ##Same with extraneous longitude lines at the bottom
         #map.drawgreatcircle(0,    0,  0, -25, linewidth=5, color="w")
         #map.drawgreatcircle(15,   0, 15, -30, linewidth=5, color="w")
         #map.drawgreatcircle(30,   0, 30, -35, linewidth=5, color="w")
         #map.drawgreatcircle(45,   0, 45, -30, linewidth=5, color="w")
         #map.drawgreatcircle(60,   0, 60, -35, linewidth=5, color="w")
         #map.drawgreatcircle(75,   0, 75, -30, linewidth=5, color="w")

         # draw the boundary of our domain -- we want great circles here
         # note that we draw in 15 degree increments.  Otherwise the lat/long grid
         # doesn't line up with the boundary 
         #Left boundary
         for lat in range(0,90,15):
            map.drawgreatcircle(0, lat, 0, lat+15, linewidth=1, color="k")
         #Right boundary
         for lat in range(0,90,15):
            map.drawgreatcircle(90, lat, 90, lat+15, linewidth=1, color="k")
         #Bottom boundary
         for lon in range(0,90,15):
            map.drawgreatcircle(lon, 0, lon+15, 0, linewidth=1, color="k")

         if templog:
            clevs = np.linspace(logtemp.min(), logtemp.max(), 11) 
            #cs = map.contourf(temp_lons, temp_lats, logtemp, clevs, latlon=True, tri=True, cmap=cm.Reds)
            cs = map.contourf(temp_lons, temp_lats, logtemp, tllevs, latlon=True, tri=True, cmap=cm.jet)
         else:
            #clevs = np.linspace(temp.min(), temp.max(), 11) 
            clevs = np.linspace(2.25e8, crit_temp, 11) 
            cs = map.contourf(temp_lons, temp_lats, temp, tlevs, latlon=True, tri=True, cmap=temp_cmap, norm=tc_norm)
            #map.contourf(ctemp_lons, ctemp_lats, ctemp, clevs, latlon=True, tri=True, cmap=cm.Greens)
         cbar = map.colorbar(cs, location='right', pad='5%', ticks=tlevs)
         cbar.set_label('Kelvin')  #cbar.set_label('temperature ($\times 10^8$ Kelvin)')

         #Plot radius histogram with temperature color-coding
         # color-coding is achieved by plotting several bars instead of using
         # ax.hist(), and we use the projection map's colorbar
         #ax_rad.hist(r, bins=1000)

         dr = 1.e5 #use a dr of 1 km, which is roughly the radial resolution in these simulations
         radii, counts, cols = self._binData(r, temp, dr, cbar)

         for r, c, col in zip(radii, counts, cols):
            ax_rad.bar(r, c, width=dr, color=col, edgecolor=col, align='center') 
         #ax_rad.bar(radii, counts, width=dr, color=(1.0, 1.0, 1.0), align='center') 
         ax_rad.set_xlabel("radius (cm)")
         if rlim:
            ax_rad.set_xlim(rlim[0], rlim[1])
         # plot the radii (CO/He interface, start of convective region, top of convective region)
         ax_rad.plot([iface, iface], [0, ax_rad.get_ylim()[1]], color="k", linestyle='--')
         ax_rad.plot([rbot, rbot], [0, ax_rad.get_ylim()[1]], color="k", linestyle='--')
         if plot_top:
            ax_rad.plot([rbot + H, rbot + H], [0, ax_rad.get_ylim()[1]], color="k", linestyle='--')

         #Annotate the interface line
         ifrac = (iface - ax_rad.get_xlim()[0]) / (ax_rad.get_xlim()[1] - ax_rad.get_xlim()[0])
         ax_rad.annotate('WD/He\ninterface', 
               xy=(ifrac,0.8), 
               xytext=(-60,-30),
               xycoords='axes fraction', textcoords='offset points',
               arrowprops=dict(facecolor='black', arrowstyle='simple',
                  connectionstyle='arc3,rad=-0.2'))

         #Annotate the convective base line
         cbfrac = (rbot - ax_rad.get_xlim()[0]) / (ax_rad.get_xlim()[1] - ax_rad.get_xlim()[0])
         ax_rad.annotate('Convective\nbase', 
               xy=(cbfrac,0.8), 
               xytext=(30,-30),
               xycoords='axes fraction', textcoords='offset points',
               arrowprops=dict(facecolor='black', arrowstyle='simple',
                  connectionstyle='arc3,rad=-0.2'))

         #Plot density histogram with color-coding
         drho6 = 0.001
         rho6_bins, rho_counts, rho_colors = self._binData(rho6, temp, drho6, cbar)
         for r6, c, col in zip(rho6_bins, rho_counts, rho_colors):
            ax_rho.bar(r6, c, width=drho6, color=col, edgecolor=col, align='center') 
         #ax_rho.hist(rho6, bins=1000)
         ax_rho.set_xlabel(r"density ($\times 10^6$ g cm$^{-3}$)")
         #ax_rho.fill_between([rho_critl, rho_critr], 0, 500, facecolor='0.9', edgecolor='1.0')
         ax_rho.plot([rho_critr, rho_critr], [0, ax_rho.get_ylim()[1]], color="k", linestyle='--')

         #Annotate the critical density line
         rhofrac = (rho_critr - ax_rho.get_xlim()[0]) / (ax_rho.get_xlim()[1] - ax_rho.get_xlim()[0])
         ax_rho.annotate(r'$\rho_{\mathrm{cr},\mathrm{WK}}$', 
               xy=(rhofrac,0.8), size=12.5,
               xytext=(30,-30),
               xycoords='axes fraction', textcoords='offset points',
               arrowprops=dict(facecolor='black', arrowstyle='simple',
                  connectionstyle='arc3,rad=-0.2'))

         #Set plot properties 
         fig.set_size_inches(15.0, 12.5)
         #This fixes a problem with mpl's pstoeps converter when using ghostscript as distiller
         #matplotlib.rc('ps', usedistiller='xpdf')
         fig.savefig("test_ig.png", bbox_inches='tight')

      else:
         #TODO: Need to implement non-inline plotting
         pass

   ##Private methods##
   def _parseHSFile(self):
      """Parse the hotspots file, return (timestep, time)."""
      #Get timestep
      li = self._hsfile.index('_plt') + 4
      ri = self._hsfile.index('.hotspots') 
      step = int(self._hsfile[li:ri])

      #Get time
      f = open(self._hsfile)
      for line in f:
         if line.strip().startswith('# time'):
            tokens = line.partition('=')
            time = float(tokens[2])
            break
      else:
         assert False, "Invalid .hotspots file, couldn't find time in header."
      f.close()

      return (step, time)

   def _buildHSArr(self):
      """Build array of hotspots on demand."""
      print('Building hotspots array for step {0:5d}'.format(self.timestamp[0]))
      sys.stdout.flush()  

      #Build hotspots array
      self._hsarray = []
      f = open(self._hsfile)
      for i, line in enumerate(f):
         #Skip comments and blanks
         blank = not line.strip()
         if line.strip().startswith('#') or blank:
            continue
         tokens = line.split()
         temp = float(tokens[SCHotspots.TEMP_COL])
         rho = float(tokens[SCHotspots.RHO_COL])
         x, y, z = float(tokens[SCHotspots.X_COL]), float(tokens[SCHotspots.Y_COL]), float(tokens[SCHotspots.Z_COL])
         lev = float(tokens[SCHotspots.LEV_COL])
         self._hsarray.append(_Hotspot(temp, rho, x, y, z, lev, dx))
         if (i % 100000 == 0): 
            print('Read 100K!, Total: {0:7d}'.format(i))
            #break
         nx = float(self.in_dict['n_cellx'])
         xmax = float(self.in_dict['prob_hi_x'].replace('d', 'e'))
         dx = xmax/(nx*2**(lev-1))

      f.close()

   def _buildIgArr(self, critT):
      """Build array of ignitors on demand."""
      from datetime import datetime

      print('Building ignitors for step {0:5d}'.format(self.timestamp[0]))
      sys.stdout.flush()  

      #Pick out the igniting cells from the hotspots array 
      #(as determined by a critical temperature).
      igset = []
      for hs in self._hsarray:
         if hs.temp > critT:
            igset.append(hs)
            if len(igset) > 5000: 
               break

      print('igset size: ', len(igset), datetime.utcnow())
      sys.stdout.flush()  
      #Build Ignitor arrays out of igniting hotspots
      #  Each array is a collection of hotspots near one another
      #  We sort igset by temperature first so that searches should happen
      #  radially outward from the center of unique igniting volumes
      #igset_srt = sorted(igset, key=lambda ele: ele.temp, reverse=True)
      self._igarrays = []

      #Build the initial sets, which will allow for multiple arrays containing
      #the same hotspot
      print('build initial set ', datetime.utcnow())
      sys.stdout.flush()  
      initial_sets = []
      for hs in igset:
         stored = False
         sti = -1
         for i in range(len(initial_sets)):
            #If hs is near any of the hs's in arr, append it
            if self._nearby(initial_sets[i], hs):
               initial_sets[i].append(hs)
               stored = True
               sti = i
         if not stored:
            #If hs didn't find a home, made a new igarray
            initial_sets.append([hs,])

      #Go through the sets and flag for merging any that have a non-empty intersection.
      print('flag for merging', datetime.utcnow())
      sys.stdout.flush()  
      merge_map = dict([(i, i) for i in range(len(initial_sets))])   #Maps merged indices to the new index
      #for hs in igset_srt:
      for hs in igset:
         stored = False
         sti = -1
         for i in range(len(initial_sets)):
            #If hs is near any of the hs's in arr, append it
            if self._nearby(initial_sets[i], hs):
               if stored:
                  merge_map[i] = sti
               else:
                  stored = True
                  sti = i
      
      #Merge intersecting sets 
      print('merge', datetime.utcnow())
      sys.stdout.flush()  
      for k in merge_map:
         v = merge_map[k]
         if (k != v):
            initial_sets[v] = list(set(initial_sets[v] + initial_sets[k]))

      #Create the final list of ignition arrays
      #by only including merged lists
      print('build final list', datetime.utcnow())
      sys.stdout.flush()  
      self._igarrays = []
      for k in merge_map:
         v = merge_map[k]
         if (k == v):
            self._igarrays.append(initial_sets[k])
      
      #Make array of _Ignitor objects
      self._Ig_arr = [_Ignitor(arr) for arr in self._igarrays]
           
   def _nearby(self, arr, hs):
      """Return true if hs is within 10dx of any hs in arr, false otherwise."""
      from numpy import sqrt
      if arr is None:
         return False
      for h in arr:
         x1, y1, z1 = h.loc[0][0], h.loc[0][1], h.loc[0][2]
         x2, y2, z2 = hs.loc[0][0], hs.loc[0][1], hs.loc[0][2]
         dist = sqrt( (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2 )
         if dist < 30.0e5: #30 km
            return True
      return False
        

   def _binData(self, data, color_data, dx, cmap):
      """Bin data with bin size dx and return corresponding colors based on the max value in 
      color_data using the given colormap."""
      import numpy as np

      #Build data bins
      dmin = data.min()
      dmax = data.max()
      bin_cnt = int(round((dmax-dmin)/dx))
      data_bins = np.linspace(dmin + dx*0.5, dmax - dx*0.5, num=bin_cnt)
 
      #Prepare arrays and variables for loop
      counts = np.zeros(bin_cnt + 1, dtype = np.int32)
      cols = []
      sorted_data_col = sorted(zip(data, color_data))
      i = 0
      upper_bound = dmin + dx
      avg_col = 0.0
      max_col = 0.0

      #Loop through parallel lists of data and color_data, sorted by data.
      #Calculate the counts, the average color_data, and the maximum color_data for each bin. 
      #Convert color_data to color
      for dat, col_dat in sorted_data_col:
         if col_dat > max_col:
            max_col = col_dat
         if dat < upper_bound:
            counts[i] += 1
            avg_col += col_dat
         else:
            #Before incrementing, append avg_col's color
            avg_col = avg_col / counts[i]
            cols.append(cmap.to_rgba([avg_col,]))
            #cols.append(cmap.to_rgba([max_col,]))
            #cols.append((0.5,0.5,0.5,1))
 
            #Increment to next bin, don't forget to count the data and color just found
            #in this loop iteration
            upper_bound += dx
            i += 1
            counts[i] += 1
            avg_col = col_dat
            max_col = 0.0

      return (data_bins, counts, cols)

class _Hotspot(object):
   """Class representing a single hotspot"""
   ##Shared class data

   ##Constructor
   def __init__(self, temp, rho, x, y, z, lev, dx):
      """_Hotspot constructor.
 
      self    --> reference to this instance of _Hotspot (implicitly passed)
      temp    --> temperature of hotspot [K]
      x, y, z --> Cartesian location of hotspot [cm]
      lev     --> refinement level of hotspot's location
      dx      --> size of hotspot's cell"""

      self.temp = temp
      self.rho = rho
      #Loc = ( (x, y, z) , (r, theta, phi) )
      self.loc = ( (x, y, z), self._getSphTuple(x, y, z) )
      self.lev = lev
      self.dx = dx

   ##Public methods

   ##Private methods
   def _getSphTuple(self, x, y, z):
      """Return tuple of (r, theta, phi[azimuth]) based on Cartesian x,y,z."""
      import numpy as np

      r = np.sqrt(x**2 + y**2 + z**2)
      theta = np.arccos(z/r)
      phi = np.arctan2(y, x)
      return (r, theta, phi)

class _Ignitor(object):
   """Class representing a spatially coherent volume of initing fluid."""
   ##Shared class data

   ##Constructor
   def __init__(self, arr):
      """_Ignitor constructor.
      
      self    --> reference to this instance of _Ignitor (implicitly passed)
      arr     --> array of _Hotspots in this _Ignitor volume"""
      import numpy as np

      self.com = _Ignitor._calcCom(arr)
      self.rho_avg = _Ignitor._calcRavg(arr)
      self.M = _Ignitor._calcMtot(arr)

      x = self.com[0]
      y = self.com[1]
      z = self.com[2]
      self.r = np.sqrt(x**2 + y**2 + z**2)
      self.theta = np.arccos(z/self.r)
      self.phi = np.arctan2(y, x)

   ##Public methods

   ##Private methods
   @staticmethod
   def _calcCom(arr):
      """Return the center of mass of an array of _Hotspot's."""

      comx = 0.0
      comy = 0.0
      comz = 0.0
      m_tot = 0.0
      for hs in arr:
         m = hs.dx**3 * hs.rho
         m_tot += m
         x, y, z = hs.loc[0][0], hs.loc[0][1], hs.loc[0][2]
         comx += x * m
         comy += y * m
         comz += z * m
      comx = comx/m_tot
      comy = comy/m_tot
      comz = comz/m_tot

      return (comx, comy, comz)

   @staticmethod
   def _calcRavg(arr):
     """Return the average density (rho) of an array of _Hotspot's."""
     import numpy as np

     rho_tot = 0.0
     for hs in arr:
        rho_tot += hs.rho
     rho_avg = rho_tot/len(arr)

     return rho_avg

   @staticmethod
   def _calcMtot(arr):
     """Return the total mass of an array of _Hotspot's."""

     m_tot = 0.0
     for hs in arr:
        m = hs.dx**3 * hs.rho
        m_tot += m

     return m_tot


class _Globals(object):
   """struct-like class containing a simulation's globals, e.g. T_peak, r(T_peak), ..."""
   Tpeak = None #peak temperature
   tploc = (None, None, None, None) #x, y, z, r
   tpvel = (None, None, None, None) #vx, vy, vz, vr
   epeak = None #peak enucdot
   eploc = (None, None, None, None) #x, y, z, r

   def __str__(self):
      ret = 'Tpeak: ' + str(self.Tpeak) + '\n'
      ret += '  loc (x,y,z,r): ' + str(self.tploc) + '\n'
      ret += '  vel (vx,vy,vz,vr): ' + str(self.tpvel) + '\n'
      ret += 'peak enucdot: ' + str(self.epeak) + '\n'
      ret += '  loc (x,y,z,r): ' + str(self.eploc) + '\n'
      return ret

class SCPltfile(object):
   """A class representing a pltfile from a sub-Chandra simulation."""
   ##Shared class data##

   ##Constructor##
   def __init__(self, stage_dir, scratch_dir, label):
      """self        --> implicitly passed reference to this instance of SCSimulation
         stage_dir   --> staging directory containing all simulations
         scratch_dir --> scratch directory where the work is done but data's purged
         label       --> this simulation's label (e.g. 12050-107-175-3levs)"""

   ##Public methods##
   def sample(self):
      """Check the scratch directory for any updated output.  If found, copy to
      the stage directory."""
      #TODO: Implement this
      pass

   ##Private methods##
   def sampleprv(self):
      """Check the scratch directory for any updated output.  If found, copy to
      the stage directory."""
      #TODO: Implement this
      pass

#################
### Functions ###
#################

