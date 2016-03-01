#!/usr/bin/env python

#Title:         Titan Python Module (titan.py)
#Author:        Adam Jacobs
#Creation Date: Aug 31, 2015

#Usage: load as a module

#Description: This module implements concrete classes of the abstract classes found
#             in the `supercomputer` module.

#TODO: 
# 1) 
# 2) 

###########################################
### Global Imports, Data, and Constants ###
###########################################
from supercomputer import HPSS, Supercomputer, SupercomputerConfig


###############
### Classes ###
###############
class TitanSC(Supercomputer):
   """Concrete subclass of Supercomputer representing an implementation of the Titan system."""
   
   ##Constructor##
   def __init__(self, myconfig=None):
      """self          --> implicitly passed reference to this instance of TitanSC
         myconfig      --> Configuration of this machine"""
      if myconfig is None:
         myconfig = SupercomputerConfig
         myconfig.default_queue = 'batch'
         myconfig.default_accnt = 'ast106'
         super(TitanSC, self).__init__(myconfig)
         
         self._initProjAndScratch()
         self.myconfig.proc_script = 'process.titan'
         self.myconfig.run_script = 'titan.run'
         self.myconfig.npn = 2
      elif isinstance(myconfig, SupercomputerConfig):
         super(TitanSC, self).__init__(myconfig)
      else:
         raise TypeError('myconfig argument is not of type SupercomputerConfig')

   ##Public Methods##
   def getProjsAndScratch(self):
      """Return a tuple of (projects base directory, scratch directory) for the current host supercomputer."""
      return (self.myconfig.projs_base, self.myconfig.scratch_base)

   def generateJobScript(self, outpath, jobconfig, template_path=None):
      """Generate a PBS script to submit to the queue
      
      outpath       --> path to save the new job script in
      jobconfig     --> A supercomputer.JobConfig object for the new job
      template_path --> Full path to the template file for the job script"""
      from os.path import isfile, join
      from os import chmod
      from datetime import time
      import supercomputer as scomp
      
      #Get job config dictionary
      if not jobconfig.accnt:
         jobconfig.accnt = self.myconfig.default_accnt
      if not jobconfig.queue:
         jobconfig.queue = self.myconfig.default_queue
      if not isinstance(jobconfig.time, time):
         raise ValueError('jobconfig.time must be of type datetime.time!')
      job_dict = scomp.JobConfig.toPBSDict(jobconfig)

      #Write the new job script 
      if template_path is None:
         raise NotImplementedError("Titan job script generation without template not implemented.")
      elif isfile(template_path):
         outfile = join(outpath, 'titan.run')
         #TODO: add error checking to make sure outfile is valid
         with open(outfile, 'w') as ofile, open(template_path, 'r') as tfile:
            for line in tfile:
               newline = line
               for key in job_dict:
                  if(line.lstrip().startswith(key)):
                     newline = key + job_dict[key] + '\n'
               ofile.write(newline)
         chmod(outfile, 0o755)  #Everyone can read and execute, only user can write
      else:
         raise ValueError("Invalid template_path sent: {}".format(template_path))

   def isQueued(self, label):
      """Is this job in the queue?
      
      label         --> name of job you are checking for"""

      return self.getQStatus(label)[0].count('inactive') == 0

   def getQStatus(self, label):
      """Get a tuple of (queue status string, the number of times this job is in the queue)
      
      label         --> name of job you are checking for"""
      from supercomputer import TermColors
      job_dict = self._getQJobDict()

      #Determine queue status and count
      queue_status = 'inactive'
      qcount = 0
      if label in job_dict:
         for id, prop_dict in job_dict[label].iteritems():
            if prop_dict['job_state'] in ['H', 'Q']:
               if queue_status != 'running!':
                  queue_status = 'queued'
               qcount += 1
            if prop_dict['job_state'] == 'R':
               qcount += 1
               if queue_status == 'running!':
                  print('WARNING! You have multiple instances of {} running!'.format(label))
               queue_status = 'running!'
   
      #Colorize
      startcol = None
      endcol = TermColors.RESET
      if queue_status == 'inactive':
         startcol = TermColors.START_RED
      elif queue_status == 'queued':
         startcol = TermColors.START_BLUE
      elif queue_status == 'running!':
         startcol = TermColors.START_GREEN
      queue_status = startcol + queue_status + endcol

      return (queue_status, qcount)

   def printUsageOverview(self):
      """Print an overview of the current user's usage including current queue status, 
      allocation's hours used, HPSS usage."""
      import getpass
      
      print('=== Queue Overview ===')
      print('- label ------------------ status --------- qcount -------------------')
      job_dict = self._getQJobDict()
      for job in job_dict:
         stat, count = self.getQStatus(job)
         count_len = 29 - len(stat)
         print('{:27.27s}'.format(job) + stat + '{:{width}d}'.format(count, width=count_len))
      print(' ')
      
      print('=== Memory Overview ===')
      print('----------------------------------------------------------------------')
      #Get user home info
      mem_dict = self._getMemDict()
      print('Not implemented yet, getting bad data from Titan')
      #tokens = out[7].split()
      #pct = float(tokens[1][:-1])/float(tokens[2][:-1]) #The [:-1] strips the 'G' from the string
      #pct = '{0:.2%}'.format(pct)
      #hm_used = tokens[1]
      #hm_total = tokens[2]
      #
      ##Get HPSS info
      ##out = !showusage -s hpss
      #tokens = out[25].split()
      #hpss_total = tokens[4] #Project quota is currently 50 TB
      #hpss_pct = float(tokens[2][:-1])/float(hpss_total[:-1])
      #hpss_pct = '{0:.2%}'.format(hpss_pct)
      #hpss_used = tokens[2]
      ##Print
      #print('User Home: ' + hm_used + '/' + hm_total + ' (' + pct + ' used)') 
      #print('Proj HPSS: ' + hpss_used + '/' + hpss_total + ' (' + hpss_pct + ' used)') 
      print(' ')
      print('=== Node-hour Overview ===')
      print('- user/proj -------------- usage (Mhr) ---- total/pcnt ---------------')
      hour_dict = self._getHourDict()
      current_user = getpass.getuser()
      current_proj = self.getAccnt()
      hour_tup = hour_dict[current_proj]
      print('{:26.26s} {:5.2f} {:16.2f}'.format(current_proj, hour_tup[0], hour_tup[1]))
      del hour_dict[current_proj]
      hour_tup = hour_dict[current_user]
      print('{:26.26s} {:5.2f} {:17.2%}'.format(current_user, hour_tup[0], hour_tup[0]/hour_tup[1]))
      del hour_dict[current_user]
      print(' ')
      for key in hour_dict:
         hour_tup = hour_dict[key]
         print('{:26.26s} {:5.2f} {:17.2%}'.format(key, hour_tup[0], hour_tup[0]/hour_tup[1]))
   
  
   def getAccnt(self):
      """Get the name of the current allocation/account"""
      return self.myconfig.default_accnt


   ##Private Methods##
   def _initProjAndScratch(self):
      import getpass

      current_user = getpass.getuser()
      projs_base   = '/ccs/home/{}/Projects/'.format(current_user)
      scratch_base = '/lustre/atlas/scratch/{}/{}/'.format(current_user, self.myconfig.default_accnt)
      
      self.myconfig.projs_base = projs_base
      self.myconfig.scratch_base = scratch_base

   def _getQStatOut(self):
      """Return the output from a `qstat -f1 -u <current user>` call."""
      from subprocess import call, Popen, PIPE, STDOUT
      from shlex import split
      import getpass

      current_user = getpass.getuser()
      args = split('qstat -f1 -u {}'.format(current_user))
      qproc = Popen(args, stdout=PIPE, stderr=PIPE)
      qproc_out, qproc_err = qproc.communicate()
      #TODO: Add error checking

      return qproc_out

   def _getQJobDict(self):
      """Return a job dictionary that relates a job's label to its various
      instances in the queue and each instance's properties.  Based on output 
      from a `qstat -f1 -u <current user>` call.
      
      ret = {'job_name': jobid_dict}
      jobid_dict = {jobid: properties_dict}
      properties_dict = {'job property': 'val'}"""
      from subprocess import call, Popen, PIPE, STDOUT
      from shlex import split
      import getpass
      #TODO: Cache the dictionary?  Avoids repeated calls to qstat.

      current_user = getpass.getuser()
      args = split('qstat -f1 -u {}'.format(current_user))
      qproc = Popen(args, stdout=PIPE, stderr=PIPE)
      qproc_out, qproc_err = qproc.communicate()
      #TODO: Add error checking

      found_new_job = False
      job_dict = {}
      cur_id_dict = {}
      cur_properties_dict = {}
      for line in qproc_out.splitlines():
         if found_new_job:
            tokens = line.partition('=')
            key = tokens[0].strip()
            val = tokens[2].strip()
            if tokens[1] == '':
               #Didn't find '=' in the string
               found_new_job = False
            else:
               cur_properties_dict[key] = val
            if key == 'Job_Name':
               if val not in job_dict:
                  job_dict[val] = {}
               job_dict[val][cur_id] = cur_properties_dict

         if line.startswith('Job Id:'):
            found_new_job = True
            cur_id = line.partition(':')[2].strip()
            cur_properties_dict = {}

      return job_dict

   def _getMemDict(self):
      """Return a memory dictionary with information about current memory use
      Based on output from `quota -Qs` and `showusage -s hpss`
      
      ret = {'user_home' | 'proj_hpss': mem_tuple}
      mem_array = (<usage in GB>, <quota in GB>)"""
      from subprocess import call, Popen, PIPE, STDOUT
      from shlex import split
      import getpass

      pass
      #current_user = getpass.getuser()
      #args = split('qstat -f1 -u {}'.format(current_user))
      #qproc = Popen(args, stdout=PIPE, stderr=PIPE)
      #qproc_out, qproc_err = qproc.communicate()
      ##TODO: Add error checking

      #found_new_job = False
      #job_dict = {}
      #cur_id_dict = {}
      #cur_properties_dict = {}
      #for line in qproc_out.splitlines():
      #   if found_new_job:
      #      tokens = line.partition('=')
      #      key = tokens[0].strip()
      #      val = tokens[2].strip()
      #      if tokens[1] == '':
      #         #Didn't find '=' in the string
      #         found_new_job = False
      #      else:
      #         cur_properties_dict[key] = val
      #      if key == 'Job_Name':
      #         if val not in job_dict:
      #            job_dict[val] = {}
      #         job_dict[val][cur_id] = cur_properties_dict

      #   if line.startswith('Job Id:'):
      #      found_new_job = True
      #      cur_id = line.partition(':')[2].strip()
      #      cur_properties_dict = {}

      #return mem_dict

   def _getHourDict(self):
      """Return a dictionary of data on allocation CPU-hour usage.
      Based on output from a `showusage -f` call.
      
      ret = {'username' | 'projid': hour_tuple}
      hour_tuple = (<CPU hours used in Mhrs>, <CPU allocation in Mhrs>)"""
      from subprocess import call, Popen, PIPE, STDOUT
      from shlex import split
      import getpass

      current_user = getpass.getuser()
      current_proj = self.getAccnt()
      args = split('showusage -f')
      suproc = Popen(args, stdout=PIPE, stderr=PIPE)
      suproc_out, suproc_err = suproc.communicate()
      #TODO: Add error checking

      hour_dict = {}
      for line in suproc_out.splitlines():
         if line.count(current_proj) > 0:
            tokens = line.split()
            proj_total = float(tokens[1])/1.e6
            proj_usage = float(tokens[3])/1.e6
            hour_dict[current_proj] = (proj_usage, proj_total)
         if line.strip().endswith('%'):
            tokens = line.split()
            username = tokens[0]
            user_usage = float(tokens[1])/1.e6
            hour_dict[username] = (user_usage, proj_total)

      return hour_dict


class TitanHPSS(HPSS):
   """Concrete subclass of HPSS representing an implementation of Blue Waters' HPSS system."""

   ##Constructor##
   def __init__(self, hpss_base_dir):
      """self          --> implicitly passed reference to this instance of HPSS
         hpss_base_dir --> the base HPSS directory all other files and directories 
                           will be relative to"""
      super(TitanHPSS, self).__init__(hpss_base_dir)
      #self.local_ep     = "ncsa#BlueWaters"     #ep = endpoint
      #self.storage_ep     = "ncsa#Nearline"
      #self.go_host    = "cli.globusonline.org"
      #self.sync_level = "3" # 0: Copy files that do not exist at the destination
      #                      # 1: Copy files if the destination size doesn't match
      #                      #      source size
      #                      # 2: Copy files if destination timestamp is older than
      #                      #      source
      #                      # 3: Copy files if source and destination checksums 
      #                      #      don't match
   
   ##Public Methods##
   def sendToHPSS(self, item):
      """Send self.loc_dir/item (a file or directory) 
      to the HPSS directory self.hpss_base_dir/self.loc_base/item"""

      ##Generate a task id
      #task_id = self._genID()

      ##Execute the transfer
      #self._transfer(task_id, item)

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
      #from subprocess import Popen, PIPE, STDOUT
      #import shlex

      #TRY_MAX = 4
      #tries = 0
      #tid = None
      #serr = ''

      #while tries < TRY_MAX:
      #   #Create process to generate id
      #   #Command to generate id:   task_id=$(ssh $go_host transfer --generate-id)
      #   gid_command = 'ssh ' + self.go_host + ' transfer --generate-id'
      #   args = shlex.split(gid_command)
      #   if DEBUG:
      #      print 'executing: ', gid_command
      #   gid_proc = Popen(args, stdout=PIPE, stderr=PIPE)

      #   #Check the output
      #   (sout, serr) = gid_proc.communicate()
      #   retcode = gid_proc.poll()
      #   if retcode == 0:
      #      tid = sout.replace('\n','')
      #      break
      #   else:
      #      tries += 1

      #if tries >= TRY_MAX:
      #   print '_genID() error output: ', serr
      #   raise RuntimeError, 'Globus failed to generate a task id! Make sure all endpoints are activated'

      #return tid

   def _transfer(self, task_id, item):
      """Use the Globus CLI to transfer item using the given task"""
      from os.path import isdir, isfile
      from subprocess import Popen, PIPE, STDOUT
      import shlex
      from warnings import warn

      #TRY_MAX = 4
      #SUCCESS = 0
      #SSH_ERROR = 255
      #tries = 0
      #
      ##Is item a file or directory?
      #dir_flag = ''
      #if isdir(item):
      #   dir_flag = '-r'
      #   if not item.endswith('/'):
      #      item = item + r'/'  #Directories must end with a / according to Globus
      #else:
      #   if not isfile(item):
      #      raise ValueError, 'The argument item is invalid! Must be local file or directory'

      ##Command to transfer: ssh $go_host transfer --verify-checksum --taskid=$task_id --label=\"$task_label\" -s $sync_level   -- $src_dst [-r]
      ##Build src/dst string for the transfer command
      #src = self.local_ep + self.loc_dir + '/' + item
      #dst = self.storage_ep + self.hpss_base_dir + '/' + self.loc_base + '/' + item
      #src_dst = src + ' ' + dst
     
      ##Make a label (remove any illegal characters [./])
      #task_label = 'Archive ' + item
      #task_label = task_label.replace('.', '-')
      #task_label = task_label.replace('/', '')


      ##Build the transfer command string
      #trans_command = 'ssh ' + self.go_host + ' transfer --verify-checksum --taskid=' + task_id
      #trans_command += r' --label=\"' + task_label + r'\" -s ' + self.sync_level 
      #trans_command += ' -- ' + src_dst + ' ' + dir_flag

      ##Execute the transfer
      #if DEBUG:
      #   print 'calling: ', trans_command
      #args = shlex.split(trans_command)
      #while tries < TRY_MAX:
      #   trans_proc = Popen(args, stdout=PIPE, stderr=PIPE)
      #   (sout, serr) = trans_proc.communicate()
      #   retcode = trans_proc.poll()

      #   if retcode == SUCCESS:
      #      if DEBUG:
      #         print 'Transfer successfully established!'
      #      #Wait for the transfer to complete
      #      #TODO: Maybe we don't want to do this? For now I like it, makes sure transfers finish
      #      wait_command = 'ssh ' + self.go_host + ' wait ' + task_id
      #      args = shlex.split(wait_command)
      #      wait_proc = Popen(args)
      #      retcode = wait_proc.wait()
      #      if retcode == SUCCESS:
      #         return
      #      else:
      #         raise RuntimeError, 'Globus wait failed on task ' + task_id
      #   elif retcode == SSH_ERROR:
      #      tries += 1
      #      warn('ssh had a connection error, retry ' + str(tries) + ' of ' + str(TRY_MAX))
      #   else:
      #      print sout
      #      print serr
      #      raise RuntimeError, 'Globus transfer command failed fatally! Make sure endpoints are activated'

      #if tries >= TRY_MAX:
      #   raise RuntimeError, 'Globus failed to transfer ' + item + ' after ' + str(tries) + ' tries!'


#################
### Functions ###
#################

#################
### Execution ###
#################
#This is only for testing.  bluewaters.py is intended to be used as a module.
if __name__== "__main__":
   if(len(sys.argv) <= 1):
      #TODO: Add any desired testing
      pass

