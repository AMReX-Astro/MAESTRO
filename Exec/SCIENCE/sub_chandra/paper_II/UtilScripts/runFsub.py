#!/usr/bin/python

# Simple Python script to execute 
# ~/Codebase/AmrPostprocessing/F_Src/MAESTRO_sub_chandra/fsubchandra*
# on a set of plotfiles and save the output to the staging directory
import subchandra as sc
import sys
import getopt  #If you have Python 2.7+, the 'argparse' module is a better solution
from os.path import isdir, join

def main():
  #Parse options and args.  
  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], 't:c:d:f', ['stage=', 'scratch=', 'do-hotspots=', 'full-star'])
  except getopt.GetoptError as err:
    #Print help, exit
    print str(err)
    usage()
    sys.exit(2)
  
  #Error check
  if not args:
    print 'No run label given!'
    usage()
    sys.exit(2)
  #if len(args) > 1:
  #  print 'Unexpected arguments given!'
  #  usage()
  #  sys.exit(2)
  
  #Set default values
  stage_dir = '/u/sciteam/ajacobs/Projects/SubChandra/Runs/'
  scratch_dir = '/u/sciteam/ajacobs/scratch/sub_chandra/'
  hpcnt = -1.0
  fstar = False
  
  #Read in parsed opts/args
  for o, v in opts:
    if o in ('-t', '--stage'):
      stage_dir = v
    elif o in ('-c', '--scratch'):
      scratch_dir = v
    elif o in ('-d', '--do-hotspots'):
      hpcnt = float(v)
    elif o in ('-f', '--full-star'):
      fstar = True
    elif o in ('-h', '--help'):
      usage()
      sys.exit(0)
    else:
      assert False, 'unrecognized option!'

  for a in args:
    run_label = a
    
    #Make sure run is valid
    if not (isdir(join(stage_dir, run_label)) and isdir(join(scratch_dir, run_label))):
      assert False, "The run you specified doesn't exist!"
    
    simout = sc.SCOutput(stage_dir, scratch_dir, run_label, None)
    simout.analyzePlts(hpcnt=hpcnt, full_star=fstar)

def usage():
  print 'usage:'
  print './runFsub [-t|--stage <staging directory>] [-c|--scratch <scratch directory>] <run label>'
  print ' '
  print 'Arguments: '
  print '  <run label>: The label for the run you want to analyze the plotfiles of, e.g. 12050-107-175-3lev'
  print 'Options: '
  print '  <staging directory>: The directory used to stage runs (default: /ccs/home/ajacobs/Projects/SubChandra/Runs/)'
  print '  <scratch directory>: The directory where runs are executed and plotfiles are dumped, but gets regularly purged (default: /tmp/work/ajacobs/sub_chandra/)'

if __name__ == "__main__":
  main()
