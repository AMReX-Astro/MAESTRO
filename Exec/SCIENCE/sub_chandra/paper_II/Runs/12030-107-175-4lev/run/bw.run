#!/bin/bash
# Name the job
#PBS -N 12030-107-175-4lev
# Specify the project you're associated with (not technically required on BW)
#PBS -A jni
# email settings (send email when job's (a)borted, (b)egins, and/or (e)nds
#PBS -m a
#PBS -M adam.jacobs@stonybrook.edu
# Merge stdout and stderr
#PBS -j oe
# Select queue
#PBS -q normal
# Specify resources to request
#PBS -l walltime=04:00:00,nodes=128:ppn=32:xe

# this script runs with 8 threads, 4 MPI tasks/node, and 128 nodes on Blue Waters
#

# NOTE: To add certain modules that you do not have added via ~/.modules,
# you can do, e.g.:
#. /opt/modules/default/init/bash
#module load craype-hugepages2M  perftools

export PSC_OMP_AFFINITY=FALSE
export OMP_NUM_THREADS=8

# These are needed for very large MAESTRO runs, such as
# inputs_3d.2304.5dr.eq.dx_3levels with 6 threads
#export MPICH_PTL_OTHER_EVENTS=16384
#export MPICH_UNEX_BUFFER_SIZE=100000000
#export MPICH_MAX_SHORT_MSG_SIZE=10000

cd $PBS_O_WORKDIR

# run the compression script to tar up the plot and checkpoint files
# as they are created.
./process.bw &
PID=$!
trap 'kill -s TERM $PID' EXIT TERM HUP XCPU KILL

# find the latest restart file -- first look for one with 6 digits then fall
# back to 5
restartFile=$(find . -maxdepth 1 -type d -name "*chk??????" -print | sort | tail -1)

# the Header is the last thing written -- check if it's there, otherwise,
# fall back to the second-to-last check file written
if [ ! -f ${restartFile}/Header ]; then
  # how many *chk?????? files are there? if only one, then skip
  nl=$(find . -maxdepth 1 -type d -name "*chk??????" -print | sort | wc -l)
  if [ $nl -gt 1 ]; then
	  restartFile=$(find . -maxdepth 1 -type d -name "*chk??????" -print | sort | tail -2 | head -1)    
  else
	  restartFile=""
  fi
fi

# if the above checks failed, then there are no valid 6-digit chk files, so
# check the 5-digit ones
if [ "${restartFile}" = "" ]; then
  restartFile=$(find . -maxdepth 1 -type d -name "*chk?????" -print | sort | tail -1)

  # make sure the Header was written, otherwise, check the second-to-last
  # file
  if [ ! -f ${restartFile}/Header ]; then
    # how many *chk????? files are there? if only one, then skip
    nl=$(find . -maxdepth 1 -type d -name "*chk?????" -print | sort | wc -l)
    if [ $nl -gt 1 ]; then
	    restartFile=$(find . -maxdepth 1 -type d -name "*chk?????" -print | sort | tail -2 | head -1)    
    else
	    restartFile=""
    fi
  fi
fi

# cut out the numerical part of the *chkXXXXX file, here we use the 'k' in 
# 'chk' as the delimiter
restartNum=`echo ${restartFile} | cut -d'k' -f2`  

# restartString will be empty if no chk files are found -- i.e. new run
if [ "${restartNum}" = "" ]; then
    restartString=""
else
    restartString="--restart ${restartNum}"
    echo "Restarting with: " ${restartString}
fi

# Blue Waters has 26864 nodes (22640 Cray XE nodes + 4224 Cray XK nodes) 
# Each XE node has:
#   2 AMD 6276 'interlagos' processors
#   = 8x2=  16 Bulldozer modules (also called computer units or compute cores)
#   = 2x8x2=32 integer cores / integer scheduling units (this is smallest unit that a thread can exist on)
#   = 4 NUMA nodes (4 Bulldozer modules = 1 NUMA (shared memory) node)
#
# Each XK node has:
#   ?
# 
# -n  is the total number of MPI tasks
# -S  is the number of MPI tasks per NUMA/logical node 
# -N  is the number of MPI tasks per physical node 
#       (used if MPI tasks encompass multiple NUMA nodes. -S can't be < 1)
# -d  is the number of OpenMP threads per MPI task (must match OMP_NUM_THREADS)
# -ss forces MPI tasks to only allocate memory in their local NUMA node.
#   This can boost performance by preventing costly remote memory I/O, though 
#   it also restricts the amount of memory available to MPI tasks.

aprun -n 512 -S 1 -d $OMP_NUM_THREADS -ss ./main.Linux.Cray.mpi.omp.exe inputs3d.256.5dr.eq.dx_4levs ${restartString}

rm -f process.pid
