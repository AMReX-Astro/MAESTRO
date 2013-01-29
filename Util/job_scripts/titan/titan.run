#!/bin/ksh
#PBS -A ast006
#PBS -N 10050-2lev
#PBS -j oe
#PBS -q batch
#PBS -l walltime=04:00:00,nodes=128

# this script runs with 8 threads per MPI task, 2 MPI tasks/node, and 128 nodes on titan 
#

export PSC_OMP_AFFINITY=FALSE
export OMP_NUM_THREADS=8

# These were needed for very large MAESTRO runs, such as
# inputs_3d.2304.5dr.eq.dx_3levels with 6 threads on jaguar.
# Not sure if they're needed on titan.
#export MPICH_PTL_OTHER_EVENTS=16384
#export MPICH_UNEX_BUFFER_SIZE=100000000
#export MPICH_MAX_SHORT_MSG_SIZE=10000

cd $PBS_O_WORKDIR

# run the compression script to tar up the plot and checkpoint files
# as they are created.
./process.titan &
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

# Titan has 18688 physical nodes, each of which has 16 cores and 2 NUMA nodes
# 
# -n  is the total number of MPI tasks (should be nodes*-S*2)
# -S  is the number of MPI tasks per NUMA node 
# -d  is the number of OpenMP threads per MPI task (must match OMP_NUM_THREADS)
# -ss forces MPI tasks to only allocate memory in their local NUMA node.
#   This can boost performance by preventing costly remote memory I/O, though 
#   it also restricts the amount of memory available to MPI tasks.

aprun -n 256 -S 1 -d 8 -ss ./main.Linux.Cray.mpi.omp.exe inputs_3d.256.5dr.eq.dx_2levs ${restartString}

rm -f process.pid
