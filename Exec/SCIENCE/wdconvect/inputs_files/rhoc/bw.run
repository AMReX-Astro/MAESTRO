#!/bin/ksh
#PBS -N wdc-rhoc2
#PBS -j oe
#PBS -l walltime=00:20:00
#PBS -l nodes=256:ppn=32:xe
###PBS -q batch

cd $PBS_O_WORKDIR


# Titan has 18688 physical nodes, each of which has 16 cores and 2 NUMA nodes
# 
# -n  is the total number of MPI tasks (should be nodes*-S*2)
# -S  is the number of MPI tasks per NUMA node 
# -d  is the number of OpenMP threads per MPI task (must match OMP_NUM_THREADS)
# -ss forces MPI tasks to only allocate memory in their local NUMA node.
#   This can boost performance by preventing costly remote memory I/O, though 
#   it also restricts the amount of memory available to MPI tasks.


# start a script that tars the output as it is produced
./process.bw &
PID=$!
trap 'kill -s TERM $PID' EXIT TERM HUP XCPU KILL


# automated restart

# check the 6-digit ones              
restartFile=$(find . -type d -name "*chk??????" -print | sort | tail -1)

# the Header is the last thing written -- check if it's there, otherwise,                                
# fall back to the second-to-last check file written                                                     
if [ ! -f ${restartFile}/Header ]; then

    # how many *chk?????? files are there? if only one, then skip                                        
    nl=$(find . -type d -name "*chk??????" -print | sort | wc -l)
    if [ $nl -gt 1 ]; then
        restartFile=$(find . -type d -name "*chk??????" -print | sort | tail -2 | head -1)
    else
        restartFile=""
    fi
fi


# if the above checks failed, then there are no valid 6-digit chk files, so                              
# check the 5-digit ones                                                                                 
if [ "${restartFile}" = "" ]; then
    restartFile=$(find . -type d -name "*chk?????" -print | sort | tail -1)

    # make sure the Header was written, otherwise, check the second-to-last                              
    # file                                                                                               
    if [ ! -f ${restartFile}/Header ]; then
        restartFile=$(find . -type d -name "*chk?????" -print | sort | tail -2 | head -1)
    fi
fi


# cut out the numerical part of the *chkXXXXX file, here we use the                                      
# 'k' in 'chk' as the delimiter                                                                          
restartNum=`echo ${restartFile} | cut -d'k' -f2`


# restartString will be empty if no chk files are found -- i.e. new run                                  
if [ "${restartNum}" = "" ]; then
    restartString=""
else
    restartString="--restart ${restartNum}"
fi


export OMP_NUM_THREADS=16

aprun -n 512 -N 2 -d 16 ./main.Linux.Cray.mpi.omp.exe inputs_3d.512.5dr.rhoc2 ${restartString}


rm -f process.pid
