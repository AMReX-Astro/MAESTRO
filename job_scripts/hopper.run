#!/bin/ksh
#PBS -A m1400
#PBS -N xrb-6cm-wt-CNO-hot
#PBS -j oe
#PBS -q regular
#PBS -l walltime=3:00:00
#PBS -l mppwidth=128


# this script runs with 6 threads on hopper -- this seems to give the best 
# performance.

export OMP_NUM_THREADS=6


cd $PBS_O_WORKDIR

# run the compression script to tar up the plot and checkpoint files
# as they are created.
./process.xrb &
PID=$!
trap 'kill -s TERM $PID' EXIT TERM HUP XCPU KILL

# find the latest restart file -- first look for one with 7 digits then fall
# back to 6 and then 5
restartFile=$(find . -type d -name "*chk???????" -print | sort | tail -1)

# the Header is the last thing written -- check if it's there, otherwise,
# fall back to the second-to-last check file written
if [ ! -f ${restartFile}/Header ]; then

    # how many *chk?????? files are there? if only one, then skip
    nl=$(find . -type d -name "*chk???????" -print | sort | wc -l)
    if [ $nl -gt 1 ]; then
	restartFile=$(find . -type d -name "*chk???????" -print | sort | tail -2 | head -1)    
    else
	restartFile=""
    fi
fi

# if the above checks failed, then there are no valid 7-digit chk files, so
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


# -n is the total number of MPI tasks
# -N is the number of MPI tasks PER compute node
# -d is the number of OpenMP threads per MPI task (must match OMP_NUM_THREADS)
#
# see http://www.nersc.gov/users/computational-systems/hopper/running-jobs/using-openmp-with-mpi/

# you need to change -n to give the number of MPI tasks you want.  Then
# that number * 6 should match the mppwidth set at the top of this script
echo "restart string: " ${restartString}

aprun -n 128 ./main.Linux.Cray.mpi.exe inputs_2d_6.0cm.toy.wide.tall.hot ${restartString}


rm -f process.pid

