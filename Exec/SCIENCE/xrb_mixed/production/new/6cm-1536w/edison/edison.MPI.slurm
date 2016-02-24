#!/bin/ksh
#SBATCH -J xrb-simplenet
#SBATCH -A m1938
#SBATCH -N 11                
#SBATCH -p regular
#SBATCH -t 04:00:00

# -N = number of nodes


# this script runs with 12 threads on edison -- this seems to give the best 
# performance.

export OMP_NUM_THREADS=1


cd $SLURM_SUBMIT_DIR

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

# Edison has 2 Intel "IvyBridge" processors each with 12 cores.  Each processor
# connects to 1/2 of the node's memory, and is called a NUMA node.  So there
# are 2 NUMA nodes per Edison node.

srun -n 256 ./main.Linux.Intel.mpi.exe inputs_2d_6.0cm ${restartString}


rm -f process.pid

