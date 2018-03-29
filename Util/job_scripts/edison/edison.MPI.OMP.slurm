#!/bin/ksh
#SBATCH -q regular
#SBATCH -J urca_rc4.5_Tc6.0_Mc1
#SBATCH --mail-user=eugene.willcox@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -N 128
#SBATCH -t 8:00:00
#SBATCH -A m3018
#SBATCH -L SCRATCH,project

# -N = number of nodes


# this script runs with 12 threads on edison -- this seems to give the best 
# performance.

export OMP_NUM_THREADS=12
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

cd $SLURM_SUBMIT_DIR

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
srun -n 256 -c ${OMP_NUM_THREADS} --cpu_bind=cores main.Linux.Intel.mpi.omp.exe inputs_3d_4level_Tc6.0e8_rhoc4.5e9 ${restartString}
