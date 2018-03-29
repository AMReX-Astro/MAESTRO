#!/bin/ksh
#SBATCH -N 60
#SBATCH -C knl,quad,cache
#SBATCH -q regular
#SBATCH -S 4
#SBATCH -J speceq_rc3.5_Tc3.0
#SBATCH --mail-user=eugene.willcox@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 12:00:00
#SBATCH -A m3018

cd $SLURM_SUBMIT_DIR

## Restart Setup
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

## Parallel Configuration
export OMP_NUM_THREADS=64
export OMP_PLACES=threads
export OMP_PROC_BIND=true

# -n is the number of MPI tasks
# -c sets the number of logical threads per MPI task
# In the line below I requested 60 KNL nodes and am running 60 MPI tasks
# Each MPI task gets 256 logical (64 physical) threads
srun -n 60 -c 256 --cpu_bind=cores ./main.Linux.Intel.mpi.omp.exe inputs_3d_4level ${restartString}
