#!/bin/ksh
#SBATCH -A m1938
#SBATCH -J process-yt
#SBATCH -N 1
#SBATCH -p debug
#SBATCH -t 00:10:00

export OMP_NUM_THREADS=24

cd $SLURM_SUBMIT_DIR

srun -n 1 -c ${OMP_NUM_THREADS} ./vol.py xrb_3d_simple_plt48398



