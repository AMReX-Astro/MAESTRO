#!/bin/ksh
#SBATCH -M escori
#SBATCH -J speceq_rc3.5_hpss-xfer
#SBATCH -A m3018
#SBATCH -p xfer
#SBATCH -t 12:00:00

cd $SLURM_SUBMIT_DIR

# do our archiving
pidfile=process.pid

./process.xrb

PID=$!
trap 'kill -s TERM $PID' EXIT TERM HUP XCPU KILL

rm -f process.pid
