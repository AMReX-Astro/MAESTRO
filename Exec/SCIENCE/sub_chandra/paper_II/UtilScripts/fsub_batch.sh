#!/bin/bash

# Name the job
#PBS -N fsub_analysis
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
#PBS -l walltime=01:00:00,nodes=1:ppn=32:xe


files=("10030-108-185-4lev" \
       "10040-108-185-4lev" \
       "11020-108-185-4lev" \
       "11030-108-185-4lev" \
       "12020-107-175-4lev" \
       "12020-108-175-4lev" \
       "12030-107-175-4lev" \
       "12030-108-175-4lev" )

for f in "${files[@]}"
do
   python runFsub.py "$f" >> "$f".out &
done

wait
