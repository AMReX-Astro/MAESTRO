#!/bin/bash

# Name the job
#PBS -N fsub_analysis
# Specify the project you're associated with (not technically required on BW)
#PBS -A jni
# Merge stdout and stderr
#PBS -j oe
# Select queue
#PBS -q normal
# Specify resources to request
#PBS -l walltime=02:00:00,nodes=1:ppn=32:xe

