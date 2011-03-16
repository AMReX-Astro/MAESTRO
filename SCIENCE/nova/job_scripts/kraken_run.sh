#!/bin/bash
#PBS -A TG-AST100037
#PBS -l size=36,walltime=24:00:00
#PBS -N resolution_stepdown
#PBS -o output/$PBS_JOBNAME.$PBS_JOBID.out
#PBS -j oe
#PBS -m a
#PBS -m e
#PBS -M brendan.krueger@stonybrook.edu
#PBS -q batch

#### -e output/$PBS_JOBNAME.$PBS_JOBID.err

#==============================================================================
# Set up ----------------------------------------------------------------------

# Assign user-definable values
inputs_file=inputs_2d_resolution_stepdown  # name of MAESTRO input file
Nsteps_per_run=70                  # maximum number of steps per run

# Parse command-line arguments
# - "--debug"
if [ "$1" = "--debug" ]; then
   debug=1
else
   debug=0
fi

# Set name of script output file
if [ "$debug" -eq 1 ]; then
   script_out=$0.test
else
   script_out=output/$PBS_JOBNAME.$PBS_JOBID.run
fi

# Change directory
if [ "$debug" -eq 1 ]; then
   echo "DEBUG MODE" > ${script_out}
else
   cd $PBS_O_WORKDIR
   echo "Move to $PBS_O_WORKDIR" > ${script_out}
fi

# Source the external functions
. ./kraken_functions.sh

#==============================================================================
# Grab MAESTRO output file prefixes -------------------------------------------

# Checkpoint files
chk_base=`get_chk_prefix $inputs_file`
if [ "${chk_base}" = "" ]; then
   chk_base=chk #MAESTRO's default value
fi

# Plot files
plt_base=`get_plt_prefix $inputs_file`
if [ "${plt_base}" = "" ]; then
   plt_base=plt #MAESTRO's default value
fi

#==============================================================================
# Find and process last checkpoint file ---------------------------------------

# Search for most-recent checkpoint file
last_chk=`get_last_chk $chk_base`

# If checkpoint file wasn't found
if [ "${last_chk}" = "" ]; then

   # Build a blank restart string
   restart_string=""
   echo "No previous checkpoint file found: starting new simulation." >> ${script_out}
else

   # Exit if lastcheck.step_num >= input_file.max_steps
   last_chk_num=${last_chk:${#chk_base}}
   max_step=`get_max_step $inputs_file`
   if [ $last_chk_num -ge $max_step ]; then
      echo "Exceeded maximum number of steps ($last_chk_num >= $max_step)" >> ${script_out}
      exit
   else
      echo "At step $last_chk_num of $max_step" >> ${script_out}
   fi

   # Exit if lastcheck.time >= input_file.max_time
   last_chk_time=`head -3 ${last_chk}/Header | tail -1 | sed 's_ *\([0-9\.]*\)[eEdD]\(.*\)_\1 * 10 ^ (0\2)_' | bc -l`
   max_time=`get_max_time $inputs_file`
   time_check=`echo "$last_chk_time >= $max_time" | bc -l`
   if [ $time_check -eq 1 ]; then
      echo "Exceeded maximum time ($last_chk_time >= $max_time)" >> ${script_out}
      exit
   else
      echo "At time $last_chk_time of $max_time" >> ${script_out}
   fi

   # Construct restart string
   restart_string="--restart ${last_chk_num}"
fi

#==============================================================================
# Construct additional MAESTRO options ----------------------------------------

# Construct args.max_steps from max_steps_per_run and lastcheck.step_num
if [ "$last_chk_num" = "" ]; then
   max_step_string=$Nsteps_per_run
else
   max_step_string=`echo "$last_chk_num + $Nsteps_per_run" | bc -l`
fi

# Construct max_steps string
max_step_string="--max_step $max_step_string --chk_int $Nsteps_per_run"

#==============================================================================
# Run MAESTRO -----------------------------------------------------------------

if [ "$debug" -eq 1 ]; then
   echo "DEBUG: run MAESTRO with command-line options of:" >> ${script_out}
   if [ "$restart_string" != "" ]; then
      echo "       $restart_string" >> ${script_out}
   fi
   echo "       $max_step_string" >> ${script_out}
else
   aprun -n 36 ./main.Linux.PathScale.mpi.exe ${inputs_file} ${restart_string} ${max_step_string}
fi

