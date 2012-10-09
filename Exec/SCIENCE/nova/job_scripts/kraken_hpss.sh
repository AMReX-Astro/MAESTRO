#!/bin/bash
#PBS -A TG-AST100037
#PBS -l walltime=1:00:00
#PBS -N resolution_stepdown
#PBS -o output/$PBS_JOBNAME.$PBS_JOBID.hpss
#PBS -j oe
#PBS -q hpss

#==============================================================================
# Set up ----------------------------------------------------------------------

# Assign user-definable values
inputs_file=inputs_2d_resolution_stepdown
HTAR=/usr/local/hsi/bin/htar
FTIME_EXE=./ftime.Linux.PathScale.exe

# Parse command-line arguments
if [ "$1" = "--debug" ]; then
   debug=1
else
   debug=0
fi

# Change directory
if [ "$debug" -eq 1 ]; then
   echo "DEBUG MODE"
else
   cd $PBS_O_WORKDIR
   echo "Move to $PBS_O_WORKDIR"
fi

# Source the external functions
. ./kraken_functions.sh

# Get prefixes
chk_prefix=`get_chk_prefix $inputs_file`
if [ "$chk_prefix" = "" ]; then
   echo "Error in \"get_chk_prefix()\"."
   chk_prefix=chk
fi
plt_prefix=`get_plt_prefix $inputs_file`
if [ "$plt_prefix" = "" ]; then
   echo "Error in \"get_plt_prefix()\"."
   plt_prefix=plt
fi
if [ "$debug" -eq 1 ]; then
   echo "Checkpoint prefix = \"$chk_prefix\""
   echo "Plotfile prefix = \"$plt_prefix\""
fi

# Set directory to archive to on HPSS
work_dir=`pwd`
HPSS_DIR=`basename $work_dir`

# Make the storage directories
#    once we process a file, we will move the plotfiles into the plotfiles/
#    directory.  This then hides them from the script, so if the system
#    later purges the files in the pltXXXXX directory and the .processed
#    file, we don't overwrite our archived data with a tarred empty
#    directory structure.  We do the same with the checkpoint files (using
#    checkfiles/)
if [ ! -d plotfiles ]; then
  mkdir plotfiles
fi
if [ ! -d checkfiles ]; then
  mkdir checkfiles
fi

#==============================================================================
# Archive diagnostics files ---------------------------------------------------

# archive any diagnostic files first -- give them a unique name, appending 
# the date string, to make sure that we don't overwrite anything
datestr=$(date +"%Y%m%d_%H%M_%S")
diag_files=$(find . -maxdepth 1 -name "*_diag.out" -print)
ftime_files=$(find . -maxdepth 1 -name "ftime.out" -print)

if [ "${diag_files}" ] || [ "${ftime_files}" ]; then
    ${HTAR} -cvf ${HPSS_DIR}/diag_files_${datestr}.tar ${diag_files} ${ftime_files} >> /dev/null
    echo "Diagnostics files archived as \"diag_files_$datestr.tar\""
fi

#==============================================================================
# Archive plot files ----------------------------------------------------------

# Take all but the final plt file -- we want to ensure they're completely 
# written to disk.  Strip out any tar files that are lying around as well 
# as pltXXXXX.processed files.
pltlist5=`ls -d ${plt_prefix}????? 2> /dev/null | sort`
pltlist6=`ls -d ${plt_prefix}?????? 2> /dev/null | sort`
pltlist="$pltlist5 $pltlist6"
if [ "$pltlist" ]; then
   nl=$(echo "$pltlist" | wc -l)
   nl=$(expr $nl - 1)
   if [ $nl -eq 0 ]; then
      pltlist=""
   else
      pltlist=$(echo "$pltlist" | head -$nl)
   fi
fi

for dir in ${pltlist}
do
   if [ -d ${dir} ]; then

      dirbase=`basename ${dir}`

      # only work on the file if there is not a .processed file in the
      # main directory or the plotfiles/ directory
      if [ ! -f ${dirbase}.processed ] && [ ! -f plotfiles/${dirbase}.processed ]; then

         # do processing

         # store the file on HPSS
         ${HTAR} -H copies=2 -cvf ${HPSS_DIR}/${dirbase}.tar ${dir} > ${dirbase}.htar

         # Ordinarily, we'd check htar's exit status (0 = successful), but 
         # on some machines (like Atlas) htar doesn't return a valid exit
         # status.  Instead we'll grep for the success line at the end of 
         # htar's output (which we piped into a file) and check the output 
         # status of grep
         grep "HTAR: HTAR SUCCESSFUL" ${dirbase}.htar >> /dev/null

         # The variable $? holds the exit status of the previous command
         if [ $? -eq 0 ]; then
 
            # mark this file as processed so we skip it next time
            date > ${dirbase}.processed

            # output the plotfile name and simulation time to ftime.out
            if [ -f ${FTIME_EXE} ] ; then
               ${FTIME_EXE} ${dir} >> ftime.out
            fi

            # remove the htar temporary file
            rm ${dirbase}.htar

            # move the plotfile into the plotfiles directory
            mv ${dir} plotfiles/

            # ..and the corresponding .processed file too.
            mv ${dirbase}.processed plotfiles/

            echo "$dirbase archived to HPSS"

         fi
      fi   # end test of whether plotfile already processed
   fi   # end test of whether plotfile is a directory (as it should be)
done

#==============================================================================
# Archive checkpoint files ----------------------------------------------------

# Take all but the final chk file -- we want to ensure they're completely 
# written to disk.  Strip out any tar files that are lying around as well 
# as chkXXXXX.processed files.
chklist5=`ls -d ${chk_prefix}????? 2> /dev/null | sort`
chklist6=`ls -d ${chk_prefix}?????? 2> /dev/null | sort`
chklist="$chklist5 $chklist6"
if [ "$chklist" ]; then
   nl=$(echo "$chklist" | wc -l)
   nl=$(expr $nl - 1)
   if [ $nl -eq 0 ]; then
      chklist=""
   else
      chklist=$(echo "$chklist" | head -$nl)
   fi
fi

for dir in ${chklist}
do
   if [ -d ${dir} ]; then

      dirbase=`basename ${dir}`

      if [ ! -f ${dirbase}.processed ] && [ ! -f checkfiles/${dirbase}.processed ]; then

         # store the file on HPSS
         ${HTAR} -H copies=2 -cvf ${HPSS_DIR}/${dirbase}.tar ${dir} > ${dirbase}.htar

         # Ordinarily, we'd check htar's exit status (0 = successful), but 
         # on some machines (like Atlas) htar doesn't return a valid exit
         # status.  Instead we'll grep for the success line at the end of 
         # htar's output (which we piped into a file) and check the output 
         # status of grep
         grep "HTAR: HTAR SUCCESSFUL" ${dirbase}.htar >> /dev/null

         # The variable $? holds the exit status of the previous command
         if [ $? -eq 0 ]; then

            # mark this file as processed so we skip it next time
            date > ${dirbase}.processed

            # remove the htar temporary file
            rm ${dirbase}.htar

            # move the checkpoint file into the checkfiles directory
            mv ${dir} checkfiles/

            # ..and the corresponding .processed file too.
            mv ${dirbase}.processed checkfiles/

            echo "$dirbase archived to HPSS"

         fi
      fi
   fi
done

