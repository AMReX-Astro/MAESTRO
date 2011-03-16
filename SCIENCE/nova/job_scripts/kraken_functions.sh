#==============================================================================
# Get the checkpoint file prefix from the inputs file
# - arguments: name of the inputs file
# - return: basename of checkpoint files

function get_chk_prefix {

   # Verify that an argument was received
   if [ $# -eq 1 ]; then
      inputs_file=$1
      chk_prefix=`grep check_base_name ${inputs_file} | grep -v --regexp "^ *\!" | sed 's| *check_base_name *= *"*\([^"]*\)"*|\1|'`
      if [ "$chk_prefix" = "" ]; then
         chk_prefix=chk # if it's not specified in the inputs, use the default
      fi
   else
      # No inputs file to search, so don't bother searching
      chk_prefix=""
   fi

   echo $chk_prefix
   
}

#==============================================================================
# Get the plot file prefix from the inputs file
# - arguments: name of the inputs file
# - return: basename of plot files

function get_plt_prefix {

   # Verify that an argument was received
   if [ $# -eq 1 ]; then
      inputs_file=$1
      plt_prefix=`grep plot_base_name ${inputs_file} | grep -v --regexp "^ *\!" | sed 's| *plot_base_name *= *"*\([^"]*\)"*|\1|'`
      if [ "$plt_prefix" = "" ]; then
         plt_prefix=chk # if it's not specified in the inputs, use the default
      fi
   else
      # No inputs file to search, so don't bother searching
      plt_prefix=""
   fi

   echo $plt_prefix
   
}

#==============================================================================
# Search for most-recent checkpoint file
# - arguments: basename of checkpoint files
# - return: name of last usable checkpoint file

function get_last_chk {

   # Verify that an argument was received
   if [ $# -eq 1 ]; then
      chk_base=$1

      # Look for 6-digit checkpoint file
      last_chk=`ls -d ${chk_base}?????? 2> /dev/null | sort | tail -1`

      # Check that the last 6-digit checkpoint file has a header
      # - If not, assume the last file was interrupted during writing and try
      #   the second-to-last
      if [ ! -f ${last_chk}/Header ]; then

         # Verify that more than one 6-digit checkpoint file exists
         Nchk=`ls -d ${chk_base}?????? 2> /dev/null | sort | wc -l`
         if [ $Nchk -gt 1 ]; then
   
            # Assume second-to-last 6-digit checkpoint file is usable
            last_chk=`ls -d ${chk_base}?????? 2> /dev/null | sort | tail -2 | head -1`
         else
   
            # No second-to-last 6-digit checkpoint file: try 5-digit
            last_chk=`ls -d ${chk_base}????? 2> /dev/null | sort | tail -1`

            # Check the header for the last file
            if [ ! -f ${last_chk}/Header ]; then
      
               # Verify that more than one 5-digit checkpoint file exists
               Nchk=`ls -d ${chk_base}????? 2> /dev/null | sort | wc -l`
               if [ $Nchk -gt 1 ]; then

                  # Assume second-to-last 5-digit checkpoint file is usable
                  last_chk=`ls -d ${chk_base}????? 2> /dev/null | sort | tail -2 | head -1`
               else
                  last_chk=""
               fi
            fi
         fi
      fi
   else

      # No chk_base supplied, so don't bother searching
      last_chk=""
   fi

   # "return"
   echo $last_chk
}

#==============================================================================
# Get maximum number of steps from inputs file
# - arguments: name of inputs file
# - return: inputs_file.max_step

function get_max_step {
   
   # Verify that the arguments were received
   if [ $# -eq 1 ]; then
      inputs_file=$1

      # Get maximum number of steps from inputs file
      max_steps=`grep max_step ${inputs_file} | grep -v --regexp "^ *\!"`
      if [ "$max_steps" = "" ]; then
         max_steps=1 # MAESTRO's default value
      else
         max_steps=`echo $max_steps | sed 's/ *max_step *= *\([^! ]\)/\1/'`
         max_steps=`echo $max_steps | sed 's_\([^eEdD]*\)[eEdD]\(.*\)_\1 * 10 ^ (0\2)_' | bc -l`
         # if not in scientific notation, sed and bc just return the value
         # if in scientific notation, parses into a number
      fi

      echo $max_steps
   else
      # No argument supplied: return an error
      echo -1
   fi
}

#==============================================================================
# Get maximum time from inputs file
# - arguments: name of inputs file
# - return: inputs_file.max_time

function get_max_time {

   # Verify that the arguments were received
   if [ $# -eq 1 ]; then
      inputs_file=$1

      # Get time from inputs file
      max_time=`grep stop_time ${inputs_file} | grep -v --regexp "^ *\!"`
      if [ "$max_time" = "" ]; then
         max_time=-1.0
      else
         max_time=`echo $max_time | sed 's/ *stop_time *= *\([^! ]*\)/\1/' | sed 's_\([0-9\.]*\)[eEdD]\(.*\)_\1 * 10 ^ (0\2)_' | bc -l`
      fi

      # Compare
      echo $max_time
   else
      # No arguments supplied: return an error
      echo -1
   fi
}

