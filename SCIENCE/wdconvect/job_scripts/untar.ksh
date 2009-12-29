#!/bin/ksh -p

plt_prefix=wd_384_6.25e8K_1.5rotate_plt

pltlist=$(find . -maxdepth 1 -type f -name "${plt_prefix}00???.tar" -print | sort)

for tarfile in ${pltlist}
do
  if [ -f ${tarfile} ]; then

      tar xf ${tarfile}
      
      # $? holds the exit status of the previous command
      if [ $? -eq 0 ]; then
	        
	  rm ${tarfile}

      fi

  fi   # end test of whether tarfile is a file (as it should be)

done

