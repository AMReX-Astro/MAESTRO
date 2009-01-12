#!/bin/bash

#
# This script moves all the plotfiles in the current directory to the directory
# supplied in $args and renames the last checkfile, e.g. chk99999, to the 
# appropriate name for restart for the next 100k timesteps using MAESTRO.
#

if [ "$#" -eq 0 ]
    then
    echo ""
    echo "Usage `basename $0` [-c, -i inputsFile] dirName"
    echo "     where dirName is the name of the directory to store plt files"
    echo "     e.g. first100k"
    echo ""
    echo "     options:"
    echo "         -c : copy chk99999/ to dirName before renaming it to"
    echo "              chk00001"
    echo "         -i : modify inputsFile to have restart = 1"
    echo ""
    exit
fi

copy_chkfile=0
modify_inputs=0
inputsFile=""

RM="rm -r"
CP="cp -r"

last_chkfile=chk99999/
first_chkfile=chk00001

# parse the options
while getopts ":ci:" Option
  do
  case $Option in 
      c ) copy_chkfile=1;;
      i ) modify_inputs=1;
	  inputsFile=$OPTARG;;
  esac
done

shift $(($OPTIND - 1))

# the remaining item in $@ should be dirName since we have shifted past all
# of the options
dirName=$1

# check to see if dirName already exists; if so, warn the user and ask for 
# verification before proceeding
if [ -e $dirName ]
    then
    echo "Directory \"$dirName\" already exists!"
    echo "All data in \"$dirName\" will be erased."

    while :
    do
      echo -n "Are you sure you want to proceed? [y/n]:"
      read choice

      if [ "$choice" == "n" -o "$choice" == "N" ]
	  then
	  echo "Please run `basename $0` again with a different dirName."
	  exit
      elif [ "$choice" == "y" -o "$choice" == "Y" ]
	  then
	  echo "Deleting data in $dirName..."
	  $RM $dirName
	  break
      fi
    done
fi

# if we get here then we are go for moving the directories
echo "Moving pltfiles and renaming chk99999/"
mkdir $dirName
mv plt* ${dirName}/

if [ $copy_chkfile -eq 1 ]
    then
    $CP $last_chkfile ${dirName}/
    echo "     Copied \"$last_chkfile\" to the \"${dirName}/\" directory..."
fi

# remove the unwanted chkfiles
$RM chk[0-8]* chk9[0-8]* chk99[0-8]* chk999[0-8]* chk9999[0-8]*

# rename chk99999/ to chk00001/ and fix the filenames within that directory to 
# numbered 00001
mv $last_chkfile $first_chkfile
cd $first_chkfile
for file in *99999
  do
  mv $file `basename $file 99999`00001
done
cd ..

echo "Done."

# update the inputs file if need be
if [ $modify_inputs -eq 1 ]
    then
    echo "Changing the inputs file, \"${inputsFile},\" to restart from 00001."
    sed -i '/^[ ]*restart =/c\ restart = 1' $inputsFile
    echo "Done."
fi

exit