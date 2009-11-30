#!/bin/bash
#
# This script renames directories to be in correct numerical order if we've 
# gone over 100k steps.  E.g. doing an 'ls -d plt*' may yield:
#
#  plt00000
#  plt05000
#  plt10000
#  plt100000
#  plt15000
#  plt150000
#
# After applying this script, the output would look like:
#
# plt000000
# plt005000
# plt010000
# plt015000
# plt100000
# plt150000
#
# This is useful for scripts that utilize the numbering of the pltfiles.  To
# work on files which have prefix different from 'plt', set the d option: e.g.
# '-d chk' for chkfiles.
#


prefix="plt"

while getopts "d:" options; do
    case $options in
	d ) prefix=$OPTARG;;
    esac
done

ls -d ${prefix}* &> /dev/null
if [ $? -ne 0 ]; then
    echo 
    echo "I found no files with prefix ${prefix}"
    echo
    exit 2
fi

# find the largest filename
max_length=0
for file in `ls -d ${prefix}*`; do
    name_length=${#file}

    if [ ${name_length} -gt ${max_length} ]; then
	max_length=${name_length}
    fi
done

echo "Maximum filename length is: ${max_length}"
echo "Padding files with shorter names..."

# fix the shorter names by inserting 0's
for file in `ls -d ${prefix}*`; do
    name_length=${#file}

    if [ ${name_length} -eq ${max_length} ]; then
	continue
    else
	nzeros=$(( ${max_length} - ${name_length} ))

	zero_str=""

	for i in $(seq 1 1 ${nzeros}); do
	    zero_str=${zero_str}"0"
	done

	outfile=${prefix}${zero_str}${file#"${prefix}"}

	mv ${file} ${outfile}
    fi
done