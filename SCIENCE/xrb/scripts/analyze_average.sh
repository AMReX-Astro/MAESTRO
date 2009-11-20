#!/bin/bash
#
# This sample script takes the average of 2d plots for various variables and 
# creates Kippenhann-like diagrams showing the evolution of the average 
# quantities.
#

# get the average data using faverage2d.f90 and store it in dataDir
faverage=/home/cmalone/MAESTRO/fParallel/data_processing/faverage2d.Linux.gfortran.exe

# CHANGE THESE
vars=( X\(C12\) X\(He4\) momentum entropy )

dataDir=./data

# if dataDir doesn't exist, we make it
[ ! -d ${dataDir} ] && mkdir ${dataDir}

# loop over ALL the pltfiles in the current directory
for pltfile in `ls -d plt*`; do

    # check to see if we have already processed this file
    # if so, skip it
    if [ -f ${dataDir}/${pltfile}.out ]; then
	continue
    fi

    # average it
    ${faverage} -p ${pltfile} \
	-o ${dataDir}/${pltfile}.out \
	-v ${#vars[*]} ${vars[*]}

    # loop over variables
    for varIndex in $(seq 0 $(( ${#vars[*]} - 1 ))); do

	# figure out which column of faverage output corresponds to the current
	# variable
	dataColumn=$(( 2*${varIndex} + 2 ))
	var=${vars[${varIndex}]}
	outputFile=${dataDir}/${var}.out

	# grab the height location and average info for current variable and
	# dump it to outputFile
	awk \
	    "/time/ {time=\$3} \
            (NR > 2) {print time, \$1, \$${dataColumn}}" \
	    ${dataDir}/${pltfile}.out >> ${outputFile}

        # insert a newline so gnuplot is happy
	echo >> ${outputFile}

    done

done

# get the time from the last (in terms of `ls`) output file from faverage2d
# this is the maximum time
maxTime=$(head -1 `ls -l ${dataDir}/plt*.out | tail -1 | awk '{print $8}'` | awk '{print $3}')

# loop over variables and do the plotting
for var in ${vars[*]}; do

    psfile=${dataDir}/${var}.ps
    pngfile=${var}.png
    outputFile=${dataDir}/${var}.out

    # set the options the color palette range
    # currently this has to be done by hand
    # also note that we are plotting the logarithm of these values
    case "${var}" in 
	X\(C12\))
	    optionsString="set palette defined (-10 'white', -7 'green', -4 'blue'); set cbrange [-10:-4]"
	    ;;
	X\(He4\))
	    optionsString="set palette defined (-6 'white', -3 'green', 0 'blue'); set cbrange [-6:0]"
	    ;;
	momentum)
	    optionsString="set palette defined (6 'white', 10 'green', 14 'blue'); set cbrange [6:14]"
	    ;;
	entropy)
	    optionsString="set palette defined (8 'white', 8.3 'green', 8.6 'blue'); set cbrange [8:8.6]"
	    ;;
    esac

    # dump this to gnuplot
    gnuplot <<EOF
set terminal post enhanced color
set pm3d map
set out "${psfile}"
unset key
${optionsString}
set xlabel "time(s)"
set ylabel "y(cm)"
set cblabel "log(${var})"
set title "${var} evolution"
set xrange [0:${maxTime}]
set yrange [0:1024]
splot "${outputFile}" u 1:2:(log10(\$3))
EOF
    
    # clean it up and make it web-friendly
    convert -trim -density 100 -rotate 90 ${psfile} ${pngfile}

done
