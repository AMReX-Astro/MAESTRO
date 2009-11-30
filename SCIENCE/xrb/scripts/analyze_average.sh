#!/bin/bash
#
# This sample script takes the average of 2d plots for various variables and 
# creates Kippenhann-like diagrams showing the evolution of the average 
# quantities.
#

# get the average data using faverage.f90 and store it in dataDir
faverage=/home/cmalone/faverage.Linux.Intel.exe

# CHANGE THESE
vars=( X\(C12\) entropy tfromp )

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
	-v ${#vars[*]} ${vars[*]} \
	-d 2

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


# useful file
lastFile=`ls -l ${dataDir}/plt*.out | tail -1 | awk '{print $NF}'`

# get the max and min height from the last plotfile - this is assumed to be
# constant for all files
minPos=`awk '(NR==3) {print \$1}' ${lastFile}`
maxPos=`tail -1 ${lastFile} | awk '{print \$1}'`

# get the time from the last output file from faverage
maxTime=`head -1 ${lastFile} | awk '{print $3}'`


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
	    optionsString="set palette defined (-10 'white', -5 'green', -1 'blue'); set cbrange [-10:-1]"
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
	tfromp)
	    optionsString="set palette defined (6 'white', 7.5 'green', 9 'blue'); set cbrange [6:9]"
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
set yrange [${minPos}:${maxPos}]
splot "${outputFile}" u 1:2:(log10(\$3))
EOF
    
    # clean it up and make it web-friendly
    convert -trim -density 100 -rotate 90 ${psfile} ${pngfile}

done
