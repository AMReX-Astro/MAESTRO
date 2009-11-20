#!/bin/bash
#
# This is a sample script to generate some movies of 2d datasets by using the
# fsnapshot routine.  There are various options pertaining to exactly what the
# user wants to do.
#


# check to see if we just want to make the movies
only_do_movies=0
convert_images=0
do_all_files=0
while [ "$1" != "" ]; do
    case $1 in 
	-m | --movies_only)
	    only_do_movies=1
	    ;;
        -c | --convert_images)
	    convert_images=1
	    ;;
	-a | --all_files)
	    do_all_files=1
	    ;;
	*)
	    break
	    ;;
	esac
    shift
done

# THESE NEED TO BE UPDATED FOR THE CURRENT MACHINE
SNAPSHOT=/home/cmalone/MAESTRO/fParallel/data_processing/fsnapshot2d.Linux.gfortran.exe

FEXTREMA=/home/cmalone/MAESTRO/fParallel/data_processing/fextrema.Linux.gfortran.exe

FTIME=/home/cmalone/MAESTRO/fParallel/data_processing/ftime.Linux.gfortran.exe

CWD=`pwd`
DEST_DIR=${CWD}/data/images

# for best results, these should be changed based on number of zones in the
# plotfile
xsize=512
ysize=1024

CONVERT_OPTIONS="-quality 100 -resize ${xsize}x${ysize} -depth 8"

MENCODER_OPTIONS="-mf w=${xsize}:h=${ysize}:fps=12:type=png -ovc lavc -lavcopts vcodec=msmpeg4:vbitrate=16000"

vars=( X\(C12\) Machnumber )
# ivars is the column of fextrema output that contains the MINIMUM value of var
# this is a bit annoying and should be fixed to be transparent to the user
ivars=( 12 52 )

# file to save the output of ftime
TIME_FILE="${CWD}/data/images/times.txt"

# let's do this
if [ "${only_do_movies}" != "1" ]; then

    echo

# get the max and min values for each variable
    echo "Calculating max/min values..."
    for index in ${!ivars[*]}; do
	imin=${ivars[${index}]}
	imax=$(( ${imin} + 1 ))
# this nasty awk string just takes the output of fextrema, looks at the 
# appropriate columns for the max and min of variable ${vars[index]} for a 
# given plotfile and finds the GLOBAL max and min for all the pltfiles 
# currently in the directory
	awk_output_array=(`${FEXTREMA} plt* | awk \
	    "BEGIN{max=0;min=1e30};\
             (NR>2){if(\$ ${imax}>max)max=\$ ${imax};\
                    if(\$ ${imin}<min)min=\$ ${imin}};\
             END{print min,max}"`)

	min_array[${index}]=${awk_output_array[0]}
	max_array[${index}]=${awk_output_array[1]}
    done

# make the dir where we will keep the images for the movie
# make sure we don't overwrite anything we dont want to
    if [ ! -e ${DEST_DIR} ]; then
	do_all_files=1
    fi
    if [ "${do_all_files}" = "1" ]; then
	if [ -e ${DEST_DIR} ]; then
	    echo "All data in \"${DEST_DIR}\" will be erased."

	    while :
	    do
		echo -n "Are you sure you want to proceed? [y/n]:"
		read choice

		if [ "$choice" == "n" -o "$choice" == "N" ]; then
		    echo "Aborting"
		    exit
		elif [ "$choice" == "y" -o "$choice" == "Y" ]; then
		    echo "Deleting data in ${DEST_DIR}"
		    rm -r ${DEST_DIR}
		    break
		fi
	    done

	fi

	mkdir ${DEST_DIR}

	for var in ${vars[*]}; do
	    mkdir ${DEST_DIR}/${var}
	done
    fi

# get the times for the plotfiles
    echo "Getting plotfiles' time..."
    ${FTIME} -f "ES15.4" `ls -d plt*` > ${TIME_FILE}

# loop over pltfiles and start making images
    echo "Taking snapshots..."
    for pltfile in `ls -d plt*`; do

	echo "Working ${pltfile}..."

	for index in ${!vars[*]}; do

	    # check to see if we need to process this file; if not, then skip
	    # it
	    if [ -f ${DEST_DIR}/${vars[${index}]}/${pltfile}.${vars[$index]}.ppm ] && [ "${do_all_files}" != "1" ]; then
		continue
	    fi

	    # set the max and min and any logarithmic flags for the snapshot
	    # routine
	    # currently these have to be set by hand
	    case "${vars[${index}]}" in
		X\(C12\))
		    minflag="-m 1.e-10"
		    maxflag="-M ${max_array[${index}]}"
		    lflag="-l"
		    ;;
		Machnumber)
		    minflag="-m ${min_array[${index}]}"
		    maxflag="-M ${max_array[${index}]}"
		    lflag=""
		    ;;
	    esac
	
	    SNAPSHOT_OPTIONS="-cname ${vars[${index}]} ${minflag} ${maxflag} ${lflag}"
	
	    # take a snapshot
	    ${SNAPSHOT} -p ${pltfile} ${SNAPSHOT_OPTIONS}
	    mv ${pltfile}.${vars[${index}]}.ppm ${DEST_DIR}/${vars[${index}]}/${pltfile}.${vars[${index}]}.ppm
	done
    done

    echo "Done taking snapshots..."

    # convert the .ppm images into something usable by mencoder
    # we also add some text regarding the time of the snapshot - this is a bit
    # cheesey looking currently, but gets the job done...
    for var in ${vars[*]}; do

	echo "Converting images for ${var}..."

	MOVIE_FILE=${var}.avi
	
	cd ${DEST_DIR}/${var}

	for image in `ls *.ppm`; do
	    echo "Working ${image}..."

	    plt=`basename ${image} .${var}.ppm`

	    image_time=$(awk "/${plt}/ {print \$2}" ${TIME_FILE})
	    convert ${CONVERT_OPTIONS} -font helvetica -fill yellow -pointsize 20 -draw "text 50,50 't = ${image_time}'" ${image} `basename ${image} .ppm`.png
	done

	echo "Encoding to movie..."

	mencoder "mf://*.png" ${MENCODER_OPTIONS} -o ${MOVIE_FILE}

	echo "Cleaning up..."

	mv ${MOVIE_FILE} ${CWD}
	cd ${CWD}
    done
    echo "Success!"

# we already have the snapshots we want, so just make the movie
# to be honest, I forget why I had this as an option in here...
else

    echo "Making movies..."
    for var in ${vars[*]}; do
	echo "   Working on ${var}..."
	MOVIE_FILE=${var}.avi

	cd ${DEST_DIR}/${var}

	if [ "${convert_images}" = "1" ]; then
	    for image in `ls *.ppm`; do
		echo "Working on image $image"
		
		plt=`basename ${image} .${var}.ppm`

		image_time=$(awk "/${plt}/ {print \$2}" ${TIME_FILE})

		convert ${CONVERT_OPTIONS} -font helvetica -fill yellow -pointsize 20 -draw "text 10,50 't = ${image_time}'" ${image} `basename ${image} .ppm`.png
	    done
	fi

	mencoder "mf://*.png" ${MENCODER_OPTIONS} -o ${MOVIE_FILE}

	echo "Cleaning up..."
	mv ${MOVIE_FILE} ${CWD}
	cd ${CWD}

    done
	
    echo "Success!"

fi