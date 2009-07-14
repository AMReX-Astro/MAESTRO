#!/bin/bash
#
# this script figures out some run information for a Maestro test problem
# by parsing through the output files
#

function usage() {
    echo 
    echo "This routine parses the output files of a Maestro run and determines"
    echo "the total amount of cpu-hours, among other things, used to do the"
    echo "calculation."
    echo
    echo "Usage:"
    echo "${0} [-g|--grid_file <grid_file>] [-n ncpus] <output_file_list>"
    echo
    echo " --  grid_file is used to determine the number of processors used."
    echo "     By default, grid_file, is set to 'my_grid_file.'  "
    echo " --  If ncpus is specified, then this number will override what "
    echo "     would be taken from grid_file."
    echo 
    echo "     Note that specifying a grid_file instead of ncpus will only be"
    echo "     accurate for a single level run."
    echo 
}

if [ "$#" -eq "0" ]; then
    usage
    exit 0
fi



# start of main program


# get the grid_file if we need it by parsing the options
grid_file="my_grid_file"
n_cpus=-1

while [ "$1" != "" ]; do
    case $1 in
	-g | --grid_file)
	    shift
	    grid_file=$1
	    ;;
        -n )
	    shift
	    n_cpus=$1
	    ;;
	*)
	    break
    esac

    shift
done

# if we need a grid file, make sure it exists and then read ncpus
if [ ${n_cpus} -lt 0 ]; then
        
    if [ -f "${grid_file}" ]; then

        # it exists and we get the number of processors
	n_cpus=`awk '(NR==2) {print \$NF}' ${grid_file}`

    else
	echo "File ${grid_file} is not a regular file."
	exit
    fi

fi

# parse the output files
awk "BEGIN{ \
    nmatches=0; \
    total_time=0; \
    ncpus= ${n_cpus} \
}; \
/Time to advance timestep/ { \
    nmatches+=1; \
    total_time+=\$5 \
}; \
END{ \
    print \"\"
    print \"Number of cpus: \" ncpus; \
    print \"Total time: \" total_time \" seconds\"; \
    print \"Number of timesteps: \" nmatches; \
    print \"Average time per timestep: \" total_time/nmatches \" seconds\"; \
    print \"Total number of cpu-hours: \" total_time*ncpus/3600.0; \
    print \"\"
}" ${@}

