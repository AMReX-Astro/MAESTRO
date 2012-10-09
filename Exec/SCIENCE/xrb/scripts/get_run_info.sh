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
    echo "${0} [-w weight] <output_file_list>"
    echo
    echo " --  weight is useful for determining the total cost at a particular"
    echo "     computing center:"
    echo "   total computational cost = total cpu_hours * weight_per_cpu_hour"
    echo 
}

if [ "$#" -eq "0" ]; then
    usage
    exit 0
fi



# start of main program


# get the grid_file if we need it by parsing the options
weight=1

while [ "$1" != "" ]; do
    case $1 in
	-w )
	    shift
	    weight=$1
	    ;;
	*)
	    break
    esac

    shift
done

# get the number of cpus from the first outputfile
n_cpus=`awk '/number of processors/ {print $5}' $1`

# parse the output files
awk "BEGIN{ \
    nmatches=0; \
    total_time=0; \
    ncpus=${n_cpus}; \
    w=${weight} \
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
    print \"Total computational cost: \" total_time*ncpus*w/3600.0 \" cpu-hours\"; \
    print \"\"
}" ${@}

