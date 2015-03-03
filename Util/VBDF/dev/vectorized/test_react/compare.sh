#!/bin/bash

if [ -z "$1" ]; then
    echo "This script is used to compare plotfiles to some reference plotfiles."
    echo "Assumes 'fcompare' is in your PATH."
    echo 
    echo "Usage: `basename $0` <pltfile prefix> <directory containing references>"
    echo 
    echo "Example: To compare all the react_rprox_* pltfiles to ./ref/react_rprox_*"
    echo "   `basename $0` react_rprox ref"
    exit
fi

for file in `ls -d $1*`; do
    echo $file: `fcompare --infile1 $file --infile2 $2/$file`
done
