#!/bin/sh

# clone Maestro and all associated repos and reset them to a particular
# date/time.  Then output the necessary environment variables to build
# with this source.

branch=development


#-----------------------------------------------------------------------------

date=$1

pwd=`pwd`

# MAESTRO
echo "cloning MAESTRO"
git clone ssh://git@github.com/BoxLib-Codes/MAESTRO

echo " "
echo "resetting to before ${date} on branch ${branch}"
cd MAESTRO
git checkout ${branch}
hash=`git rev-list -n 1 --before="$date" ${branch}`
git checkout ${hash}

cd ..

# BoxLib
echo " "
echo "cloning BoxLib"
git clone ssh://git@github.com/BoxLib-Codes/BoxLib

echo " "
echo "resetting to before ${date} on branch ${branch}"
cd BoxLib
git checkout ${branch}
hash=`git rev-list -n 1 --before="$date" ${branch}`
git checkout ${hash}

cd ..


# output the necessary environment changes
if [ -f exports.sh ]; then
    echo "removing old exports.sh"
    rm -f exports.sh
fi

cat >> exports.sh << EOF 
export MAESTRO_HOME="${pwd}/MAESTRO"
export BOXLIB_HOME="${pwd}/BoxLib"
EOF

echo " "
echo "source exports.sh to setup the environment for building"




