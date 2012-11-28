#!/bin/sh

# clone Maestro and all associated repos and reset them to a particular
# date/time.  Then output the necessary environment variables to build
# with this source.


username=mzingale
server=gamera.lbl.gov:/usr/local/gitroot

#-----------------------------------------------------------------------------

date=$1

pwd=`pwd`

# MAESTRO
echo "cloning MAESTRO"
git clone ${username}@${server}/MAESTRO

echo " "
echo "resetting to before ${date}"
cd MAESTRO
hash=`git rev-list -n 1 --before="$date" master`
git reset --hard ${hash}

cd ..

# BoxLib
echo " "
echo "cloning BoxLib"
git clone ${username}@${server}/BoxLib

echo " "
echo "resetting to before ${date}"
cd BoxLib
hash=`git rev-list -n 1 --before="$date" master`
git reset --hard ${hash}

cd ..

# AstroDev
echo " "
echo "cloning AstroDev"
git clone ${username}@${server}/AstroDev

echo " "
echo "resetting to before ${date}"
cd AstroDev
hash=`git rev-list -n 1 --before="$date" master`
git reset --hard ${hash}

cd ..

# MAESTRO_Exec
echo " "
echo "cloning MAESTRO_Exec"
git clone ${username}@${server}/MAESTRO_Exec

echo " "
echo "resetting to before ${date}"
cd MAESTRO_Exec
hash=`git rev-list -n 1 --before="$date" master`
git reset --hard ${hash}

cd ..


# output the necessary environment changes
if [ -f exports.sh ]; then
    echo "removing old exports.sh"
    rm -f exports.sh
fi

cat >> exports.sh << EOF 
export MAESTRO_HOME="${pwd}/MAESTRO"
export BOXLIB_HOME="${pwd}/BoxLib"
export ASTRODEV_DIR="${pwd}/AstroDev"
EOF

echo " "
echo "source exports.sh to setup the environment for building"




