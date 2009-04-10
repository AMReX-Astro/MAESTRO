#!/bin/ksh

parameters=$(grep namelist probin.f90 | awk '{print $3}' | sort)

for i in ${parameters}
do
    value=$(grep $i probin.f90 | grep = | grep -v "len=" | cut -d'=' -f2)

    echo "\\verb= " $i "=" " &   & " ${value}
done

