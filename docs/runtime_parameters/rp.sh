#!/bin/ksh

parameters=$(grep namelist probin.f90 | awk '{print $3}' | sort)

for i in ${parameters}
do
    value=$(grep -m 1 "$i =" probin.f90 | cut -d'=' -f2 | cut -d"!" -f1)

    echo "\\verb= " $i "=" " &   & " ${value} "//"
done

