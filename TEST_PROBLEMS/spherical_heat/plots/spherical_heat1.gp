set term post eps enhanced
set output 'spherical_heat1.eps'

set format "%3.1l x 10^{%L}"

set size 0.55, 0.67;

set xlabel "radius (cm)"
set ylabel "{/Symbol r}_0 (g/cm^3)"

set rmargin 4

set origin 0.0, 0.0
set xrange [0:2e8]
set xtics 1e8
plot 'base.orig' using 1:2 w l lw 3 ti "t=0";

