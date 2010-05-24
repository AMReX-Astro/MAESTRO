set term post eps enhanced
set output 'spherical_heat2.eps'
set format "%3.1l x 10^{%L}"

set size 0.55, 0.67;

set xlabel "radius (cm)"
set ylabel "p_0 (dyne/cm^2)"

set rmargin 4

set origin 0.0, 0.0
set xrange [0:2e8]
set xtics 1e8
set ytics 4e26
plot 'base.orig' using 1:4 w l lw 3 ti "t=0";
