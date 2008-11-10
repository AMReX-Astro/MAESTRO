set format "%3.1l x 10^{%L}"
set term post portrait enhanced color 11

set output 'dipole.ps';

set multiplot;

set size 1, 0.333;

set rmargin 4

# make a plot of the radial velocity statistics

set origin 0.0, 0.666;

set xlabel "time (s)";
set ylabel "velocity (cm/s)" 2,0

set title "average radial velocity velocity"

plot "wdconvect_radvel_diag.out" using 1:2 title "x" with lines, "wdconvect_radvel_diag.out" using 1:3 title "y" with lines, "wdconvect_radvel_diag.out" using 1:4 title "z" with lines


# make a plot of the density-weighted radial velocity statistics

set origin 0.0, 0.333;

set xlabel "time (s)";
set ylabel "velocity (cm/s)"

set title "density-weighted average radial velocity"

plot "wdconvect_radvel_diag.out" using 1:7 title "x" with lines, "wdconvect_radvel_diag.out" using 1:8 title "y" with lines, "wdconvect_radvel_diag.out" using 1:9 title "z" with lines


# make a plot of the maximum radial velocity 

set origin 0.0, 0.0;

set xlabel "time (s)";
set ylabel "velocity (cm/s)"

set title "peak radial velocity"

plot "wdconvect_radvel_diag.out" using 1:6 notitle with lines

set nomultiplot;
set term x11;
