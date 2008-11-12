set term post portrait enhanced color 11
set output 'dipole.ps';



#set xtics 200
set rmargin 4

#==============================================================================
# radial velocity component plots
#==============================================================================

set format y "%3.2l x 10^{%L}"

set multiplot;

set size 1, 0.5;

#------------------------------------------------------------------------------
# make a plot of the radial velocity statistics
#------------------------------------------------------------------------------

set origin 0.0, 0.5;

set xlabel "time (s)";
set ylabel "velocity" 2,0

set title "average radial velocity velocity"

plot "wdconvect_radvel_diag.out" using 1:2 title "x" with lines, \
     "wdconvect_radvel_diag.out" using 1:3 title "y" with lines, \
     "wdconvect_radvel_diag.out" using 1:4 title "z" with lines

#------------------------------------------------------------------------------
# make a plot of the density-weighted radial velocity statistics
#------------------------------------------------------------------------------

set origin 0.0, 0.0;

set xlabel "time (s)";
set ylabel "velocity (cm/s)"

set title "density-weighted average radial velocity"

plot "wdconvect_radvel_diag.out" using 1:7 title "x" with lines, \
     "wdconvect_radvel_diag.out" using 1:8 title "y" with lines, \
     "wdconvect_radvel_diag.out" using 1:9 title "z" with lines


set nomultiplot;


#==============================================================================
# radial velocity direction plots
#==============================================================================

set multiplot;

set size 1, 0.5

#------------------------------------------------------------------------------
# theta, phi for average radial velocity
#------------------------------------------------------------------------------

set format y " %g"

set origin 0.0, 0.5;
set angle degrees;

set xlabel "time (s)";
set ylabel "angle" 0,0

set title "dipole angles from average velocity"

plot "wdconvect_radvel_diag.out" using 1:(atan2($3,$2)) title "phi" with lines, \
     "wdconvect_radvel_diag.out" using 1:(atan2(sqrt($2**2 + $3**2),$4)) title "theta" with lines


#------------------------------------------------------------------------------
# theta, phi for density-weighted average radial velocity
#------------------------------------------------------------------------------

set format y " %g"

set origin 0.0, 0.0;
set angle degrees;

set xlabel "time (s)";
set ylabel "angle" 0,0

set title "dipole angles from density-weighted average velocity"

plot "wdconvect_radvel_diag.out" using 1:(atan2($8,$7)) title "phi" with lines, \
     "wdconvect_radvel_diag.out" using 1:(atan2(sqrt($7**2 + $8**2),$9)) title "theta" with lines




set nomultiplot;


#==============================================================================
# make a plot of the maximum radial velocity 
#==============================================================================

set format y "%3.2l x 10^{%L}"

set size 1, 0.5;

set origin 0.0, 0.333;

set xlabel "time (s)";
set ylabel "velocity (cm/s)" 2,0

set title "peak radial velocity"

plot "wdconvect_radvel_diag.out" using 1:6 notitle with lines



#==============================================================================
# peak temperature plot
#==============================================================================

set format y "%3.2l x 10^{%L}"

set size 1, 0.5;

set origin 0.0, 0.333;

set xlabel "time (s)"
set ylabel "T (K)" 2,0

set title "peak temperature"

plot "wdconvect_temp_diag.out" using 1:2 notitle with lines



#------------------------------------------------------------------------------

set term x11;