
#set xtics 200
set rmargin 4

#set key reverse

#==============================================================================
# radial velocity component plots
#==============================================================================

set term post eps enhanced color 11 solid
set output 'wdconvect-dipole.eps';

set format y "%3.2l x 10^{%L}"

set multiplot;

set size 1, 1;

#------------------------------------------------------------------------------
# make a plot of the density-weighted radial velocity statistics
#------------------------------------------------------------------------------

#set origin 0.0, 0.333;

set xlabel "time (s)";
set ylabel "velocity (cm/s)"

set yrange [-4.e5:4.e5]

set title "density-weighted average radial velocity"

plot "wdconvect_radvel_diag.out" using 1:7 title "< {/Symbol r} v_r >_x / <{/Symbol r}>" with lines, \
     "wdconvect_radvel_diag.out" using 1:8 title "< {/Symbol r} v_r >_y / <{/Symbol r}>" with lines, \
     "wdconvect_radvel_diag.out" using 1:9 title "< {/Symbol r} v_r >_z / <{/Symbol r}>" with lines, \
     "wdconvect_radvel_diag.out" using 1:(sqrt($7**2 + $8**2 + $9**2)) title "< {/Symbol r} v_r > / <{/Symbol r}>" with lines


set yrange [*:*]

set nomultiplot;


#==============================================================================
# radial velocity direction plots
#==============================================================================

set term post eps enhanced color 11 solid
set output 'wdconvect-angles.eps';

set multiplot;

set size 1, 1;

#------------------------------------------------------------------------------
# theta, phi for density-weighted average radial velocity
#------------------------------------------------------------------------------

set format y " %g"

#set origin 0.0, 0.333;
set angle degrees;

set xlabel "time (s)";
set ylabel "angle" 0,0

set yrange [-180:180]
set ytics 30

set title "dipole angles from density-weighted average velocity"

plot "wdconvect_radvel_diag.out" using 1:(atan2($8,$7)) title "phi" with lines, \
     "wdconvect_radvel_diag.out" using 1:(atan2(sqrt($7**2 + $8**2),$9)) title "theta" with lines



set nomultiplot;

set yrange [*:*]


#==============================================================================
# make a plot of the maximum radial velocity 
#==============================================================================

set term post eps enhanced color 11 solid
set output 'wdconvect-max_velr.eps';

set ytics autofreq;

set format y "%3.2l x 10^{%L}"

set size 1, 1;

#set origin 0.0, 0.333;

set xlabel "time (s)";
set ylabel "velocity (cm/s)" 2,0

set yrange [5.e6:5.e7]

set title "peak radial velocity"

plot "wdconvect_radvel_diag.out" using 1:6 notitle with lines


#==============================================================================
# peak temperature plot
#==============================================================================

set term post eps enhanced color 11 solid
set output 'wdconvect-max_temp.eps';

set format y "%3.2l x 10^{%L}"

set size 1, 1;

#set origin 0.0, 0.333;

set xlabel "time (s)"
set ylabel "T (K)" 2,0

set yrange [6.e8:8.e8]

set title "peak temperature"

plot "wdconvect_temp_diag.out" using 1:2 notitle with lines



#------------------------------------------------------------------------------

set term x11;