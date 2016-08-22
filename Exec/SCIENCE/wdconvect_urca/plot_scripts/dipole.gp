set term post portrait enhanced color 11
set output 'dipole.ps';

#set xtics 200
set rmargin 4

#set key reverse

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

set yrange [-8.e4:8.e4]

set title "average radial velocity velocity"

plot "wdconvect_radvel_diag.out" using 1:2 title "<v_r>_x" with lines, \
     "wdconvect_radvel_diag.out" using 1:3 title "<v_r>_y" with lines, \
     "wdconvect_radvel_diag.out" using 1:4 title "<v_r>_z" with lines, \
     "wdconvect_radvel_diag.out" using 1:(sqrt($2**2 + $3**2 + $4**2)) title "<v_r>" with lines     

#------------------------------------------------------------------------------
# make a plot of the density-weighted radial velocity statistics
#------------------------------------------------------------------------------

set origin 0.0, 0.0;

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

set yrange [-180:180]
set ytics 30

set title "dipole angles from average velocity"

# an alternate way to do phi, so it ranges from 0 to 360
# inspiration came from here: http://www.gnuplot.info/faq/faq.html#SECTION00082000000000000000
# plot "< cat wdconvect_radvel_diag.out | sed 's/#.*//' | awk 'BEGIN { pi = 3.14159} (NF > 0) { t = $1; phi = atan2 ( $3,$2 )*180/pi; if (phi < 0) phi = phi + 360; print t, phi } '" with lines

plot "wdconvect_radvel_diag.out" using 1:(atan2($3,$2)) title "phi" with lines, \
     "wdconvect_radvel_diag.out" using 1:(atan2(sqrt($2**2 + $3**2),$4)) title "theta" with lines

set yrange [*:*]


#------------------------------------------------------------------------------
# theta, phi for density-weighted average radial velocity
#------------------------------------------------------------------------------

set format y " %g"

set origin 0.0, 0.0;
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

set ytics autofreq;

set format y "%3.2l x 10^{%L}"

set size 1, 0.5;

set origin 0.0, 0.333;

set xlabel "time (s)";
set ylabel "velocity (cm/s)" 2,0

set yrange [*:*]

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