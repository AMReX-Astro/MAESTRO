set term post portrait enhanced color 11 solid
set output 'base_expand.ps';

#set xtics 200
set rmargin 4

#set key reverse

#==============================================================================
# density and temperature changes
#==============================================================================

set format y "%3.2l x 10^{%L}"

set multiplot;

set size 1, 0.5;

#------------------------------------------------------------------------------
# make a plot of the change in density
#------------------------------------------------------------------------------

set origin 0.0, 0.5;

set xlabel "radius (cm)";
set ylabel "density (g / cm^3)" 2,0

set xrange [0:2.5e8]
set yrange [1.e2:1.e10]
set logscale y

set title "density change"

plot "base.orig"      using 1:2 title "initial model" w l, \
     "base.new.1.e-4" using 1:2 title "{/Symbol r}_{cutoff} = 10^{-4}" w l, \
     "base.new.1.e5"  using 1:2 title "{/Symbol r}_{cutoff} = 10^{5}" w l, \
     "base.new.1.e6"  using 1:2 title "{/Symbol r}_{cutoff} = 10^{6}" w l 


#------------------------------------------------------------------------------
# make a plot of the temperature
#------------------------------------------------------------------------------

set origin 0.0, 0.0;

set xlabel "radius (cm)";
set ylabel "temperature (K)" 2,0

set xrange [0:2.5e8]
set yrange [1.e7:5.e9]
set logscale y

set title "temperature change"

plot "base.orig"      using 1:3 title "initial model" w l, \
     "base.new.1.e-4" using 1:3 title "{/Symbol r}_{cutoff} = 10^{-4}" w l, \
     "base.new.1.e5"  using 1:3 title "{/Symbol r}_{cutoff} = 10^{5}" w l, \
     "base.new.1.e6"  using 1:3 title "{/Symbol r}_{cutoff} = 10^{6}" w l 



set yrange [*:*]

set nomultiplot;


set term x11;