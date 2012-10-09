set term post eps color
set output 'plot-planar.eps'

set multiplot;

set size 0.5, 0.33;

set xlabel "x";
set xrange [0:2.5e8]
set xtics 1e8

set logscale y

set style line 1 lw 1
set style line 2 lw 1
set style line 3 lw 1

set origin 0.0, 0.0
set yrange [1e-5:1e10]
set ylabel "density";
plot 'base.orig' using 1:2 ti "t=0" with lines ls 1,\
     'base.new'  using 1:2 ti "t=5" with lines ls 2,\
     'compressible-1.e17-t=5.00_planar.out' using 1:2 ti "compressible" with lines ls 3;

set origin 0.0, 0.33
set yrange[4e6:1e10]
set ylabel "temp";
plot 'base.orig' using 1:3 ti "t=0" with lines ls 1,\
     'base.new'  using 1:3 ti "t=5" with lines ls 2,\
     'compressible-1.e17-t=5.00_planar.out' using 1:3 ti "compressible" with lines ls 3;

set origin 0.0, 0.66
set yrange[1e13:1e28]
set ylabel "pres";
plot 'base.orig' using 1:4 ti "t=0" with lines ls 1,\
     'base.new'  using 1:4 ti "t=5" with lines ls 2,\
     'compressible-1.e17-t=5.00_planar.out' using 1:4 ti "compressible" with lines ls 3;

unset multiplot;
set term x11;
