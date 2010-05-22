set term post eps color
set output 'plot-spherical.eps'

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
set xrange [0:1.5e7]
set yrange [2.3e9:2.7e9]
set ylabel "density";
plot 'base.orig' using 1:2 ti "t=0" with lines ls 1,\
     'base.new'  using 1:2 ti "t=2" with lines ls 2,\
     'castro_2sec_den' using 1:2 ti "compressible" with lines ls 3;

set origin 0.0, 0.33
set xrange [0:4e7]
set yrange [4e8:9.5e8]
set ylabel "temp";
plot 'base.orig' using 1:3 ti "t=0" with lines ls 1,\
     'base.new'  using 1:3 ti "t=2" with lines ls 2,\
     'castro_2sec_temp' using 1:2 ti "compressible" with lines ls 3;

set origin 0.0, 0.66
set xrange [0:1.5e7]
set yrange [1.5e27:1.8e27]
set ylabel "pres";
plot 'base.orig' using 1:4 ti "t=0" with lines ls 1,\
     'base.new'  using 1:4 ti "t=2" with lines ls 2,\
     'castro_2sec_pres' using 1:2 ti "compressible" with lines ls 3;

unset multiplot;
set term x11;
