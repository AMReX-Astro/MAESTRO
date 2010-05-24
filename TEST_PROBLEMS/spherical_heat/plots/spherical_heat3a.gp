set term post eps enhanced
set output 'spherical_heat3a.eps'
set format "%3.1l x 10^{%L}"

set size 0.55, 0.67;

set xlabel "radius (cm)"
set ylabel "temperature (K)"

set rmargin 4

set origin 0.49, 0.0
set xrange [0:3.2e7]
set yrange [5e8:9.5e8]
set xtics 1.5e7
set ytics 1e8
plot 'base.orig'        using 1:3 w l lw 3 lt 2 ti "t=0", \
     'base.new'         using 1:3 w l lw 3 lt 1 ti "t=2, 1D MAESTRO", \
     'castro_2sec_temp' using 1:2 w p pt 4 ti "t=2, 1D CASTRO", \
     'model_cc_00065'   using 2:9 w p pt 2 ti "t=2, 3D MAESTRO";
