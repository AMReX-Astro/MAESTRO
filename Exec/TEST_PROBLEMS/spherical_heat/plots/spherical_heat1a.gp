set term post eps enhanced
set output 'spherical_heat1a.eps'

set format "%3.1l x 10^{%L}"

set size 0.55, 0.67;

set xlabel "radius (cm)"
set ylabel "{/Symbol r}_0 (g/cm^3)"

set rmargin 4

set origin 0.49, 0.0
set xrange [0:1.5e7]
set xtics 7.5e6
set yrange [2.35e9:2.67e9]
set ytics 1e8
plot 'base.orig'       using 1:2 w l lw 3 lt 2 ti "t=0", \
     'base.new'        using 1:2 w l lw 3 lt 1 ti "t=2, 1D MAESTRO", \
     'castro_2sec_den' using 1:2 w p pt 4 ti "t=2, 1D CASTRO", \
     'model_cc_00065'  using 2:3 w p pt 2 ti "t=2, 3D MAESTRO";
