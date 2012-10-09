#!/bin/bash
#
# NOTE this script assumes that gnuplot, convert, and spiffify are in your path
#

gnuplot_file="make_time_evol_plots.gn"

spiffify_options="-0 1:0:0 -1 0:0:1"

convert_options="-rotate 90 -resize 800x800"

gnuplot ${gnuplot_file}
spiffify ${spiffify_options} `ls 0.5*.ps`

# fix the colors of the labels
for file in `ls 0.5*.spiffy`; do
    awk 'BEGIN{ \
y1_color="1 0 0 setrgbcolor\n"; \
y2_color="0 0 1 setrgbcolor\n"; \
title_color="1 0 1 setrgbcolor\n" \
}; \
/max/ {\
       if (/y/) {print y2_color $0} \
       else {print y1_color $0} \
}; \
/Temporal/ {print title_color $0}; \
{print}' ${file} > temp

    mv temp ${file}

    convert ${convert_options} ${file} `basename ${file} .ps.spiffy`.png

    rm ${file}

done

