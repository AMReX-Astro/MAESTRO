#!/bin/sh
slice-field.py $1 -res 8192 -w 1.5e8 -axis x
slice-field.py $1 -res 8192 -w 1.5e8 -axis y
slice-field.py $1 -res 8192 -w 1.5e8 -axis z
slice-vel.py $1 -res 8192 -w 1.5e8 -axis x
slice-vel.py $1 -res 8192 -w 1.5e8 -axis y
slice-vel.py $1 -res 8192 -w 1.5e8 -axis z
slice-enucdot.py $1 -res 8192 -w 1.5e8 -axis x -logmax 11
slice-enucdot.py $1 -res 8192 -w 1.5e8 -axis y -logmax 11
slice-enucdot.py $1 -res 8192 -w 1.5e8 -axis z -logmax 11
plot-edot.py $1 -minr 2.5 -maxr 1500.0 -n 100
