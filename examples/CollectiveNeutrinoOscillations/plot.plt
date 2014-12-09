#!/usr/bin/env gnuplot 
#set key box
if ( GPVAL_VERSION >= 4.4 && strstrt(GPVAL_TERMINALS, 'wxt') > 0 ) set terminal wxt persist
if ( GPVAL_VERSION >= 4.4 && strstrt(GPVAL_TERMINALS, 'wxt') == 0 ) print "wxt terminal not available, proceeding with default"
if ( GPVAL_VERSION < 4.4 ) print "gnuplot is too old to check for available terminals" ; print "attempting to use wxt terminal and hoping for the best" ; set terminal wxt persist
set key opaque

#set yrange [-1.5:1.5]
set xrange [-2:2]

set style line 1 lt 1 lc rgb "red" lw 3
set style line 7 lt 3 lc rgb "black" lw 3

set xlabel "w"
set ylabel "Swap factor"

plot "collective.dat" u 1:($2/$3) w l ls 1  title "Swap Factor" , "collective.dat" u 1:($3) w l ls 7 title "Initial Spectrum" 

set terminal postscript eps enhanced color
set output "Collective.eps"
plot "collective.dat" u 1:($2/$3) w l ls 1  title "Swap Factor" , "collective.dat" u 1:($3) w l ls 7 title "Initial Spectrum" 




