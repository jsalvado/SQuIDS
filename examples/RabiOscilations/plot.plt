#!/usr/bin/env gnuplot
if ( GPVAL_VERSION >= 4.4 && strstrt(GPVAL_TERMINALS, 'wxt') > 0 ) set terminal wxt persist
if ( GPVAL_VERSION >= 4.4 && strstrt(GPVAL_TERMINALS, 'wxt') == 0 ) print "wxt terminal not available, proceeding with default"
if ( GPVAL_VERSION < 4.4 ) print "gnuplot is too old to check for available terminals" ; print "attempting to use wxt terminal and hoping for the best" ; set terminal wxt persist
set key box
set key opaque
set multiplot layout 2,1
set yrange [-1.5:1.5]
set xrange [0:120]
set xlabel "time (a.u.)"
plot "rabi.dat" u 1:2 w l title "dipole" , "rabi.dat" u 1:3 w l ls 7 title "Occupation GS" ,  "rabi.dat" u 1:4 w l title "Occupation ES"
plot "rabi_detuned.dat" u 1:2 w l title "dipole", "rabi_detuned.dat" u 1:3 w l ls 7 title "Occupation GS",  "rabi_detuned.dat" u 1:4 w l title "Occupation ES"
unset multiplot

set terminal postscript enhanced eps color
set output "rabi.eps"


set key box
set key opaque
set multiplot layout 2,1
set yrange [-1.5:1.5]
set xrange [0:120]
set xlabel "time"
plot "rabi.dat" u 1:2 w l title "dipole" , "rabi.dat" u 1:3 w l ls 7 title "Occupation GS" ,  "rabi.dat" u 1:4 w l title "Occupation ES"
plot "rabi_detuned.dat" u 1:2 w l title "dipole", "rabi_detuned.dat" u 1:3 w l ls 7 title "Occupation GS",  "rabi_detuned.dat" u 1:4 w l title "Occupation ES"
unset multiplot


