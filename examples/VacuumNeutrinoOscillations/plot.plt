#!/usr/bin/env gnuplot
#set key box
if ( GPVAL_VERSION >= 4.4 && strstrt(GPVAL_TERMINALS, 'wxt') > 0 ) set terminal wxt persist
if ( GPVAL_VERSION >= 4.4 && strstrt(GPVAL_TERMINALS, 'wxt') == 0 ) print "wxt terminal not available, proceeding with default"
if ( GPVAL_VERSION < 4.4 ) print "gnuplot is too old to check for available terminals" ; print "attempting to use wxt terminal and hoping for the best" ; set terminal wxt persist
set key opaque

#set yrange [-1.5:1.5]
set logscale x

set xlabel "Energy(GeV)"
set ylabel "Flux"

plot "oscillations.dat" u 1:2 w l  title "nu_e" , "oscillations.dat" u 1:3 w l ls 7 title "nu_mu" ,  "oscillations.dat" u 1:4 w l title "nu_tau"

set terminal postscript eps enhanced color
set output "VacOscillation.eps"
plot "oscillations.dat" u 1:2 w l  title "{/Symbol n}_e" , "oscillations.dat" u 1:3 w l ls 7 title "{/Symbol n}_{/Symbol m}" ,  "oscillations.dat" u 1:4 w l title "{/Symbol n}_{/Symbol t}"



