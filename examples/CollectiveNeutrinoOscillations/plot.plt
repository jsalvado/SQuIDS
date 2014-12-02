#!/usr/bin/env gnuplot 
#set key box
set terminal wxt persist
set key opaque

#set yrange [-1.5:1.5]
set xrange [-2:2]


set xlabel "w"
set ylabel "Swap factor"

plot "collective.dat" u 1:($2/$3) w l  title "" , "collective.dat" u 1:($3) w l ls 7 title "nu_mu" 

set terminal postscript eps enhanced color
set output "Collective.eps"
plot "collective.dat" u 1:($2/$3) w l  title "" , "collective.dat" u 1:($3) w l ls 7 title "nu_mu" 




