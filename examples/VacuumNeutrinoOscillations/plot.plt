#!/usr/bin/gnuplot -persist
set key box
set key opaque

#set yrange [-1.5:1.5]
set xrange [0.0005:10]
set logscale x

plot "oscillations.dat" u 1:2 w l  title "nu_e" , "oscillations.dat" u 1:3 w l ls 7 title "nu_mu" ,  "oscillations.dat" u 1:4 w l title "nu_tau"



