#!/usr/bin/gnuplot -persist
#set key box
set key opaque

#set yrange [-1.5:1.5]
set xrange [0.001:5]
set logscale x

set xlabel "Energy(GeV)"
set ylabel "Flux"

plot "oscillations.dat" u 1:2 w l  title "nu_e" , "oscillations.dat" u 1:3 w l ls 7 title "nu_mu" ,  "oscillations.dat" u 1:4 w l title "nu_tau"

set terminal postscript eps enhanced color
set output "VacOscillation.eps"
plot "oscillations.dat" u 1:2 w l  title "{/Symbol n}_e" , "oscillations.dat" u 1:3 w l ls 7 title "{/Symbol n}_{/Symbol m}" ,  "oscillations.dat" u 1:4 w l title "{/Symbol n}_{/Symbol t}"



