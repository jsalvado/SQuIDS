#!/usr/bin/gnuplot -persist
set key box
set key opaque
set multiplot layout 2,1
set yrange [-1.5:1.5]
set xrange [0:120]
plot "rabi.dat" u 1:2 w l title "dipole" , "rabi.dat" u 1:3 w l ls 7 title "Occupation GS" ,  "rabi.dat" u 1:4 w l title "Occupation ES"
plot "rabi_detuned.dat" u 1:2 w l title "dipole", "rabi_detuned.dat" u 1:3 w l ls 7 title "Occupation GS",  "rabi_detuned.dat" u 1:4 w l title "Occupation ES"
unset multiplot

# set terminal postscript color 
# #set terminal png large
# set output "rabi_plot.png"

# set style line 2 lt 1 lw 2
# set style line 7 lt 1 lw 2

# set key box
# set key opaque
# set multiplot layout 2,1
# set yrange [-1.5:1.5]
# set xrange [0:120]
# plot "rabi.dat" u 1:2 w l title "dipole" , "rabi.dat" u 1:3 w l ls 7 title "Ocupation GS" ,  "rabi.dat" u 1:4 w l title "Ocupation ES"
# plot "rabi_detuned.dat" u 1:2 w l title "dipole", "rabi_detuned.dat" u 1:3 w l ls 7 title "Ocupation GS",  "rabi_detuned.dat" u 1:4 w l title "Ocupation ES"
# unset multiplot

