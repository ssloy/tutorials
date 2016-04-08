set terminal pngcairo size 1024,768 enhanced font 'Verdana,10'
set output 'theory.png'

set ylabel "m   |   m/s   |   V"
set xlabel "time, s"
show xlabel
show ylabel

plot 'theory.dat' using 1:2 title "x(t), m" with lines, '' using 1:3 title "v(t), m/s" with lines, '' using 1:4 title "u(t), V" with lines, '' using 1:5 title "eps" with lines

