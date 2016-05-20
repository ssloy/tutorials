set terminal pngcairo size 700,524 enhanced font 'Verdana,10'
set output 'real.png'

set ylabel "position, cm"
set xlabel "time, s"
show xlabel
show ylabel

plot 'real.dat' using 1:2 title "theoretic x(t)" with lines, '' using 1:3 title "real x(t) without friction compensatio" with lines, '' using 1:4 title "real x(t) with friction compensation" with lines

