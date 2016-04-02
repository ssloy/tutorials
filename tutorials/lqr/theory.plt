set terminal pngcairo size 700,524 enhanced font 'Verdana,10'
set output 'theory.png'

set ylabel "mm   |   mm/sec   |   V"
set xlabel "time, ms"
show xlabel
show ylabel

plot 'theory.dat' using 1:2 title "x(t), mm" with lines, '' using 1:3 title "v(t), mm/sec" with lines, '' using 1:4 title "u(t), V" with lines

