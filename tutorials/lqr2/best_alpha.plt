set terminal pngcairo size 700,524 enhanced font 'Verdana,10'
set output 'best_alpha.png'
plot 'best_alpha.dat' using 1:2 title "real data" with lines, '' using 1:3 title "fitted curve" with lines

