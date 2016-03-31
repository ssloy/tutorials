set terminal pngcairo size 700,524 enhanced font 'Verdana,10'
set output 'alpha.png'
plot 'alpha.dat' using 1:2 title "energy" with lines
#plot 'alpha.dat' using 0:1 title "real data" with lines, '' using 0:2 title "fitted curve" with lines

