set terminal pngcairo size 700,524 enhanced font 'Verdana,10'
set output 'z.png'
plot 'z.dat' using 1:2 title "real data" with lines, '' using 1:3 title "fitted curve" with lines, '' using 1:4 title "real data" with lines, '' using 1:5 title "fitted curve" with lines, '' using 1:6 title "real data" with lines, '' using 1:7 title "fitted curve" with lines, '' using 1:8 title "real data" with lines, '' using 1:9 title "fitted curve" with lines

