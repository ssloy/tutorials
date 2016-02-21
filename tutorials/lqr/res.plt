set terminal pngcairo size 700,524 enhanced font 'Verdana,10'
set output 'res.png'
plot 'res.dat' using 0:1 title "error" with lines

