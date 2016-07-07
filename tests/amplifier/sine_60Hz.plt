set terminal pngcairo size 1024,768 enhanced font 'Verdana,10'
set output 'sine_60Hz.png'

set ylabel "m   |   m/s   |   V"
set xlabel "time, ms"
show xlabel

plot 'sine_60Hz.dat' using 1:2 title "pwm" with lines, '' using 1:3 title "current sense" with lines

