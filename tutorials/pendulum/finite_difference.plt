set terminal pngcairo size 1024,768 enhanced font 'Verdana,10'
set output 'finite_difference.png'

set xlabel "time, seconds"
show xlabel

plot 'finite_difference.txt' using 1:2 title "pendulum angle, radians" with lines, '' using 1:3 title "speed estimation 1, radians/second" with lines, '' using 1:4 title "speed estimation 2, radians/second" with lines


