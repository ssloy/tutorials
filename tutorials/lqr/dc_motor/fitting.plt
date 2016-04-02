set terminal pngcairo size 700,524 enhanced font 'Verdana,10'
set output 'fitting.png'

set ylabel "cart speed, mm/sec"
set xlabel "time, ms"
show xlabel
show ylabel

plot 'fitting.dat' using 1:2 title "24V measurement" with lines, '' using 1:3 title "v[k+1] = .97 v[k] + 0.218*24" with lines, '' using 1:4 title "18V measurement" with lines, '' using 1:5 title "v[k+1] = .97 v[k] + 0.218*18" with lines, '' using 1:6 title "12V measurement" with lines, '' using 1:7 title "v[k+1] = .97 v[k] + 0.218*12" with lines, '' using 1:8 title "6V measurement" with lines, '' using 1:9 title "v[k+1] = .97 v[k] + 0.218*6" with lines

