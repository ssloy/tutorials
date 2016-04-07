set terminal pngcairo size 1024,768 enhanced font 'Verdana,10'
set output 'limits.png'

set ylabel "max cart speed, m/sec"
set xlabel "voltage"
show xlabel
show ylabel

b =  0.001
c = -0.001
w(x) = b*x + c
fit w(x) "limits.dat" using 1:2 via b,c

plot 'limits.dat' using 1:2 title "w_0(u)", w(x) title "linear regression" with lines


