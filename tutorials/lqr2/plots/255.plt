set terminal pngcairo size 700,524 enhanced font 'Verdana,10'
set output '255.png'

w(x) = w0*(1-exp(alpha*x))
w0 = 100
alpha = -.1
fit w(x) "255.dat" using 1:2 via w0,alpha

plot "255.dat" using 1:2 title "run 1" with lines, "" using 1:3 title "run2" with lines,  w(x) title "w(t) fit"


