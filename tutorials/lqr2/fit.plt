set terminal pngcairo size 700,524 enhanced font 'Verdana,10'
set output 'fit.png'

w(x) = w0*(1-exp(alpha*x))
w0 = 100
alpha = -.1
fit w(x) "fit.dat" using 1:2 via w0,alpha

plot "fit.dat" using 1:2 title "" with lines, w(x) title "Best-Fit Curve"



