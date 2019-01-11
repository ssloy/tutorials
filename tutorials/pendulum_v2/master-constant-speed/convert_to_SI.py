#!/usr/bin/python3

import numpy as np
import sys

[time, refcurrent, current, x, theta] = np.loadtxt(sys.argv[1], delimiter=',', skiprows=1, unpack=True)

time = time / 1000000.
refcurrent = refcurrent / 1000.
current = current / 1000.
x = x / 1000000.
theta = theta * (3.14159265/180./1000.)

print("time(s),reference current(A),current(A),cart position(m),pendulum angle(rad)")
for i in range(len(time)):
    print (time[i],",",refcurrent[i],",",current[i],",",x[i],",",theta[i],sep='')


