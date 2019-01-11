#!/usr/bin/python3

import numpy as np
import sys

use_target = 0
if (use_target):
    [time, refcurrent, current, x, theta, target] = np.loadtxt(sys.argv[1], delimiter=',', skiprows=1, unpack=True)
else:
    [time, refcurrent, current, x, theta] = np.loadtxt(sys.argv[1], delimiter=',', skiprows=1, unpack=True)

time = time / 1000000.
refcurrent = refcurrent / 1000.
current = current / 1000.
x = x / 1000000.
theta = theta * (3.14159265/180./1000.)


if (use_target):
    target = target / 1000000.
    print("time(s),reference current(A),current(A),cart position(m),pendulum angle(rad),target(m)")
else:
    print("time(s),reference current(A),current(A),cart position(m),pendulum angle(rad)")
for i in range(len(time)):
    if (use_target):
        print (time[i],",",refcurrent[i],",",current[i],",",x[i],",",theta[i],",",target[i],sep='')
    else:
        print (time[i],",",refcurrent[i],",",current[i],",",x[i],",",theta[i],sep='')


