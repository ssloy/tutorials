#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import sys

def finite_difference(time, x):
    x_dot = []
    for i in range(len(x)-1):
        x_dot.append((x[i+1]-x[i])/(time[i+1]-time[i]))
    x_dot.append(x_dot[-1])
    return x_dot



[time, refcurrent, current, x, theta] = np.loadtxt(sys.argv[1], delimiter=',', skiprows=1, unpack=True)
[time2, refcurrent2, current2, x2, theta2] = np.loadtxt(sys.argv[2], delimiter=',', skiprows=1, unpack=True)

x_dot = finite_difference(time, x)
#x_dot.append(0)
x_dot2 = finite_difference(time2, x2)
#x_dot2.append(0)



fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.set_title("Pendulum run")
ax1.set_xlabel('Time, sec')

#ax1.plot(time, refcurrent, color='green',   label='reference current, A')
ax1.plot(time, current,    color='red',    label='measured current, A')
ax1.plot(time, x,          color='blue',    label='cart position, m')
ax1.plot(time, theta,      color='magenta', label='pendulum angle, rad')

ax1.plot(time2, refcurrent2, color='green',   label='reference current, A')
ax1.plot(time2, current2,    color='blue',    label='measured current, A')
ax1.plot(time2, x2,          color='blue',    label='cart position, m')
ax1.plot(time2, theta2,      color='magenta', label='pendulum angle, rad')
ax1.plot(time2, x_dot2,      color='magenta', label='pendulum angle, rad')
ax1.plot(time, x_dot,      color='magenta', label='pendulum angle, rad')

ax1.legend()
plt.show()

