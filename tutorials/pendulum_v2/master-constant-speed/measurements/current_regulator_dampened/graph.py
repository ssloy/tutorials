#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import sys

[time, refcurrent, current, x, theta] = np.loadtxt(sys.argv[1], delimiter=',', skiprows=1, unpack=True)

fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.set_title("Pendulum run")
ax1.set_xlabel('Time, sec')

ax1.plot(time, refcurrent, color='green',   label='reference current, A')
ax1.plot(time, current,    color='red',    label='measured current, A')
ax1.plot(time, x,          color='blue',    label='cart position, m')
ax1.plot(time, theta,      color='magenta', label='pendulum angle, rad')

ax1.legend()
plt.show()

