#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import sys

[time, adc, pwm] = np.loadtxt(sys.argv[1], delimiter=',', skiprows=1, unpack=True)

fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.set_title("Pendulum run")
ax1.set_xlabel('Time, sec')

ax1.plot(time, adc, color='green',   label='reference current, A')
ax1.plot(time, pwm,    color='red',    label='measured current, A')

ax1.legend()
plt.show()

