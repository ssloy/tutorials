#!/usr/bin/python3

import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt

l = 0.3882958276132531
m = 0.350
M = 1.2015786486289466
I = 0.006752705243320856
Ki = 13.125615548624943
fricpos =  4.6934679845449807
fricneg = -3.3514493255262714

def derivative(state, amps):
    f = Ki*amps
    friction = 0
    if abs(state[1])<.001:
        if f>0:
            friction = -fricpos
        else:
            friction = -fricneg

        if f*friction<0 and abs(f)<abs(friction):
            friction = 0
    else:
         if state[1]>0:
            friction = -fricpos
         else:
            friction = -fricneg
    f = f + friction

    A = np.matrix([[M+m, l*m*np.cos(state[2])],[l*m*np.cos(state[2]), I+l*l*m]])
    B = np.matrix([[f+l*m*np.sin(state[2])*(state[3]**2)],[9.81*l*m*np.sin(state[2])]])
    C = inv(A)*B
    return np.array([state[1], C[0,0], state[3], C[1,0]])

def rk4_step(yi, i, time, current): 
    h = time[i+1]-time[i]
    k1 = derivative(yi, current[i])*h
    k2 = derivative(yi+k1/2., (current[i]+current[i+1])/2.)*h
    k3 = derivative(yi+k2/2., (current[i]+current[i+1])/2.)*h
    k4 = derivative(yi+k3, current[i+1])*h
    return yi+k1/6.+k2/3.+k3/3.+k4/6.

def predict(time, current):
    sol = []
    yi = [0, 0, -3.14159, 0] # x xdot theta thetadot
    for i in range(len(time)-1):
        yi = rk4_step(yi, i, time, current)
        sol.append(yi)
    sol.append(yi)
    return np.transpose(sol)





def finite_difference(time, x):
    x_dot = []
    for i in range(len(x)-1):
        x_dot.append((x[i+1]-x[i])/(time[i+1]-time[i]))
    x_dot.append(x_dot[-1])
    return x_dot

[time, refcurrent, current, x, theta] = np.loadtxt("copley-amplifier-chirp-2.5A-5s-2Hz-9Hz_b_si.csv", delimiter=',', skiprows=1, unpack=True)

fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.set_title("Pendulum dynamics compared")
ax1.set_xlabel('Time, sec')

synthetic = predict(time, current)
ax1.plot(time, current,  label='current, amp')
ax1.plot(time, finite_difference(time, x),  label='cart speed, m/s')
ax1.plot(time, theta,   label='pendulum angle, rad')
ax1.plot(time, synthetic[1,:], label='synthetic cart speed, m/s')
ax1.plot(time, synthetic[2,:], label='synthetic pendulum angle, rad')

ax1.legend()
plt.show()

