#!/usr/bin/python

import numpy as np
from numpy import sin, cos
from numpy.linalg import inv
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import scipy.integrate as integrate

g = 9.8   # free fall acceleration, m/s^2
l =  .4   # distance from the pivot to the center of mass of the pendulum, m
M =  .800 # mass of cart, kg
m =  .350 # mass of pendulum, kg
I =  .054 # moment of inertia, kg m^2

[time, current, x, theta] = np.loadtxt('current-chirp-2A-10s-2Hz-10Hz.csv', delimiter=',', skiprows=1, unpack=True)

def blur(x):
    tmp = x[:]
    for i in range(1,len(x)-1):
        x[i] = (tmp[i-1]+2*tmp[i]+tmp[i+1])/4

def finite_difference(time, x):
    x_dot = []
    for i in range(len(x)-1):
        x_dot.append((x[i+1]-x[i])/(time[i+1]-time[i]))
    x_dot.append(x_dot[-1])
    return x_dot


def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or abs(value - array[idx-1]) < abs(value - array[idx])):
        return current[idx-1]
    else:
        return current[idx]

def differentiate(state, t):
    start_freq, end_freq, duration = 2., 15., 10.
    f = 4.6*find_nearest(time,t)

#    f = 30*cos(2*3.14159*((start_freq + ((end_freq-start_freq)*t)/(2.*duration))*t))
    A = np.matrix([[M+m, l*m*cos(state[2])],[l*m*cos(state[2]), I+l*l*m]])
    B = np.matrix([[f+l*m*sin(state[2])*(state[3]**2)],[g*l*m*sin(state[2])]])
    C = inv(A)*B
    return [state[1], C[0,0], state[3], C[1,0]]



# initial state [x v theta w]
state = [0, 0, -3.14159, 0]

y = integrate.odeint(differentiate, state, time)

#for i in range(5):
#    blur(x)
#    blur(theta)
#    blur(current)

x_dot     = finite_difference(time, x)
theta_dot = finite_difference(time, theta)

fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.set_title("Pendulum parameters fitting")
ax1.set_xlabel('Time, sec')

ax1.plot(time, x,      color='blue', label='cart position, m')
#ax1.plot(time, x_dot,  color='teal', label='cart speed, m/s')

ax1.plot(time, theta,     color='red', label='pendulum angle, rad')
#ax1.plot(time, theta_dot, color='magenta', label='pendulum speed, rad/s')

ax1.plot(time, current, color='black', label='current, A')

ax1.plot(time, y[:,0], color='green', label='synthetic position')
ax1.plot(time, y[:,2], color='green', label='synthetic angle')

ax1.legend()
plt.show()

