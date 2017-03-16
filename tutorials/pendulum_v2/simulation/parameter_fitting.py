#!/usr/bin/python

import numpy as np
from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit
from numpy import sin, cos
from numpy.linalg import inv
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import scipy.integrate as integrate

#g = 9.8   # free fall acceleration, m/s^2
#l =  .36  # distance from the pivot to the center of mass of the pendulum, m
#M =  .800 # mass of cart, kg
#m =  .350 # mass of pendulum, kg
#I =  .001 # moment of inertia, kg m^2
#C = 5.6

#[time, refcurrent, current, x, theta] = np.loadtxt('current_0.7A.csv', delimiter=',', skiprows=1, unpack=True)
[time, refcurrent, current, x, theta] = np.loadtxt('current-chirp-2A-5s-2Hz-9Hz.csv', delimiter=',', skiprows=1, unpack=True)

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


def interpolate_current(t):
    idx = np.searchsorted(time, t, side="left")
    if idx > 0 and (idx == len(time) or abs(t-time[idx-1]) < abs(t-time[idx])):
        return current[idx-1]
    else:
        return current[idx]

def differentiate(state, t, param_dict):
    globals().update(param_dict)

    f = interpolate_current(t)
    f = a*f

    friction = 0
    fr_eps = .001

    fr_amp = 0
    if state[1]>0:
        fr_amp = -b
    else:
        fr_amp = c

    if abs(state[1])<fr_eps:
        friction = state[1]*fr_amp/fr_eps
    else:
        friction = fr_amp

    f = f + friction

    m = .35
    M = .8

    A = np.matrix([[M+m, l*m*cos(state[2])],[l*m*cos(state[2]), I+l*l*m]])
    B = np.matrix([[f+l*m*sin(state[2])*(state[3]**2)],[9.81*l*m*sin(state[2])]])
    C = inv(A)*B
    return [state[1], C[0,0], state[3], C[1,0]]


def predict(params, x, data):
    state = [0, 0, -3.14159, 0]
    y = integrate.odeint(differentiate, state, x, args=(params.valuesdict(),))
    z = np.concatenate((y[:,1],y[:,2],y[:,3])) - data
    print params.valuesdict(),"\n", reduce(lambda x,y: x+y*y,z),"\n"
    return z

for i in range(3):
#    blur(x)
    blur(theta)
#    blur(current)

time = time[:len(time)/2]
x = x[:len(x)/2]
theta = theta[:len(theta)/2]

x_dot     = finite_difference(time, x)
theta_dot = finite_difference(time, theta)



params = Parameters()
#params.add('l', value= 0.44, min=.2,max=2)
#params.add('m', value= .20, min= .2,max=1)
#params.add('M', value= 0.53, min=.5,max=1.5)
#params.add('I', value= 0.002, min=.00001)
#params.add('a', value= 5.6, min=3,max=20)
#params.add('b', value= 2.10, min=.0001,max=7)
#params.add('c', value= 1.505, min=.0001,max=7)

params.add('l', value= 0.45, min=.2,max=2)
#params.add('m', value= .35, min= .2,max=1)
#params.add('M', value= 0.8, min=.5,max=1.5)
params.add('I', value= 0.0007, min=.00001)
params.add('a', value= 8.75, min=3,max=20)
params.add('b', value= 3.17, min=.0001,max=7)
params.add('c', value= 2.37, min=.0001,max=7)

minner = Minimizer(predict, params, fcn_args=(time, np.concatenate((x_dot,theta,theta_dot))))
result = minner.minimize(method='powell')
report_fit(result)

fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.set_title("Pendulum parameters fitting")
ax1.set_xlabel('Time, sec')


ax1.plot(time, x,      color='blue', label='cart position, m')
ax1.plot(time, x_dot,  color='cyan', label='cart speed, m/s')

ax1.plot(time, theta,     color='red', label='pendulum angle, rad')
ax1.plot(time, theta_dot, color='magenta', label='pendulum speed, rad/s')

#ax1.plot(time, current, color='black', label='current, A')

t = time


state = [0, 0, -3.14159, 0]
z = params.valuesdict()
y = integrate.odeint(differentiate, state, time, args=(z,))

ax1.plot(t, y[:,0], color='green', label='synthetic position')
ax1.plot(t, y[:,1], color='green', label='synthetic speed')
ax1.plot(t, y[:,2], color='green', label='synthetic position')
ax1.plot(t, y[:,3], color='green', label='synthetic speed')




ax1.legend()
plt.show()

