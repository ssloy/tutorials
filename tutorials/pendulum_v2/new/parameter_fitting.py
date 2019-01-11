#!/usr/bin/python3

import sys
import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt
from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit
from functools import reduce

#g = 9.8   # free fall acceleration, m/s^2
#l =  .36  # distance from the pivot to the center of mass of the pendulum, m
#M =  .800 # mass of cart, kg
#m =  .350 # mass of pendulum, kg
#I =  .001 # moment of inertia, kg m^2
#C = 5.6

start_angle = -3.14159

def finite_difference(time, x):
    x_dot = []
    for i in range(len(x)-1):
        x_dot.append((x[i+1]-x[i])/(time[i+1]-time[i]))
    x_dot.append(x_dot[-1])
    return x_dot

def derivative(state, amps):
    f = a*amps

    friction = 0

    if abs(state[1])<.001:
        if f>0:
            friction = -b
        else:
            friction = c

        if f*friction<0 and abs(f)<abs(friction):
            friction = 0
    else:
         if state[1]>0:
            friction = -b
         else:
            friction = c

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

def predict(params, time, current):
    globals().update(params.valuesdict())
    sol = []
    yi = [0, 0, start_angle, 0]
    for i in range(len(time)-1):
        yi = rk4_step(yi, i, time, current)
        sol.append(yi)
    sol.append(yi)
    return np.transpose(sol)

def energy(params, data):
    E = 0.
    for run in data:
        y = predict(params, run[0], run[1])
        for i in range(len(run[2])):
            E = E + (y[1,i]-run[2][i])**2 + (y[2,i]-run[3][i])**2
    print(params.valuesdict(),"\n", E)
    return E

params = Parameters()
#l = 0.4019118459932508 #.410
#m = 0.341 #.350
#M = 1.0982773610765633 #.800
#I = 0.0061883985651790879
#a= 12.39922476235734
#b=  4.5145764064555642
#c= 4.0050559941348025

l = 0.3882958276132531
m = 0.350
M =  1.2015786486289466
I = 0.006752705243320856
a= 13.125615548624943
b= 4.6934679845449807
c= 3.3514493255262714

#l = 0.210
#m = 0.350
#M = 0.468
#I = 0.017
#a = 6.865
#b = 2.436
#c = 1.676


#OrderedDict([('l', 0.3882958276132531), ('M', 1.2015786486289466), ('I', 0.006752705243320856), ('a', 13.125615548624943), ('b', 4.6934679845449807), ('c', 3.3514493255262714)]) 
# 22.4647572352



#params.add('l', value=l, min=.2,max=2)
#params.add('m', value=m, min= .1,max=2)
params.add('M', value=M, min=.1,max=5)
params.add('I', value=I, min=0,max=.1)
params.add('a', value=a, min=1,max=50)
params.add('b', value=b, min=.0001,max=7)
params.add('c', value=c, min=.0001,max=7)


runs = []

runs = runs + ["copley-amplifier-square-wave-2Hz-800mA_a_si.csv"]
#runs = runs + ["copley-amplifier-chirp-2.5A-5s-2Hz-9Hz_a_si.csv"]
#runs = runs + ["copley-amplifier-rand-5s_c_si.csv"]


data = []
for run in runs:
    [time, refcurrent, current, x, theta] = np.loadtxt(run, delimiter=',', skiprows=1, unpack=True)
    x_dot = finite_difference(time, x)
    data.append([time, current, x_dot, theta])


minner = Minimizer(energy, params, fcn_args=[data])
if (0):
    result = minner.minimize(method='powell')
    report_fit(result)
    params=result.params

#ax1.plot(time, x,      color='blue', label='cart position, m')
#ax1.plot(time, current,  color='cyan', label='cart speed, m/s')
#ax1.plot(time, refcurrent,  color='cyan', label='cart speed, m/s')

fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.set_title("Pendulum parameters fitting")
ax1.set_xlabel('Time, sec')

energy(params, data)
for run in data:
    start_angle = run[3][0]
    synthetic = predict(params, run[0], run[1])
    ax1.plot(run[0], run[1],  label='current')
    ax1.plot(run[0], run[2],  label='cart speed, m/s')
    ax1.plot(run[0], run[3],   label='pendulum angle, rad')
    ax1.plot(run[0], synthetic[1,:], label='synthetic speed')
    ax1.plot(run[0], synthetic[2,:], label='synthetic angle')

ax1.legend()
plt.show()

