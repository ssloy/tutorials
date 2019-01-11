#!/usr/bin/python3

import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt

l = 0.35
m = 0.30
M = 0.85
J = .0484
Ki = 9.27
fricpos =  3.98
fricneg = -2.87

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
#    f = f + friction

    A = np.matrix([[M+m, l*m*np.cos(state[2])],[l*m*np.cos(state[2]), J]])
    B = np.matrix([[f+l*m*np.sin(state[2])*(state[3]**2)],[9.81*l*m*np.sin(state[2])]])
    C = inv(A)*B
    return np.array([state[1], C[0,0], state[3], C[1,0]])

def energy(state):
    [x,xdot,theta,thetadot] = state
    vx,vy = thetadot*l*np.cos(theta),thetadot*l*np.sin(theta)
    vx = vx + xdot
#    return 9.81*l*m*np.cos(theta) + 0.5*m*(vx**2+vy**2) + .5*M*xdot**2

    return .5*J*thetadot**2 + m*9.8*l*(np.cos(theta)-1.)


xi = [0., 0., -3.14159, 0.] # x xdot theta thetadot
synthetic = [xi[:]]
E = [energy(xi)]
current = [0]
current_prev = 0.
current_new  = 0.

def sign(val):
    if val<0.:
        return -1.
    else:
        return 1.

time = np.linspace(0, 15, 15*200)
for i in range(len(time)-1):
    current_prev = current_new

    e = energy(xi)-.2
    [x,xdot,theta,thetadot] = xi
    E.append(e)


    current_new = sign(np.cos(xi[2])*xi[3]) * e * .7
    if current_new>.7:
        current_new = .7

    if current_new<-.7:
        current_new = -.7
    if (current_new>0 and x>.15) or (current_new<0 and x<-.15):
        current_new = 0


    current.append(current_new)

    h = time[i+1]-time[i]
    k1 = derivative(xi, current_prev)*h
    k2 = derivative(xi+k1/2., (current_prev+current_new)/2.)*h
    k3 = derivative(xi+k2/2., (current_prev+current_new)/2.)*h
    k4 = derivative(xi+k3, current_new)*h
    xi = xi+k1/6.+k2/3.+k3/3.+k4/6.

    synthetic.append(xi[:])
synthetic = np.transpose(synthetic)



fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.set_title("Pendulum dynamics")
ax1.set_xlabel('Time, sec')

ax1.plot(time, synthetic[0,:], label='cart position, m')
#ax1.plot(time, synthetic[1,:], label='cart speed, m/s')
ax1.plot(time, synthetic[2,:], label='pendulum angle, rad')
ax1.plot(time, E, label='energy')
ax1.plot(time, current, label='current')

ax1.legend()
plt.show()



