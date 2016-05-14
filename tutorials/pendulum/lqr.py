#!/usr/bin/python

import numpy as np
import scipy.linalg
import matplotlib.pyplot as plt

def dlqr(A,B,Q,R):
    """
    Solve the discrete time lqr controller.
    x[k+1] = A x[k] + B u[k]
    cost = sum x[k].T*Q*x[k] + u[k].T*R*u[k]
    """
    # first, solve the ricatti equation
    P = np.matrix(scipy.linalg.solve_discrete_are(A, B, Q, R))
    # compute the LQR gain
    K = np.matrix(scipy.linalg.inv(B.T*P*B+R)*(B.T*P*A))
    return -K

l = .2 # rod length 2l
m = (2*l)*(.006**2)*(3.14/4)*7856 # rod 6 mm diameter, 40cm length, 7856 kg/m^3
M = .4
dt = .02 # 20 ms
g = 9.8

A = np.matrix([[1, dt, 0, 0],[0,1, -(3*m*g*dt)/(7*M+4*m),0],[0,0,1,dt],[0,0,(3*g*(m+M)*dt)/(l*(7*M+4*m)),1]])
B = np.matrix([[0],[7*dt/(7*M+4*m)],[0],[-3*dt/(l*(7*M+4*m))]])


print A,B

Q = np.matrix("1 0 0 0; 0 .001 0 0 ; 0 0 1 0; 0 0 0 .001")
R = np.matrix(".01")

K = dlqr(A,B,Q,R)
print K

nsteps = 1000
time = np.linspace(0, 2, nsteps, endpoint=True)
xk = np.matrix(".1 ; -.1 ; .2 ; .2")

X = []
V = []
T = []
W = []
F = []

for t in time:
    uk = K*xk
    X.append(xk[0,0])
    V.append(xk[1,0])
    T.append(xk[2,0])
    W.append(xk[3,0])
    F.append(uk[0,0])
    xk = A*xk + B*uk

plt.plot(time, X,label="x")
plt.plot(time, V,label='v')
plt.plot(time, T,label='t')
plt.plot(time, W,label='w')
plt.legend(loc='upper right')
#plt.plot(time, F)
plt.show()

