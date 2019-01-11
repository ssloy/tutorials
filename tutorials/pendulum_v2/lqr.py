#!/usr/bin/python3

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

l = .39
m = .35
M = 1.2

dt = .002
g = 9.8

A = np.matrix([[1, dt, 0, 0],[0,1, -(m*g*dt)/M,0],[0,0,1,dt],[0,0,g*(m+M)*dt/(l*M),1]])
B = np.matrix([[0],[dt/M],[0],[-dt/(l*M)]])


print(A,B)

Q = np.matrix("3 0 0 0; 0 1 0 0 ; 0 0 1 0; 0 0 0 .1")
R = np.matrix(".006")



K = dlqr(A,B,Q,R)
print(K)
print("double c[] = {%f, %f, %f, %f};" % (K[0,0], K[0,1], K[0,2], K[0,3]))

nsteps = 500
time = np.linspace(0, 2, nsteps, endpoint=True)
xk = np.matrix("0 ; 0 ; .2 ; 0")

X = []
T = []
U = []

for t in time:
    uk = K*xk
    X.append(xk[0,0])
    T.append(xk[2,0])
    v = xk[1,0]
    force = uk[0,0]
#    u = (force + np.sign(v)*2.5)/8.75
    u = (force)/12
    U.append(u)
    xk = A*xk + B*uk

plt.plot(time, X, label="cart position, meters")
plt.plot(time, T, label='pendulum angle, radians')
plt.plot(time, U, label='control current, ampers')

plt.legend(loc='upper right')
plt.grid()
plt.show()

