from numpy import sin, cos
from numpy.linalg import inv
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation

g = 9.8   # free fall acceleration, m/s^2
l =  .4   # distance from the pivot to the center of mass of the pendulum, m
M =  .800 # mass of cart, kg
m =  .350 # mass of pendulum, kg
I =  .054 # moment of inertia, kg m^2

def differentiate(state, t):
    start_freq, end_freq, duration = 2., 15., 10.
    f = 30*cos(2*3.14159*((start_freq + ((end_freq-start_freq)*t)/(2.*duration))*t))
    A = np.matrix([[M+m, l*m*cos(state[2])],[l*m*cos(state[2]), I+l*l*m]])
    B = np.matrix([[f+l*m*sin(state[2])*(state[3]**2)],[g*l*m*sin(state[2])]])
    C = inv(A)*B
    return [state[1], C[0,0], state[3], C[1,0]]

dt = 0.01
t = np.arange(0.0, 10, dt)

# initial state [x v theta w]
state = [0, 0, 3.14159, 0]

y = integrate.odeint(differentiate, state, t)

fig = plt.figure()
ax = fig.add_subplot(111, autoscale_on=False, xlim=(-1, 1), ylim=(-1, 1))
ax.grid()

line, = ax.plot([], [], 'o-', lw=2)
time_template = 'time = %.1fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text


def animate(i):
    thisx = [y[i,0], l*sin(y[i,2]) + y[i,0]]
    thisy = [0,      l*cos(y[i,2])         ]

    line.set_data(thisx, thisy)
    time_text.set_text(time_template % (i*dt))
    return line, time_text

ani = animation.FuncAnimation(fig, animate, np.arange(1, len(y)), interval=1, blit=True, init_func=init)
#ani.save('double_pendulum.mp4', fps=15)
plt.show()

