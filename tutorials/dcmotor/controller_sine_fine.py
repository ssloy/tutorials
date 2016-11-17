import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

I0 = 3.56
F0 = 61.04*2*np.pi

def sine_current(x, tau):
    return [(-tau*np.cos(F0*t) + 1/F0*np.sin(F0*t) + tau*np.exp(-t/tau))*I0*F0/((tau*F0)**2 + 1) for t in x]

data = np.genfromtxt('controller_sine_61.04Hz.csv', delimiter=',', names=['t', 'J', 'I'])

tau = curve_fit(sine_current, data['t'], data['I'])[0]
print(tau)

fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.set_title("Time constant fitting")
ax1.set_xlabel('Time, sec')
ax1.set_ylabel('Current (A)')

ax1.plot(data['t'], data['J'], color='b', label='input signal')
ax1.plot(data['t'], data['I'], color='g', label='measured current')
model=sine_current(data['t'], tau)
ax1.plot(data['t'], model, color='r', label='fitted curve')
ax1.legend()

plt.show()

