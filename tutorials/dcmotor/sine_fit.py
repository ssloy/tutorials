import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

U0 = 16.8
F0 = 61.04*2*np.pi

def sine_current(x, R, L):
    return [(-L*np.cos(F0*t) + R/F0*np.sin(F0*t) + L*np.exp(-t*R/L))*U0*F0/((L*F0)**2 + R**2) for t in x]

data = np.genfromtxt('sine_16.8V_61.04Hz.csv', delimiter=',', names=['t', 'V', 'A'])

[R, L] = curve_fit(sine_current, data['t'], data['A'], p0=[5,.001])[0]
print(R, L)

fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.set_title("Resistance/inductance fitting")
ax1.set_xlabel('Time, sec')
ax1.set_ylabel('Current (A), tension (V)')

ax1.plot(data['t'], data['V'], color='b', label='input tension')
ax1.plot(data['t'], data['A'], color='g', label='measured current')
model=sine_current(data['t'], R, L)
ax1.plot(data['t'], model, color='r', label='fitted curve')
ax1.legend()

plt.show()

