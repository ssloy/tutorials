import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

U0 = 19.2

def unit_step_current(x, R, L):
    return [U0/R - U0/R*np.exp(-t*R/L) for t in x]

data = np.genfromtxt('unit_step_19.2V.csv', delimiter=',', names=['t', 'V', 'A'])

[R, L] = curve_fit(unit_step_current, data['t'], data['A'])[0]
print(R, L)

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)

ax1.set_title("Resistance/inductance fitting")
ax1.set_xlabel('Time, seconds')
ax1.set_ylabel('Current (A), tension (V)')

ax1.plot(data['t'], data['V'], color='b', label='input tension')
ax1.plot(data['t'], data['A'], color='g', label='measured current')
model=unit_step_current(data['t'], R, L)
ax1.plot(data['t'], model, color='r', label='fitted curve')
ax1.legend()

plt.show()

