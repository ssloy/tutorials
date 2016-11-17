import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

I0 = 4.08

def unit_step_current(x, tau):
    return [I0 - I0*np.exp(-t/tau) for t in x]

data = np.genfromtxt('controller_unit_step_4.08A.csv', delimiter=',', names=['t', 'J', 'I'])

tau = curve_fit(unit_step_current, data['t'], data['I'])[0]
print(tau)

fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.set_title("Time constant fitting")
ax1.set_xlabel('Time, sec')
ax1.set_ylabel('Current (A)')

ax1.plot(data['t'], data['J'], color='b', label='input signal')
ax1.plot(data['t'], data['I'], color='g', label='measured current')
model=unit_step_current(data['t'], tau)
ax1.plot(data['t'], model, color='r', label='fitted curve')
ax1.legend()

plt.show()

