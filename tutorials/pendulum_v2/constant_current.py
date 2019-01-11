import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

#I0 = 3.56
#F0 = 61.04*2*np.pi

data = np.genfromtxt('current1.2A.csv', delimiter=',', names=['t', 'X', 'J', 'I'])
def friction(x, a):
    integral = [0]
    for i in range(len(x)-1):
        integral.append(integral[-1] + (a)*(x[i+1]-x[i]))
    integral2 = [0]
    for i in range(len(x)-1):
        if (x[i]>0.4):
            integral2.append(data['X'][i])
        else:
            integral2.append(integral2[-1] + integral[i]*(x[i+1]-x[i]))
    return integral2



V = [0]
for i in range(len(data['t'])-1):
    V.append((data['X'][i+1]-data['X'][i])/(data['t'][i+1]-data['t'][i]))


print(V)

[a] = curve_fit(friction, data['t'], data['X'])[0]
print(a)

ffig = plt.figure()
aax1 = fig.add_subplot(111)

aax1.set_title("Time constant fitting")
aax1.set_xlabel('Time, sec')
aax1.set_ylabel('Current (A), cart position (m)')

ax1.plot(data['t'], data['I'], color='r', label='real current')
aax1.plot(data['t'], data['X'], color='g', label='cart position')
ax1.plot(data['t'],       V  , color='b', label='cart speed')

mmodel=friction(data['t'], a)
aax1.plot(data['t'], model, color='r', label='fitted curve')
aax1.legend()

plt.show()

