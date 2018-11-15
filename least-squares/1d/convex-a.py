import matplotlib.pyplot as plt

x = [0.] * 32
x[0] = x[31] = 1.
x[18] = 2.

plt.plot(x, drawstyle='steps-mid')

for iter in range(1000):
    for i in range(len(x)):
        if i in [0,18,31]: continue
        x[i] = (x[i-1]+x[i+1]+11./32**2)/2.

plt.plot(x, drawstyle='steps-mid')
plt.show()

