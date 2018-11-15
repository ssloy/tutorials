import matplotlib.pyplot as plt

f = [0.] * 32
f[0] = f[31] = 1.
f[18] = 2.

plt.plot(f, drawstyle='steps-mid')

for iter in range(1000):
    for i in range(len(f)):
        if i in [0,18,31]: continue
        f[i] = (f[i-1]+f[i+1])/2.

plt.plot(f, drawstyle='steps-mid')
plt.show()

