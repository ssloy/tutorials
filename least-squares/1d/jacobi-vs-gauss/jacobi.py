import matplotlib.pyplot as plt

f = [0.] * 32
f[0] = f[-1] = 1.
f[18] = 2.

plt.plot(f, drawstyle='steps-mid')

for iter in range(1):
    g = [0.]*32
    g[0] = f[1]
    for i in range(1, len(f)-1):
        g[i] = (f[i-1]+f[i+1])/2. 
    g[-1] = g[-2]
    f = g[:]

plt.plot(f, drawstyle='steps-mid')
plt.show()

