import matplotlib.pyplot as plt

f = [0.] * 32
f[0] = f[-1] = 1.
f[18] = 2.

plt.plot(f, drawstyle='steps-mid')

for iter in range(1):
    f[0] = f[1]
    for i in range(1, len(f)-1):
        f[i] = (f[i-1]+f[i+1])/2. 
    f[-1] = f[-2]

plt.plot(f, drawstyle='steps-mid')
plt.show()

