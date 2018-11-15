import numpy as np
import matplotlib.pyplot as plt

n=1000
x = [0.] * n
x[0] = x[-1] = 1.
m = n*57//100
x[m] = 2.

A = np.matrix(np.zeros((n, n)))
for i in range(1,n-2):
    A[i, i-1] = -1.
    A[i, i]   =  2.
    A[i, i+1] = -1.
A = A[1:-2,1:-2]
A[m-2,m-1] = 0
A[m-1,m-2] = 0

b = np.matrix(np.zeros((n-3, 1)))
b[0,0] = x[0]
b[m-2,0] = x[m]
b[m-1,0] = x[m]
b[-1,0] = x[-1]
for i in range(n-3):
    b[i,0] += 11./n**2

#print(A,b)
x2 = ((np.linalg.inv(A)*b).transpose().tolist()[0]) 
x2.insert(0, x[0])
x2.insert(m, x[m])
x2.append(x[-1])

plt.plot(x2, drawstyle='steps-mid')

"""
plt.plot(x, drawstyle='steps-mid')
for iter in range(1000):
    for i in range(len(x)):
        if i in [0,m,n-1]: continue
        x[i] = (x[i-1]+x[i+1]+11./n**2)/2.

plt.plot(x, drawstyle='steps-mid')
"""

plt.show()

