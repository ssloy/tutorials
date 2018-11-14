import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

fig, ax = plt.subplots()

n = 100
f = [0.] * n 
lock = [0, n*57//100, n-1]
lval = [1., 2., 1.]
for i in range(len(lock)):
    f[lock[i]] = lval[i]

A = np.matrix(np.zeros((n, n)))
b = np.matrix(np.zeros((n, 1)))
"""
for i in range(n):
    if i in lock:
        b[i, 0] = lval[lock.index(i)]
        A[i, i] = 1.
    elif i<n-1 and i>0:
        A[i, i-1] = -.5
        A[i, i]   =  1.
        A[i, i+1] = -.5
"""        

for i in range(n):
    if i in lock:
        A[i, i] = 100.
        b[i, 0] = lval[lock.index(i)]*100.
    elif i<n-1 and i>0:
        A[i, i-1] = -.5
        A[i, i]   =  1.
        A[i, i+1] = -.5

print(A)

print(A.transpose()*A)
f = ((np.linalg.inv(A.transpose()*A)*A.transpose()*b).transpose().tolist()[0]) 

ax.plot(f, drawstyle='steps-mid')
plt.show()

