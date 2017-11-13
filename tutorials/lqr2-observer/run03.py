import matplotlib.pyplot as plt
import numpy as np
import random
	
np.set_printoptions(precision=3,suppress=True,linewidth=128,threshold=100000)

def build_Abar(A, N):
	n = A.shape[0]
	Abar = np.matrix( np.zeros(((N+1)*n,n)) )
	for i in range(N+1):
		Abar[i*n:(i+1)*n,:] = A**i
	return Abar

def build_Bbar(A, B, N):
	(n,m) = B.shape
	Bbar = np.matrix( np.zeros(((N+1)*n,N*m)) )
	for i in range(N):
		for j in range(N-i):
			Bbar[(N-i)*n:(N-i+1)*n,(N-i-j-1)*m:(N-i-j)*m] = (A**j)*B
	return Bbar

A = np.matrix([[1,1],[0,1]])
B = np.matrix([[0],[1]])
(n,m) = B.shape
N=60

Abar = build_Abar(A,    N)
Bbar = build_Bbar(A, B, N)

K=-(Bbar.transpose()*Bbar+np.identity(N*m)*256).I*Bbar.transpose()*Abar

for i in range(333):
	X0 = np.matrix([[random.uniform(-3,3)],[random.uniform(-1,1)]])
	X = (Bbar*K+Abar)*X0
	plt.plot(X[0::2], X[1::2],'-')

plt.xlabel('x (m)')
plt.ylabel('v (m/s)')
plt.grid(True)
plt.show()
