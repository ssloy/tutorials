import matplotlib.pyplot as plt
import numpy as np

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
print("K:",K)

X0=np.matrix([[K[0,0]],[K[0,1]]])
U=K*X0
X=Abar*X0 + Bbar*U

plt.xlabel('x (m)')
plt.ylabel('v (m/s)')
plt.plot(K[:,0],K[:,1],'o-')
plt.plot(X[0::2],X[1::2],'o-')
plt.grid(True)
plt.xlim(-1, 1)
plt.ylim(-1, 1)
plt.gca().set_aspect('equal', adjustable='box')

plt.show()
