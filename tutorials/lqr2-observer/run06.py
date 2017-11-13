import matplotlib.pyplot as plt
import numpy as np

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

def build_Mbar(A, N):
	n = A.shape[0]
	M = np.matrix(np.zeros( ((N+1)*n,N*n) ))
	for i in range(N):
		for j in range(N-i):
			M[(N-i)*n:(N-i+1)*n,(N-i-j-1)*n:(N-i-j)*n] = A**j
	return M

A = np.matrix([[1,1],[0,1]])
B = np.matrix([[0],[1]])
(n,m) = B.shape
N=60

Abar = build_Abar(A,    N)
Bbar = build_Bbar(A, B, N)
Mbar = build_Mbar(A,    N)

K=-(Bbar.transpose()*Bbar+np.identity(N*m)*256).I*Bbar.transpose()*Abar
X0=np.matrix([[3.1],[0.5]])
U=K*X0
X=Abar*X0 + Bbar*U

Xorig = X.copy()
X = np.matrix(np.around(X*5)/5.)

plt.plot(range(N+1), X[0::2], label="measured x(t)", color='red')
plt.plot(range(N+1), X[1::2], label="measured v(t)", color='green')

L = -Mbar*(np.identity(N*n) + Mbar.transpose()*Mbar).I*Mbar.transpose()
Y = L*(Bbar*U+Abar*X0-X)
Z = Abar*X0 + Bbar*U + Y

print("cost:",((Z-Xorig).transpose()*(Z-Xorig))[0,0])

plt.plot(range(N+1), Z[0::2], label="x(t), offline estimation", color='orange')
plt.plot(range(N+1), Z[1::2], label="v(t), offline estimation", color='magenta')

plt.xlabel('time (s)')
plt.legend()
plt.show()
