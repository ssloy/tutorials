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

def gaussian(x, mu, sig):
	return np.exp(-((x-mu)/sig)**2)/(sig*np.sqrt(2*np.pi))

A = np.matrix([[1,1],[0,1]])
B = np.matrix([[0],[1]])
(n,m) = B.shape
N=60

Abar = build_Abar(A,    N)
Bbar = build_Bbar(A, B, N)

K=-(Bbar.transpose()*Bbar+np.identity(N*m)*256).I*Bbar.transpose()*Abar
X0=np.matrix([[3.1],[0.5]])
U=K*X0
X=Abar*X0 + Bbar*U

Xorig = X.copy()

X = np.matrix(np.around(X*5)/5.)


plt.plot(range(N+1), X[0::2], label="measured x(t)", color='red')
plt.plot(range(N+1), X[1::2], label="measured v(t)", color='green')
#plt.plot(range(N),   U,       label="u(t)", color='blue')

Z = X.copy()
#Z[0,0] = 2
#Z[1,0] = -1
for i in range(N):
	Zi = Z[i*n:(i+1)*n,0]
	Xi = X[i*n:(i+1)*n,0]
	Z[(i+1)*n:(i+2)*n,0] = A*Zi + B*U[i*m:(i+1)*m,0] + np.matrix([[-4./3.,0],[-5./9.,0]])*(Zi-Xi)
#	Z[(i+1)*n:(i+2)*n,0] = A*Zi + B*U[i*m:(i+1)*m,0] + np.matrix([[-0.594,-0.673],[-0.079,-0.774]])*(Zi-Xi)

print("cost:",((X-Xorig).transpose()*(X-Xorig))[0,0])
print("cost:",((Z-Xorig).transpose()*(Z-Xorig))[0,0])
	
plt.plot(range(N+1), Z[0::2], label="filtered x(t)", color='orange')
plt.plot(range(N+1), Z[1::2], label="filtered v(t)", color='magenta')

plt.xlabel('time (s)')
plt.legend()
plt.show()
