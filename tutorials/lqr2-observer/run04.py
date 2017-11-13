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

slope = gaussian(np.linspace(-5, 5, 1000), 2, .4)/6

'''
b = [slope[0]]
for item in slope[1:]:
	b.append(b[-1] - item)

plt.ylabel('elevation')
plt.xlabel('x')
plt.plot(np.linspace(-5, 5, 1000), b)
plt.show()
'''

for i in range(N):
	Xi = X[i*n:(i+1)*n,0]
	U[i*m:(i+1)*m,0] = K[0:m,:]*Xi
	
	idxslope = int((Xi[0,0]+5)/10.*1000.)
	if (idxslope<0 or idxslope>=1000):
		idxslope = 0
	X[(i+1)*n:(i+2)*n,0] = A*Xi + B*U[i*m:(i+1)*m,0] + B*slope[idxslope]

#plt.plot(range(N+1)[0:N//3], X[0::2][0:N//3], label="x(t)", color='red')
plt.plot(range(N+1), X[0::2], label="x(t)", color='red')
plt.plot(range(N+1), X[1::2], label="v(t)", color='green')
plt.plot(range(N),   U,       label="u(t)", color='blue')
plt.xlabel('time (s)')
plt.legend()
plt.show()
