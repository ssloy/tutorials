import matplotlib.pyplot as plt
import numpy as np
import random

N=60

A=np.matrix([[1,1],[0,1]])
B=np.matrix([[0],[1]])

G=np.matrix(np.zeros(((N+1)*2,N)))
for i in range(N):
    T=B
    for j in range(N-i):
        G[N*2-i*2  ,N-i-1-j] = T[0,0]
        G[N*2-i*2+1,N-i-1-j] = T[1,0]
        T = A*T

H=np.matrix(np.zeros(((N+1)*2,2)))
T=np.identity(2)
for i in range(N+1):
    for j in range(2):
      for k in range(2):
        H[i*2+j,k] = T[j,k]
    T=A*T

K=-((G.transpose()*G+np.identity(N)*256).I*G.transpose()*H)

X0=np.matrix([[3.1],[0.5]])
X=(G*K+H)*X0
plt.plot(range(N+1), X[0::2], label="x(t)", color='red')
plt.plot(range(N+1), X[1::2], label="v(t)", color='green')
plt.plot(range(N),   K*X0,    label="u(t)", color='blue')
plt.legend()
plt.show()


#print("K:",K)
#print("X0:",X0)


#print("U:",K*X0)

#
#print(K[0,:]*X[4:6,0])

X0=np.matrix([[3.1],[0.5]])
for i in range(120):
	X=(G*K+H)*X0
	plt.plot(X[0::2],X[1::2],'o-')
	break
	X0=np.matrix([[random.uniform(-3,3)],[random.uniform(-1,1)]])


plt.plot(K[:,0], K[:,1], 'o-')

#X0=K[0,:].transpose()
#X=(G*K+H)*X0
#plt.plot(X[0::2],X[1::2],'o-')

#K=K*10
#plt.plot(-K[:,0],-K[:,1],'-')
#plt.plot( K[:,1], K[:,0],'-')
#plt.plot(-K[:,1],-K[:,0],'-')
#plt.show()


K0=K[0,:]
X0=np.matrix([[3.1],[0.5]])
M1=np.matrix([[0.,0.],[0.,0.]])
M1[:,0]=X0
M1[:,1]=K0.transpose()
M2=np.matrix([[0.,0.],[0.,0.]])
M2[:,1]=X0
M2[:,0]=K0.transpose()
M=M1*M2.I

K=K*M.transpose()
#print(K)

plt.plot(K[:,0],K[:,1],'o-')
plt.grid(True)
plt.xlim(-1, 1)
plt.ylim(-1, 1)
plt.gca().set_aspect('equal', adjustable='box')
plt.show()
