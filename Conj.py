import numpy as np
from numpy.linalg import eig, det
import matplotlib.pyplot as plt

def Sym_Asym(n):
    S,A=[],[]
    for i in range(n):
        A1,A2=[],[]
        for j in range(n):
            if j<=i:
                A1.append(0)
                A2.append(0)
            else:
                a=np.random.randint(0,2)
                A1.append(a)
                if a==0:
                   A2.append(0)
                else: 
                    b=np.random.randint(0,2)
                    A2.append((2*b-1))
        S.append(A1)
        A.append(A2)

    return S+np.transpose(S),A-np.transpose(A)

n=10
I=np.arange(0,100,1)
L1,L2=[],[]

for i in I: 
    E=Sym_Asym(n)
    e1=eig(np.array(E[0]))
    e2=eig(np.array(E[1]))
    E0=e1[0]
    E1=np.imag(e2[0])
    e1max=np.sort(E0)
    e2max=np.sort(E1)
    L1.append(e1max[-1])
    L2.append(e2max[-1])

plt.plot(I,L1,'ro')
plt.plot(I,L2,'go')

plt.plot(I,I*0,'b')
plt.show()

