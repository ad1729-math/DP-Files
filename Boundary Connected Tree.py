import math as m 
import numpy as np 
import matplotlib.pyplot as plt 

#2-Tree (coord no-3) with connected Boundary 
def Z(b,J,J1,n):
    v,v1=np.tanh(b*J),np.tanh(b*J1) 
    Z1,C1,CB1=1+v**2*v1,v*(1+v1),v1+v**2
    Z,C,CB=[Z1],[C1],[CB1]
    for i in range(1,n):
        p,q,r=Z1,C1,CB1
        Z1=p**2+v**2*q**2*v1
        C1=q*p*v**2+r*v1*q*v
        CB1=q**2*v**2+r**2*v1
        Z.append(Z1)
        C.append(C1)
        CB.append(CB1)

    return Z[-1],C[-1],CB[-1]

#Symmetric tree of degree k
def Zsym(b,J,n,k):
    v=np.tanh(b*J)
    Z1,C1=1/2*((1+v**2)+(1-v**2)),1/2*((1+v**2)-(1-v**2))
    Z,C=[Z1],[C1]
    for i in range(1,n):
        p,q=Z1,C1
        Z1=1/2*((p+v**2*q)**k+(p-v**2*q)**k)
        C1=1/2*((p+v**2*q)**k-(p-v**2*q)**k)
        Z.append(Z1)
        C.append(C1)
    return Z[-1],C[-1]

J,J1,n=0.01,0.01,20
B=np.linspace(1,10,200)
#plt.plot(B,np.log(Z(B,J,J1,n)[0]),'b+')
plt.plot(B,np.log(Z(B,J,J1,n)[1]),'r+')
plt.plot(B,np.log(Z(B,J,J1,n)[2]),'g+')
#plt.plot(B,np.log(np.log(Z(B,J,J1,n)[1])),'g+')
#plt.plot(B,np.log(Z(B,J,J1,n+2))/2**(n+2),'g+')
plt.show()
    