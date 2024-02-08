import numpy as np 
import matplotlib.pyplot as plt
from scipy.integrate import quad
from numpy.linalg import eig, det

def D(x,y,z,v1,v2,v3):
    return 1+(v1*v2)**2+(v2*v3)**2+(v3*v1)**2-2*((1-v1**2)*v2*v3*np.cos(x)+(1-v2**2)*v3*v1*np.cos(y)+(1-v3**2)*v1*v2*np.cos(z))

J=1

# def Eig(n,b,J):
#     v1,v2,v3=np.tanh(J*b),np.tanh(J*b),1
#     E0=[]
#     for i in range(-n,n+1):
#         for j in range(-n,n+1):
#             x,y=2*np.pi*i/n,2*np.pi*j/n 
#             E0.append(D(x,y,2*np.pi-(x+y),v1,v2,v3)/(v1*v2*v3))
#             E0.append(-D(x,y,2*np.pi-(x+y),v1,v2,v3)/(v1*v2*v3))
#     return E0


def Eig(n,b,J):
    j=complex(0,1)
    E=[]
    w1,w2,w3=1/np.tanh(J*b),1/np.tanh(J*b),1/np.tanh(J*b)
    for r in range(1,n+1):
        for s in range(1,n+1):
            p,q=np.exp(2*j*np.pi*r/n),np.exp(2*j*np.pi*s/n)
            A=[[0,-1,1,0,-w1/p,0],[1,0,-1,0,0,w2/q],[-1,1,0,-w3,0,0],[0,0,w3,0,-1,1],[w1*p,0,0,1,0,-1],[0,-w2*q,0,-1,1,0]]
            E0=eig(np.array(A))[0]
            for e in E0:
                E.append(np.imag(e))

    return E



    

# def Z(b):
#     v1,v2,v3=np.tanh(J*b),np.tanh(J*b),1

#     def D1(y):
#         def D0(x):
#             return np.log(D(x,y,2*np.pi-(x+y),v1,v2,v3))
#         return quad(D0,0,2*np.pi)[0]

    
#     return quad(D1,0,2*np.pi)[0]


B=np.linspace(0.1,5,100)
# # Z0=[]
# # for b in B:
# #     Z0.append(Z(b))

# # plt.plot(B,Z0,'r')

# # n,n1=30,10
# # E,E1=[],[]
n1=30
E1,Z=[],[]

for b in B:
 #   E.append(Eig(n,b,J))
    Eig1=Eig(n1,b,J)
    E1.append(Eig1)
    
    s=0
    for e in  Eig1:
        if e>=0:
            s+=np.log(e)

  #  Z.append(s+n1**2*np.log(2)+2*n1*(n1-1)*np.log(np.sinh(b*J)))
    #Z.append((s+(2*n1**2+4*n1-4)*np.log(2)+(2*n1*(n1-1)+n1**2+4*n1-4)*np.log(np.sinh(b*J)))/np.log(2))
    Z.append(s)

#plt.plot(B,E,'r+')
#plt.plot(B,E1,'g+')
plt.plot(B,Z,'g+')
plt.plot(B,B*0,'b')
plt.show()

# n1,b=20,0.001
# E=Eig(n1,b,J)
# s=0
# for e in  E:
#     if e>=0:
#         s+=np.log(e)

# print((s+(2*n1**2+4*n1-4)*np.log(2)+(2*n1*(n1-1)+n1**2+4*n1-4)*np.log(np.sinh(b*J)))/np.log(2))

