import numpy as np 
import matplotlib.pyplot as plt 
from numpy.linalg import eig

def Sier(b,J,n): 
    v=np.tanh(b*J)
    w=1/v
    A0=[[0,-w,0,0,0,0],[w,0,1,0,-1,0],[0,-1,0,w,1,0],[0,0,-w,0,0,0],[0,1,-1,0,0,w],[0,0,0,0,-w,0]]
    A1=A0
    for l in range(n):
        m=len(A1)
        A2,B=[],[]

        for i in range(m):
            B0=[]
            for j in range(m):
                if i==0 and j==m-1: 
                    B0.append(w)
                elif i==m-1 and j==0:
                    B0.append(-w)
                else:
                    B0.append(0)
            B.append(B0)

        for i in range(3*m):
            Am=[]
            for j in range(3*m):
                if i<m: 
                    if j<m:
                       Am.append(A1[i][j])
                    elif m<=j<2*m:
                        Am.append(B[i][j-m])
                    else:
                        Am.append(B[i][j-2*m])
                elif m<=i<2*m: 
                    if j<m:
                       Am.append(B[i-m][j])
                    elif m<=j<2*m:
                        Am.append(A1[i-m][j-m])
                    else:
                        Am.append(B[i-m][j-2*m])
                else:
                    if j<m:
                       Am.append(B[i-2*m][j])
                    elif m<=j<2*m:
                        Am.append(B[i-2*m][j-m])
                    else:
                        Am.append(A1[i-2*m][j-2*m])

            A2.append(Am)

        A1=A2 

    return A1
    
n,J=5,1
B=np.linspace(0.01,5,100)
E,Z=[],[]
for b in B:
    e0=eig(np.array(Sier(b,J,n)))[0]
    e0c=np.imag(e0)
    E.append(e0c)
    
    Z0=(3**(n+1))*np.log(np.sinh(b*J))+(4*3**n-3*(3**n-1)/2)*np.log(2)
    s0=0
    for e in e0c:
        if e==0:
           s0=0
           In=0
           break 
        elif e>0:
           In=1
           s0+=np.log(e)

    Z.append((s0+Z0)*In)
    

plt.plot(B,E,'r+')
#plt.plot(B,Z,'g+')
plt.plot(B,B*0,'b')
plt.show()



