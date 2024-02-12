import math as m 
import numpy as np 
import matplotlib.pyplot as plt 
from numpy.linalg import eig, det

def A(b,n,J,J1): 
    v,v1=np.tanh(b*J),np.tanh(b*J1)
    
    A1=[[0,-1,1],[1,0,1],[-1,-1,0]]
    O=[[0,-1,1],[1,0,1],[-1,-1,0]]

    for g in range(1,n):

        s=len(A1)
        B,B1,C=[],[],[]

        for i in range(1,4):
            Bm,Bm3=[],[]
            for j in range(1,s+1):
                if i==1 and j==int(s/2)+1:
                    Bm.append(1/v)
                    Bm3.append(0)
                elif i==3 and j==int(s/2)+1:
                    Bm.append(0)
                    Bm3.append(1/v)
                else:
                    Bm.append(0)
                    Bm3.append(0)
            B.append(Bm)
            B1.append(Bm3)

        for i in range(1,s+1):
            Cm=[]
            for j in range(1,s+1):
                if i==s and j==1:
                    Cm.append(1/v1)
                else:
                    Cm.append(0)    

            C.append(Cm)
   
        BT,B1T,CT=np.transpose(B),np.transpose(B1),np.transpose(C)

        s1=3+2*s
        A10=[]
        for i in range(1,s+1):
            H=[]
            for j in range(1,s1+1):
                if j<=s:
                   H.append(A1[i-1][j-1])
                elif s<j<=s+3:
                    H.append(-BT[i-1][j-s-4])
                else: 
                    H.append(C[i-1][j-s-4])

            A10 += [H]

        for i in range(s+1,s+4):
            H=[]
            for j in range(1,s1+1):
                if j<=s:
                   H.append(B[i-s-1][j-1])
                elif 3<j<=s+3:
                    H.append(O[i-s-1][j-s-1])
                else: 
                    H.append(B1[i-s-1][j-s-4])

            A10 += [H]

        for i in range(s+4,s1+1):
            H=[]
            for j in range(1,s1+1):
                if j<=s:
                   H.append(-CT[i-s-4][j-1])
                elif 3<j<=s+3:
                    H.append(-B1T[i-s-4][j-s-1])
                else: 
                    H.append(A1[i-s-4][j-s-4])

            A10 += [H]
    
        A1=A10

    s0=len(A1)
    # Am=[]
    # for i in range(1,s0+2):
    #     K=[]
    #     if i==1:
    #        for j in range(1,s0+2):
    #            if j==int(s0/2)+1:
    #               K.append(1/v)
    #            else:
    #                K.append(0)
    #     elif i==int(s0/2)+1:
    #         for j in range(1,s0+2):
    #            if j==1:
    #               K.append(-1/v)
    #            else:
    #                K.append(A1[i-2][j-2])
    #     else: 
    #         for j in range(1,s0+2):
    #            if j==1:
    #               K.append(0)
    #            else:
    #                K.append(A1[i-2][j-2])
    #     Am.append(K)

    Am=[]
    for i in range(1,s0+4):
        K=[]
        if i==1:
           for j in range(1,s0+4):
               if j==4: 
                  K.append(1/v)
               else: K.append(0)
        elif i==2: 
           for j in range(1,s0+4):
               if j==3+int(s0/2)+1: 
                  K.append(1/v)
               else: K.append(0)
        else: 
            if i==3:
                for j in range(1,s0+4):
                    if j==3+s0:
                        K.append(1/v)
                    else: K.append(0)
            else: 
                if i==4: 
                   for j in range(1,s0+4): 
                       if j==1: 
                          K.append(-1/v)
                       elif j<4:
                           K.append(0)
                       else: K.append(A1[i-4][j-4]) 
                elif i==3+int(s0/2)+1: 
                    for j in range(1,s0+4):
                        if j==2: 
                           K.append(-1/v)
                        elif j<4:
                            K.append(0)
                        else: K.append(A1[i-4][j-4])
                else: 
                    if i==3+s0: 
                       for j in range(1,s0+4):
                        if j==3: 
                           K.append(-1/v)
                        elif j<4:
                            K.append(0)
                        else: 
                            K.append(A1[i-4][j-4])
                    else: 
                        for j in range(1,s0+4):
                            if j<4: 
                               K.append(0)
                            else:
                                K.append(A1[i-4][j-4])
                                
        Am.append(K)
    A=Am
#Need to modify the two other boundary points as well
     
    return A
    
n,J,J1=8,1,1
B=np.linspace(0.02,5,100)
E,Z=[],[]
for b in B:
    e0=eig(np.array(A(b,n,J,J1)))[0]
    e0c=np.imag(e0)
    E.append(e0c)
    
    Z0=(2**n+1)*np.log(np.sinh(b*J))+(2**(n-1)-1)*np.log(np.sinh(b*J1))+(2**n+2)*np.log(2)
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

# # print(np.imag(eig(A(8.7,5,1,1))[0]))

# E=A(10,3,10,10)
# print(E)
# e0=eig(np.array(E))[0]
# e0c=np.imag(e0)
# p=1
# for e in e0c:
#     if e>=0:
#         p=p*e 

# print(p)

# # v=len(E)
# # for i in range(v):
# #     for j in range(v):
# #         if E[i][j]+E[j][i]!=0:
# #             print(i,j)