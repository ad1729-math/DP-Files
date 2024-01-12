import numpy as np 
import matplotlib.pyplot as plt 

b=2
q=2
a,w1,w2=q/(q**2+1),1/q**2*(1-1/q**3),1/q**3*(1-1/q)

def P(t,p):
    L=[[1,1]]
    for i in range(0,t):
        L0=L[-1]
        v1=((p*w1+(1-p))*L0[0]+(p*w2+a*(1-p))*L0[1])**b
        v2=((p*w1+a*(1-p))*L0[0]+(p*w2+a*(1-p))*L0[1])**b
        L.append([v1,v2])
    return L 

def S(t,p):
    S=[a]
    S1=[np.log(a)]
    P0=P(t,p)
    for i in range(t):
        v=S[-1]*(a*(1-p)+w2*p)*(((1-p)*a+w1*p)*P0[i][0]+((1-p)*a+p*w2)*P0[i][1])**(b-1)
        S.append(v)
        S1.append(np.log(v))
    return S,S1


Pl=np.linspace(0,1,100)
Sl1,Sl2=[],[]
for p in Pl:
    Sl1.append(S(10,p)[1][-1])
    Sl2.append(S(10,p)[1][5])


# T=np.arange(0,11,1)
# S0=S(10,0.28)[1]
# plt.plot(T,S0,'ro')
plt.plot(Pl,Sl1,'ro')   #Phase transition is within p=0.27 ~ 0.28
plt.plot(Pl,Sl2,'bo')
plt.show()





        


    
    
