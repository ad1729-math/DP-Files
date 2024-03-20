import numpy as np 
import matplotlib.pyplot as plt 
import cmath as cm 
from numpy.linalg import eig, det
from scipy.integrate import dblquad

#Ferromagnetic Ising model with plus-minus field on one but all plaquettes.

h0=10**-3

def Zl(b,J,h): 
    a=2**(-1/4)
    d,g,t,e=np.exp(b*J),np.exp(-b*J/2),a*np.exp(-b*h/4),np.sqrt(a)*np.exp(b*h/8)

    def A(x,y): 
        j=complex(0,1)
        # w1,w2=np.exp(2*j*np.pi*r/m),np.exp(2*j*np.pi*s/n)
        w1,w2=np.exp(j*x),np.exp(j*y)
        A0=[[0,-g,g,0,0, 0,-d/w1,0,0,0, 0,0,0,0,0, 0,0,0,0,0],[g,0,-1,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,-e],[-g,1,0,0,0, 0,0,0,0,0, 0,0,e,0,0, 0,0,0,0,0],
        [0,0,0,0,-g, g,0,0,0,d/w2, 0,0,0,0,0, 0,0,0,0,0],[0,0,0,g,0, -1,0,0,0,0, 0,0,0,-e,0, 0,0,0,0,0],[0,0,0,-g,1, 0,0,0,0,0, 0,0,0,0,-e, 0,0,0,0,0 ],
        [d*w1,0,0,0,0, 0,0,-g,g,0, 0,0,0,0,0, 0,0,0,0,0],[0,0,0,0,0, 0,g,0,-1,0, 0,0,0,0,0, e,0,0,0,0],[0,0,0,0,0, 0,-g,1,0,0, 0,0,0,0,0, 0,-e,0,0,0],
        [0,0,0,-d*w2,0, 0,0,0,0,0, -g,g,0,0,0, 0,0,0,0,0],[0,0,0,0,0, 0,0,0,0,g, 0,-1,0,0,0, 0,0,e,0,0],[0,0,0,0,0, 0,0,0,0,-g, 1,0,0,0,0, 0,0,0,-e,0],
        [0,0,-e,0,0, 0,0,0,0,0, 0,0,0,t,0, 0,0,0,0,0],[0,0,0,0,e, 0,0,0,0,0, 0,0,-t,0,0, 0,0,0,0,0],[0,0,0,0,0, e,0,0,0,0, 0,0,0,0,0, -t,0,0,0,0],
        [0,0,0,0,0, 0,0,-e,0,0, 0,0,0,0,t, 0,0,0,0,0],[0,0,0,0,0, 0,0,0,e,0, 0,0,0,0,0, 0,0,-t,0,0],[0,0,0,0,0, 0,0,0,0,0, -e,0,0,0,0, 0,t,0,0,0 ],
        [0,0,0,0,0, 0,0,0,0,0, 0,e,0,0,0, 0,0,0,0,t],[0,e,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,-t,0]]
             
            
        return 0.5*np.log(abs(det(np.array(A0))))

    # s0=0
    # for r in range(m):
    #     for s in range(n):
    #        m0=abs(det(np.array(A(r,s))))
    #        s0+=1/2*np.log(m0)

    x_lower = 0
    x_upper = 2*np.pi
    y_lower = 0
    y_upper = 2*np.pi

    s0=dblquad(A, x_lower, x_upper, y_lower, y_upper)[0]

    return s0

B=np.linspace(0.1,1,300)
Mag,Mag1=[],[]
for b in B:
#    m=(Zl(b,1,h0,100,100)-Zl(b,1,0,100,100))/h0
   m1=(Zl(b,1,h0)-Zl(b,1,0))/h0
#    Mag.append(m/10**4)
   Mag1.append(m1)


a=np.log(1+np.sqrt(2))/2
plt.plot(B,Mag1,'g')
plt.axvline(x=a, color='r', linestyle='--')
plt.xlabel("$1/T$ --->")
plt.ylabel("Magnetization --->")
plt.plot()
plt.show()
   
