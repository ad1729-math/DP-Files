import numpy as np 
import matplotlib.pyplot as plt
from scipy.integrate import quad

def D(x,y,z,v1,v2,v3):
    return 1+(v1*v2)**2+(v2*v3)**2+(v3*v1)**2-2*((1-v1**2)*v2*v3*np.cos(x)+(1-v2**2)*v3*v1*np.cos(y)+(1-v3**2)*v1*v2*np.cos(z))

J=1

def Z(b):
    v1,v2,v3=np.tanh(J*b),np.tanh(J*b),1

    def D1(y):
        def D0(x):
            return np.log(D(x,y,2*np.pi-(x+y),v1,v2,v3))
        return quad(D0,0,2*np.pi)[0]
    
    return quad(D1,0,2*np.pi)[0]

B=np.linspace(0.01,10,1000)
Z0=[]
for b in B:
    Z0.append(Z(b))

plt.plot(B,Z0,'r')
plt.show()