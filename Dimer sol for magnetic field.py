import numpy as np 
import matplotlib.pyplot as plt 
import cmath as cm 
from numpy.linalg import eig, det
from scipy.integrate import dblquad

#Ferromagnetic Ising model with plus-minus field on one but all plaquettes.

h0=10**-3
x0,y0=0,0
def Zl(b,J,h): 
    a=2**(-1/4)
    d,g,t,e=np.exp(b*J),np.exp(-b*J/2),a*np.exp(-b*h/4),np.sqrt(a)*np.exp(b*h/8)

    def A(x,y): 
        j=complex(0,1)
       # w1,w2=np.exp(2*j*np.pi*r/m),np.exp(2*j*np.pi*s/n)
        w1,w2=np.exp(2*np.pi*j*x),np.exp(2*np.pi*j*y)
        A0=[[0,-g,g,0,0, 0,-d*np.conjugate(w1),0,0,0, 0,0,0,0,0, 0,0,0,0,0],[g,0,-1,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,-e],[-g,1,0,0,0, 0,0,0,0,0, 0,0,e,0,0, 0,0,0,0,0],
        [0,0,0,0,-g, g,0,0,0,d*np.conjugate(w2), 0,0,0,0,0, 0,0,0,0,0],[0,0,0,g,0, -1,0,0,0,0, 0,0,0,-e,0, 0,0,0,0,0],[0,0,0,-g,1, 0,0,0,0,0, 0,0,0,0,-e, 0,0,0,0,0 ],
        [d*w1,0,0,0,0, 0,0,-g,g,0, 0,0,0,0,0, 0,0,0,0,0],[0,0,0,0,0, 0,g,0,-1,0, 0,0,0,0,0, e,0,0,0,0],[0,0,0,0,0, 0,-g,1,0,0, 0,0,0,0,0, 0,-e,0,0,0],
        [0,0,0,-d*w2,0, 0,0,0,0,0, -g,g,0,0,0, 0,0,0,0,0],[0,0,0,0,0, 0,0,0,0,g, 0,-1,0,0,0, 0,0,e,0,0],[0,0,0,0,0, 0,0,0,0,-g, 1,0,0,0,0, 0,0,0,-e,0],
        [0,0,-e,0,0, 0,0,0,0,0, 0,0,0,t,0, 0,0,0,0,0],[0,0,0,0,e, 0,0,0,0,0, 0,0,-t,0,0, 0,0,0,0,0],[0,0,0,0,0, e,0,0,0,0, 0,0,0,0,0, -t,0,0,0,0],
        [0,0,0,0,0, 0,0,-e,0,0, 0,0,0,0,t, 0,0,0,0,0],[0,0,0,0,0, 0,0,0,e,0, 0,0,0,0,0, 0,0,-t,0,0],[0,0,0,0,0, 0,0,0,0,0, -e,0,0,0,0, 0,t,0,0,0 ],
        [0,0,0,0,0, 0,0,0,0,0, 0,e,0,0,0, 0,0,0,0,t],[0,e,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,-t,0]]
             
        Ei=eig(np.array(A0))[0]
        E0=np.imag(Ei)
        E=[]
        for e1 in E0:
            if e1>=0:
                E.append(e1)

        return np.sort(E)[0]  #0.5*np.log(abs(det(np.array(A0)))) #E0 #np.sort(E)[0] 

    # s0=0
    # for r in range(m):
    #     for s in range(n):
    #        m0=abs(det(np.array(A(r,s))))
    #        s0+=1/2*np.log(m0)

    x_lower = 0
    x_upper = 1
    y_lower = 0
    y_upper = 1

    s0=dblquad(A, x_lower, x_upper, y_lower, y_upper)[0]

#     s0=A(x0,y0)
#     # E=[]
#     # for r in range(m):
#     #     for s in range(n):
#     #         for e in A(r,s):
#     #             E.append(e)

    return s0 #E

# a=np.log(1+np.sqrt(2))/2

# B=np.linspace(a-10**-2,a+10**-2,1000)
# Mag,Mag1=[],[]
# for b in B:
# #    m=(Zl(b,1,h0,100,100)-Zl(b,1,0,100,100))/h0
#    m=Zl(b,1,0)
# #    m1=Zl(b,1,0)  #(Zl(b,1,h0)-Zl(b,1,0))/h0
#    Mag.append(m)
# #    Mag1.append(m1)

# # E=[]
# # for b in B:
# #     E.append(Zl(b,1,2,30,30))

# # plt.plot(B,E,'r+')
# plt.plot(B,Mag,'b')
# plt.axhline(y=0, color='b', linestyle='--')
# plt.xlabel("$1/T$ --->")
# plt.ylabel("Magnetization --->")
# # plt.plot(B,E,'r+')
# # plt.plot(B,E1,'g+')
# plt.show() 

def integrand(phi1, phi2, beta, J, h):
    return np.log(1/4*(8 + 4*np.exp(-4*beta*h) + 2*np.exp(2*beta*h) + 2*np.exp(4*beta*h)*np.cosh(2*beta*J) - 8*(np.cos(phi1) + np.cos(phi2))*np.cosh(beta*J)*np.exp(2*beta*h) + 4*(np.cos(phi1-phi2) + np.cos(phi1+phi2))*(np.cos(beta*J)-1)))

def Z(beta, J, h):
    integrand_wrapper = lambda phi1, phi2: integrand(phi1, phi2, beta, J, h)
    integral, _ = dblquad(integrand_wrapper, 0, 2*np.pi, lambda x: 0, lambda x: 2*np.pi)
    return 1/(4*np.pi**2)*integral

# Example usage:
# B=np.linspace(1,5,10)
# h=10**-2
# E=[]
# for b in B:
#     a=(Z(b,1,h)-Z(b,1,0))/h
#     E.append(a)

# plt.plot(B,E,'r')
# plt.show()

print(Z(1,0,0))
print(Zl(1,0,0))
