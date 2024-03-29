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

        return 0.5*np.log(abs(det(np.array(A0)))) #E0 #np.sort(E)[0] 

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

a0=np.log(1+np.sqrt(2))/2

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


def integrand1(beta,J,h,phi1,phi2):

    term1 = -beta*16*np.exp(-4*beta*h) + 8*beta*np.exp(4*beta*h) + 8*beta*np.exp(4*beta*h)*np.cosh(2*beta*J)
    term2 = -8*(np.cos(phi1) + np.cos(phi2))*(2*beta*np.exp(2*beta*h)*np.cosh(beta*J) + 2*beta*np.exp(-2*beta*h))
    result=1/4*(term1+term2)

    return result


def integrand(beta, J, h, phi1, phi2):

    term1 = 8 + 4*np.exp(-4*beta*h) + 2*np.exp(4*beta*h) + 2*np.exp(4*beta*h)*np.cosh(2*beta*J)
    term2 = -8*(np.cos(phi1) + np.cos(phi2))*(np.exp(2*beta*h)*np.cosh(beta*J) - np.exp(-2*beta*h))
    term3 = -4*(np.cos(phi1 - phi2) + np.cos(phi1 + phi2))*(np.cosh(beta*J)-1)
    result=1/4*(term1+term2+term3)

    r1=integrand1(beta,J,h,phi1,phi2)

    return result


def Z(beta, J, h):
    integrand_wrapper = lambda phi1, phi2: integrand(beta, J, h, phi1, phi2)
    integral, _ = dblquad(integrand_wrapper, 0, 2*np.pi, lambda x: 0, lambda x: 2*np.pi)
    return 1/(8*np.pi**2)*integral/(2*beta)


# Example usage:
B=np.linspace(0.1,3,200)
h=0
E=[]
for b in B:
    a=integrand(b,1,h,0,0)
    E.append(a)

plt.plot(B,E,'g')
plt.axvline(x=a0, color='b', linestyle='--')
plt.axhline(y=0, color='r', linestyle='-')
plt.xlabel("$1/T$ --->")
plt.ylabel("Magnetization --->")
plt.show()




