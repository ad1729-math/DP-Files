import math as m
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from numpy.linalg import eig, det

#Total number points upto generation g
def cmax(n,g):
    if g==-1:
        return 0
    elif g==0:
        return n
    else:
        n0,n1=0,n
        c0=n
        for i in range(1,g+1):
            a,b=n0,n1
            n0=b
            n1=(n-4)*b-a
            c0+=n0+n1
    return c0

#Special treatment for g=0
def sigma(g):
    if g==0:
        return 0
    else:
        return 1

#Differing the looping part
def f(n,c,j):
    if j==0:
       if c==n:
           return 1
       else:
           return c+1
    else:
        if c<cmax(n,j):
            return c+1
        else:
            return cmax(n,j-1)+1
        
#Determining signs layer by layer
def Dir(n,s):
    a=int((n-3)/2)
    if s==1:
        if a%2==0:
           v=a
        else:
           v=a+1
    else:
        if a%2==0:
           v=a-1
        else:
           v=a

    return v

#Construction of the (n,3) Hyperbolic graph with Pfaffian orientation
        
def Hyp(n, g):
    E,Eo = [],[]  #Edges and edges with orientations
    vertices_list=[[]]
    c = 0

    for i in range(1,n+1):  #This only works for n=odd. For even n we just need to tweak a little. The orientation in that case must be 
        if i < n:           #similarly given as we have written for higher layers, i.e. odd number of edges have '1' orientation. The 
            c += 1          #rest will follow similarly.
            vertices_list[0].append([i,c])
            E.append((c,c+1))
            Eo.append([(c,c+1),1])
        else:
            c += 1
            vertices_list[0].append([i,c])
            E.append((c,1))
            Eo.append([(c,1),1])

    for j in range(1,g+1):
        new_vertices = []
        for v in vertices_list[-1]:
            c0=v[-1]
            if v[-2]!=0:
                if v[-2]!=(n-4)*sigma(j-1):
                    for i in range(0, n-3):
                        c += 1
                        if i==0:
                           E.append((c0,c))
                           Eo.append([(c0,c),1])
                        if c<cmax(n,j):
                            nv = v[:-1]  # Remove the last element of v
                            nv.append(i)
                            nv.append(c)
                            new_vertices.append(nv)
                            E.append((c,c+1)) # Putting conditions in order to construct Pfaffian orientation
                            sgn=Eo[E.index((c0,f(n,c0,j-1)))][1]
                            if i<Dir(n,sgn): #int((n-3)/2)+(1+sgn)/2:
                               Eo.append([(c,c+1),1])
                            else:
                               Eo.append([(c,c+1),-1])
                        else:
                            nv = v[:-1]  # Remove the last element of v
                            nv.append(i)
                            nv.append(c)
                            new_vertices.append(nv)
                            E.append((c,cmax(n,j-1)+1))
                            Eo.append([(c,cmax(n,j-1)+1),-1])
                else:
                    for i in range(0, n-4):  #These are the points just to the left of an n0 point
                        c += 1
                        if i==0:
                            E.append((c0,c))
                            Eo.append([(c0,c),1])
                        if c<cmax(n,j):
                            nv = v[:-1]  # Remove the last element of v
                            nv.append(i)
                            nv.append(c)
                            new_vertices.append(nv)
                            E.append((c,c+1))
                            sgn=Eo[E.index((c0,f(n,c0,j-1)))][1]
                            if i<Dir(n,sgn): #int((n-3)/2)+(1+sgn)/2:
                                Eo.append([(c,c+1),1])
                            else:
                                Eo.append([(c,c+1),-1])
                        else:
                            nv=v[:-1]  # Remove the last element of v
                            nv.append(i)
                            nv.append(c)
                            new_vertices.append(nv)
                            E.append((c,cmax(n,j-1)+1))
                            Eo.append([(c,cmax(n,j-1)+1),-1])
                            
        vertices_list.append(new_vertices)

    return vertices_list, E , Eo

#Picture of the hyperbolic lattice

# a=int(input("Enter number of generations"))
# E,Eo=Hyp(7,a)[1],Hyp(7,a)[2]
# G=nx.Graph(E)
# G.add_edges_from(E)
# nx.draw(G, with_labels=True)
# print(E)
# #print(Hyp(7,a)[0])
# plt.show()

#Building the Pfaffian matrix and calculate the eigenspectrum (May be Sturm sequnce method)

def A(b, I ,n , g):
    K=1/np.tanh(b*I)
    v=cmax(n,g)
    L=Hyp(n,g)
    Ver,Edges,Eo=L[0],L[1],L[2]
    A=[]

    Vertices=[]
    for j in range(0,g+1):
        for vals in Ver[j]:
            Vertices.append(vals)
    
    def fl(i,j):
        if j==0:
           if i==1:
               return n
           else: return i-1
        else: 
            if i==cmax(n,j-1)+1:
                return cmax(n,j)
            else:
                return i-1
    
    def fr(i,j):
        if j==0:
           if i<n:
               return i+1
           else: 
               return 1
        else: 
            if i==cmax(n,j):
                return cmax(n,j-1)+1
            else:
                return i+1

# #Modifying the graph with Fisher construction

    for gen in range(0, g):
        for i in range(cmax(n, gen-1)+1, cmax(n, gen)+1):
            B1, B2, B3 = [], [], []
        
            for j in range(1, v + 1):
                if j <= cmax(n, g - 1):  # i,j both not last layers
                    if Vertices[i - 1][-2]!= 0:
                        if j==fl(i, gen):
                            w=K*Eo[Edges.index((j, i))][1]
                            B1 += [0, 0, -w]
                            B2 += [0, 0, 0]
                            B3 += [0, 0, 0]
                          
                        elif j==fr(i, gen):
                            w=K*Eo[Edges.index((i, j))][1]
                            B1 += [0, 0, 0]
                            B2 += [0, 0, 0]
                            B3 += [w, 0, 0]
                           
                        else:
                            if (i, j) in Edges:
                                w=K*Eo[Edges.index((i, j))][1]
                                B1 += [0, 0, 0]
                                B2 += [0, w, 0]
                                B3 += [0, 0, 0]

                            else:
                                if j == i:
                                    B1 += [0, 1, -1]
                                    B2 += [-1, 0, 1]
                                    B3 += [1, -1, 0]
                                    
                                else:
                                    B1 += [0, 0, 0]
                                    B2 += [0, 0, 0]
                                    B3 += [0, 0, 0]
                                    
                    else:
                        if j == fl(i, gen):
                            w = K * Eo[Edges.index((j, i))][1]
                            B1 += [0, 0, -w]
                            B2 += [0, 0, 0]
                            B3 += [0, 0, 0]
                           
                        elif j == fr(i, gen):
                            w = K * Eo[Edges.index((i, j))][1]
                            B1 += [0, 0, 0]
                            B2 += [0, 0, 0]
                            B3 += [w, 0, 0]
                            
                        else:   
                            if (j, i) in Edges:
                                w = K * Eo[Edges.index((j,i))][1]
                                B1 += [0, 0, 0]
                                B2 += [0, -w, 0]
                                B3 += [0, 0, 0]
                               
                            else:
                                if j == i:
                                    B1 += [0, -1, 1]
                                    B2 += [1, 0, -1]
                                    B3 += [-1, 1, 0]
                                    
                                else:
                                    B1 += [0, 0, 0]
                                    B2 += [0, 0, 0]
                                    B3 += [0, 0, 0]
                else:  # j is in the last layer
                    if Vertices[j - 1][-2]!= 0:
                        B1 += [0, 0]
                        B2 += [0, 0]
                        B3 += [0, 0]
                    else:
                        if (i, j) in Edges:
                            w=K*Eo[Edges.index((i,j))][1]
                            B1 += [0, 0, 0]
                            B2 += [0, w, 0]
                            B3 += [0, 0, 0]
                        else:
                            B1 += [0, 0, 0]
                            B2 += [0, 0, 0]
                            B3 += [0, 0, 0]
                  
            A += [B1, B2, B3]

    for i in range(cmax(n,g-1)+1,cmax(n,g)+1):
        if Vertices[i - 1][-2] == 0:
            B1, B2, B3 = [], [], []
            for j in range(1, v + 1):
                if j == fl(i, g):
                    w = K*Eo[Edges.index((j, i))][1]
                    B1 += [0, -w]
                    B2 += [0, 0]
                    B3 += [0, 0]
                elif j == fr(i, g):
                    w = K * Eo[Edges.index((i, j))][1]
                    B1 += [0, 0]
                    B2 += [0, 0]
                    B3 += [w, 0]
                else:
                    if (j, i) in Edges:
                        w = K * Eo[Edges.index((j, i))][1]
                        B1 += [0, 0, 0]
                        B2 += [0, -w, 0]
                        B3 += [0, 0, 0]
                    else:
                        if j == i:
                            B1 += [0, -1, 1]
                            B2 += [1, 0, -1]
                            B3 += [-1, 1, 0]
                        else:
                            if j>cmax(n,g-1):
                                if Vertices[j-1][-2]==0:
                                   B1 += [0, 0, 0]
                                   B2 += [0, 0, 0]
                                   B3 += [0, 0, 0]
                                else:
                                    B1 += [0, 0]
                                    B2 += [0, 0]
                                    B3 += [0, 0]
                            else:
                                B1 += [0, 0, 0]
                                B2 += [0, 0, 0]
                                B3 += [0, 0, 0]

            A += [B1, B2, B3]

        else:
            B1, B2 = [], []
            for j in range(1, v + 1):
                if Vertices[j-1][-2]==0:
                    if j==fl(i, g):
                        w=K*Eo[Edges.index((j, i))][1]
                        B1 += [0, 0, -w]
                        B2 += [0, 0,  0]
                    elif j==fr(i, g):
                        w=K*Eo[Edges.index((i, j))][1]
                        B1 += [0, 0, 0]
                        B2 += [w, 0, 0]
                    else:
                        B1+=[0,0,0]
                        B2+=[0,0,0]
                else:
                    if j==fl(i, g):
                        w=K*Eo[Edges.index((j, i))][1]
                        B1 += [0, -w]
                        B2 += [0, 0]
                    elif j==fr(i,g):
                        w=K*Eo[Edges.index((i, j))][1]
                        B1 += [0,  0]
                        B2 += [w, 0]
                    else:
                        if j==i:
                            B1+=[0,-1]
                            B2+=[1,0]
                        else:
                            if j>cmax(n,g-1):
                                B1 += [0,0]
                                B2 += [0,0]
                            else:
                                B1 += [0,0,0]
                                B2 += [0,0,0]
            A += [B1, B2]

    return A

#Hermitization of the matrix

# def Ah(b,I,n,g):
#     Ah=[]
#     L=A(b,I,n,g)
#     for i in range(len(L)):
#         Ah1=[]
#         for j in range(len(L)):
#             Ah1.append(complex(0,1)*L[i][j])
#         Ah.append(Ah1)
#     return Ah

# L=A(1,1,7,3)
# v=len(L)
# c=0
# MM=[]
# for i in range(v):
#     if len(L[i])!=len(L):
#         c+=1
#         MM.append([i,len(L[i])])
# print(c, MM)

#Plotting th eigenspectrum

B=np.linspace(0.1,10,100)
E,E1,Pf,Pf1=[],[],[],[]
for b in B:
    Pfaff=np.array(A(b,1,7,2))
    Pfaff1=np.array(A(b,1,7,1))
    e0=eig(Pfaff)[0]
    e01=eig(Pfaff1)[0]
    e0c=np.imag(e0)
    e0c1=np.imag(e01)

    p,p1=1,1
    for e in e0c:
        if e>=0:
            p=p*e
        else:
            p=p

    for e1 in e0c1:
        if e1>=0:
            p1=p1*e1
        else:
            p1=p1

    E.append(e0c)
    E1.append(e0c1)
    Pf.append(np.log(p))
    Pf1.append(np.log(p1))

    
# plt.plot(B,E,'r+')
# plt.plot(B,E1,'r+')
plt.plot(B,Pf,'g+')
plt.plot(B,Pf1,'r+')
plt.plot(B,B*0,'b')
plt.xlabel("Beta--->")
plt.ylabel("Eigenspectrum--->")
# plt.ylim([-10,10])
plt.legend()
plt.show()

