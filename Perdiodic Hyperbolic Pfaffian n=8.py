import math as m
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from numpy.linalg import eig, det

def Metric(gen):
    # if gen==g:
    #     return 1
    # else: 
    #     return 10**-4

    return 1/(1-np.tanh(gen/2)**2) #Scaling the interaction strength with metric
#Total number points upto generation g


def Layers(n,g):
    if g==-1:
        return 0,0,0,0
    elif g==0:
        return 0,n,0,n
    else:
        n0,n1,N0,N1=0,n,0,n
        for i in range(1,g+1):
            a,b=n0,n1
            n0=b
            n1=(n-4)*b-a
            N0+=n0
            N1+=n1
        return n0,n1,N0,N1

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

def s(n):
    if (n/2)%2!=0:
        return n/2
    else:
        return n/2+1
    
#Construction of the (n,3) Hyperbolic graph with Pfaffian orientation
        
def Hyp(n, g):
    E,Eo = [],[]  #Edges and edges with orientations
    vertices_list=[[]]
    c = 0

    if n%2!=0:
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

    else:
        for i in range(1,n+1):   
            if i < n: 
                if i<=s(n):         
                    c += 1         
                    vertices_list[0].append([i,c])
                    E.append((c,c+1))
                    Eo.append([(c,c+1),1])
                else:
                    c += 1         
                    vertices_list[0].append([i,c])
                    E.append((c,c+1))
                    Eo.append([(c,c+1),-1])

            else:
                c += 1
                vertices_list[0].append([i,c])
                E.append((c,1))
                Eo.append([(c,1),-1])
        

    for j in range(1,g+1):
        Ver=vertices_list[j-1]
        new_vertices = []
        En=[]

        for v in Ver:
            c0=v[-1]
            if v[-2]!=0:
                if Ver[f(n,c0,j-1)-cmax(n,j-2)-1][-2]!=0:
                    for i in range(0, n-3):
                        c += 1
                        if i==0:
                           E.append((c0,c))
                           Eo.append([(c0,c),1])
                           En.append(c)
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

# a=2 #int(input("Enter number of generations"))
# E,Eo=Hyp(7,a)[1],Hyp(7,a)[2]
# G=nx.Graph(E)
# G.add_edges_from(E)
# nx.draw(G, with_labels=True)
# print(Eo)
# #print(Hyp(7,a)[0])
# plt.show()

#Building the Pfaffian matrix and calculate the eigenspectrum (May be Sturm sequnce method)

def A(b, I ,n , g):

    def K(a):
        return 1/np.tanh(b*I) #1/np.tanh(b*I*Metric(a))


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
            
    def Num(n,g):
        C1,C0=1,0
        if g==0:
            v=1
        else: 
            for i in range(g-1):
                c,d=C1,C0
                C1=(n-4)*c+1
                C0=c
            v=C1+C0
        return (cmax(n,g-1)+1)+v+1
            
    n0=Num(n,g)

    L_half, R_half=[n0],[n0+1]
    i0,i1=n0,n0+1
    for j in range(0,int((cmax(n,g)-cmax(n,g-1))/2)-1):
        if Vertices[fr(i1,g)-1][-2]!=0:
            R_half.append(fr(i1,g))
        i1=fr(i1,g)

        if Vertices[fl(i0,g)-1][-2]!=0:
            L_half.append(fl(i0,g))
        i0=fl(i0,g)


    i0,i1=L_half[0],R_half[0]
    Ed, Dir= [(i0,i1)],[1]

    for j in range(1, len(L_half)):

        Ed.append([L_half[j],R_half[j]])
        p1,q1=L_half[j],R_half[j]
        p0,q0=L_half[j-1],R_half[j-1]
        s1,s2=0,0

        if p1==fl(p0,g):
            if Eo[Edges.index((p1, p0))][1]==1:
                s1+=1
            
            if Eo[Edges.index((q0, q1))][1]==1:
                s2+=1

        else: 
            if Eo[Edges.index((fl(p0,g), p0))][1]==1:
                s1+=1
                if Eo[Edges.index((p1, fl(p0,g)))][1]==1:
                    s1+=1
            elif Eo[Edges.index((p1, fl(p0,g)))][1]==-1 and Eo[Edges.index((p1, fl(p0,g)))][1]==1:
                    s1+=1
            else: 
                s1+=0 

            if Eo[Edges.index((q0, fr(q0,g)))][1]==1:
                s2+=1
                if Eo[Edges.index((fr(p0,g),p1))][1]==1:
                    s2+=1
            elif Eo[Edges.index((q0,fr(q0,g)))][1]==-1 and Eo[Edges.index((fl(q1,g),q1))][1]==1:
                    s2+=1
            else: 
                s2+=0 
        
        s=Dir[-1]
        if (s1+s2+s)%2==0:
            Dir.append(-1)
        else: 
            Dir.append(1)
                
        
    def Tor_dir(i,j):
        
        if i in L_half and j in R_half:
            if L_half.index(i)==R_half.index(j):
                return Dir[L_half.index(i)]

        elif i in R_half and j in L_half: 
            if R_half.index(i)==L_half.index(j):
                return -Dir[R_half.index(i)]

        else: 
            return 0
    
          
# #Modifying the graph with Fisher construction

    for gen in range(0, g):
        for i in range(cmax(n, gen-1)+1, cmax(n, gen)+1):
            B1, B2, B3 = [], [], []
        
            for j in range(1, v + 1):
                if j <= cmax(n, g - 1):  # i,j both not last layers
                    if Vertices[i - 1][-2]!= 0:
                        if j==fl(i, gen):
                            w=K(gen)*Eo[Edges.index((j, i))][1]
                            B1 += [0, 0, -w]
                            B2 += [0, 0, 0]
                            B3 += [0, 0, 0]
                          
                        elif j==fr(i, gen):
                            w=K(gen)*Eo[Edges.index((i, j))][1]
                            B1 += [0, 0, 0]
                            B2 += [0, 0, 0]
                            B3 += [w, 0, 0]
                           
                        else:
                            if (i, j) in Edges:
                                w=K(gen)*Eo[Edges.index((i, j))][1]
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
                            w = K(gen)*Eo[Edges.index((j, i))][1]
                            B1 += [0, 0, -w]
                            B2 += [0, 0, 0]
                            B3 += [0, 0, 0]
                           
                        elif j == fr(i, gen):
                            w = K(gen)*Eo[Edges.index((i, j))][1]
                            B1 += [0, 0, 0]
                            B2 += [0, 0, 0]
                            B3 += [w, 0, 0]
                            
                        else:   
                            if (j, i) in Edges:
                                w = K(gen)*Eo[Edges.index((j,i))][1]
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
                        B1 += [0, 0, 0]
                        B2 += [0, 0, 0] 
                        B3 += [0, 0, 0]
                        
                    else:
                        if (i, j) in Edges:
                            w= K(gen)*Eo[Edges.index((i,j))][1]
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
                    w = K(g)*Eo[Edges.index((j, i))][1]
                    B1 += [0, 0, -w]
                    B2 += [0, 0, 0]
                    B3 += [0, 0, 0] 

                elif j == fr(i, g):
                    w = K(g)*Eo[Edges.index((i, j))][1]
                    B1 += [0, 0, 0]
                    B2 += [0, 0, 0]
                    B3 += [w, 0, 0]
                else:
                    if (j, i) in Edges:
                        w = K(g)*Eo[Edges.index((j, i))][1]
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

            A += [B1, B2, B3]
        else: 
            B1,B2,B3=[],[],[]
            for j in range(1,v+1):
                if Vertices[j-1][-2]==0:
                    if j==fl(i, g):
                        w= K(g)*Eo[Edges.index((j, i))][1]
                        B1 += [0, 0, -w]
                        B2 += [0, 0,  0]
                        B3 += [0, 0, 0]


                    elif j==fr(i, g):
                        w= K(g)*Eo[Edges.index((i, j))][1]
                        B1 += [0, 0, 0]
                        B2 += [0, 0, 0]
                        B3 += [w ,0, 0]

                    else:
                        B1+=[0,0,0]
                        B2+=[0,0,0]
                        B3+=[0,0,0]

                else:
                    if j==fl(i,g):
                        if j==L_half[0]:
                            w= K(g)*Eo[Edges.index((j, i))][1]
                            w0=K(g)*Tor_dir(i,j)
                            B1 +=[0,0,-w]
                            B2 +=[0,w0,0]
                            B3 +=[0,0,0]

                        elif j==R_half[-1]:
                            w= K(g)*Eo[Edges.index((j, i))][1]
                            w0=K(g)*Tor_dir(i,j)
                            B1 +=[0,0,-w]
                            B2 +=[0,w0,0]
                            B3 +=[0,0,0]
                        else: 
                            w= K(g)*Eo[Edges.index((j, i))][1]
                            B1 +=[0,0,-w]
                            B2 +=[0,0,0]
                            B3 +=[0,0,0]

                    
                    elif j==fr(i,g):
                        if i==L_half[0]:
                            w= K(g)*Eo[Edges.index((i, j))][1]
                            w0=K(g)*Tor_dir(i,j)
                            B1 +=[0,0,0]
                            B2 +=[0,w0,0]
                            B3 +=[w,0,0]

                        elif i==R_half[-1]:
                            w= K(g)*Eo[Edges.index((i, j))][1]
                            w0=K(g)*Tor_dir(i,j)
                            B1 +=[0,0,0]
                            B2 +=[0,w0,0]
                            B3 +=[w,0,0]

                        else:
                            w= K(g)*Eo[Edges.index((i,j))][1]
                            B1 +=[0,0,0]
                            B2 +=[0,0,0]
                            B3 +=[w,0,0]

                    else:
                        if i in L_half and j in R_half:
                            if L_half.index(i)==R_half.index(j):
                                w0=K(g)*Tor_dir(i,j)
                                B1 +=[0,0,0]
                                B2 +=[0,w0,0]
                                B3 +=[0,0,0]

                            else:
                                B1 +=[0,0,0]
                                B2 +=[0,0,0]
                                B3 +=[0,0,0]

                        elif i in R_half and j in L_half: 
                            if R_half.index(i)==L_half.index(j):
                                w0=K(g)*Tor_dir(i,j)
                                B1 +=[0,0,0]
                                B2 +=[0,w0,0]
                                B3 +=[0,0,0]

                            else:
                                B1 +=[0,0,0]
                                B2 +=[0,0,0]
                                B3 +=[0,0,0]


                        else:
                            if j==i:
                                B1 += [0, -1, 1]
                                B2 += [1, 0, -1]
                                B3 += [-1, 1, 0]
                            else:
                                B1 +=[0,0,0]
                                B2 +=[0,0,0]
                                B3 +=[0,0,0]
            
            A +=[B1,B2,B3]

    return A

n,g,I=8,2,1
# n1,g1,I1=6,3,1
b=1

N=cmax(n,g)+Layers(n,g)[2]
# N1=cmax(n1,g1)+Layers(n1,g1)[2]

B=np.linspace(0.4,4,100)
E,Pf=[],[]
E1,Pf1=[],[]

for b in B:
    Pfaff=np.array(A(b,I,n,g))
   # Pfaff1=np.array(Ah(b,I1,n1,g1))
    e0=eig(Pfaff)[0]
   # e01=eig(Pfaff1)[0]
    e0c=np.imag(e0)
    # e0c1=np.imag(e01)

    s0=0
    for e in e0c:
        if e>=0:
           s0 +=np.log(e)
  
    # s1=0
    # for e1 in e0c1: 
    #     if e1>=0:
    #        s1 +=np.log(e1)

    E.append(e0c)
   # E1.append(e0c1)
    Pf.append((s0+N*np.log(np.sinh(b*I))+cmax(n,g)*np.log(2)))
   # Pf1.append(s1+N1*np.log(np.sinh(b*I1))+cmax(n1,g1)*np.log(2))

plt.plot(B,E,'g+')
#plt.plot(B,E1,'g+')
#plt.plot(B,Pf,'g+')
#plt.plot(B,Pf1,'r+')
plt.plot(B,B*0,'b')
plt.xlabel("Beta--->")
plt.ylabel("Eigenspectrum--->")
# plt.ylim([-10,10])
#plt.legend()
plt.show()


