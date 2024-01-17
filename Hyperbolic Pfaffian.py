import math as m
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt


def cmax(n,g):
    if g==0:
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


def sigma(g):
    if g==0:
        return 0
    else:
        return 1

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

def Hyp(n, g):
    E,Eo = [],[]  #Edges and edges with orientations
    vertices_list=[[]]
    c = 0

    for i in range(1,n+1):
        if i < n:
            c += 1
            vertices_list[0].append([i,c])
            E.append((c,c+1))
            Eo.append([(c,c+1),1])
        else:
            c += 1
            vertices_list[0].append([i,c])
            E.append((c,1))
            Eo.append([(c,1),1])

    j = 1
    while j <= g:
        new_vertices = []
        for v in vertices_list[-1]:
            c0=v[-1]
            if v[-2] != 0:
                if v[-2]!=(n-4)*sigma(j-1):
                    for i in range(0, n-3):
                        c += 1
                        if i==0:
                           E.append((c0,c))
                           Eo.append([(c0,c),1])
                        if c<cmax(n,g):
                            nv = v[:-1]  # Remove the last element of v
                            nv.append(i)
                            nv.append(c)
                            new_vertices.append(nv)
                            E.append((c,c+1)) #Need to put conditions 
                            sgn=Eo[E.index((c0,f(n,c0,j-1)))][1]
                            if i<int((n-3)/2)+(1+sgn)/2:
                               Eo.append([(c,c+1),1])
                            else:
                               Eo.append([(c,c+1),-1])
                            
                        else:
                            nv = v[:-1]  # Remove the last element of v
                            nv.append(i)
                            nv.append(c)
                            new_vertices.append(nv)
                            E.append((c,cmax(n,g-1)+1))
                            Eo.append([(c,cmax(n,g-1)+1),-1])
                            break
                else:
                    for i in range(0, n-4):
                        c += 1
                        if i==0:
                           E.append((c0,c))
                           Eo.append([(c0,c),1])
                        if c==cmax(n,g):
                            nv = v[:-1]  # Remove the last element of v
                            nv.append(i)
                            nv.append(c)
                            new_vertices.append(nv)
                            E.append((c,cmax(n,g-1)+1))
                            Eo.append([(c,cmax(n,g-1)+1),-1])
                            break
                        else:
                            v = v[:-1]  # Remove the last element of v
                            nv.append(i)
                            nv.append(c)
                            new_vertices.append(nv)
                            E.append((c,c+1))
                            if i<int((n-3)/2)+(1+sgn)/2:
                                Eo.append([(c,c+1),1])
                            else:
                               Eo.append([(c,c+1),-1])
        vertices_list.append(new_vertices)
        j += 1
    return vertices_list, E , Eo

# Example usage
# Print the vertices of the second layer (j=1)


E,Eo=Hyp(7,2)[1],Hyp(7,2)[2]
#E=Hyp(7,2)[1]
G=nx.Graph(E)
G.add_edges_from(E)
nx.draw(G, with_labels=True)
print(Eo)
plt.show()

#Almost correct, but a little correction is required. And need to construct Pfaffian orientation. Doit recurrsively as well.

#Building the Pfaffian matrix and calculate the eigenspectrum (May be Sturm sequnce method)

def A(b, I ,g ,n):
    v=cmax(n,g)
    Eo=Hyp(n,g)[2]
    A=[]
    Edges=[]
    for k in range(len(Eo)):
        Edges.append(Eo[k][0])
        
    for i in range(1,v+1):
        B=[]
        for j in range(1,v+1):
            if (i,j) in Edges:
                B.append(1/np.tanh(b*I)*Eo[Edges.index((i,j))][1])
            else: 
                if (j,i) in Edges:
                    B.append(-1/np.tanh(b*I)*Eo[Edges.index((j,i))][1])
                else:
                    B.append(0)
        A.append(B)
    return A 


#print(A(1,1,1,7))

