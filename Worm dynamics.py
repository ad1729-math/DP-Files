import networkx as nx
import matplotlib.pyplot as plt
from numpy import random
import numpy as np 

m, n = 100, 100

# Generate grid points
P = [(i, j) for i in range(m) for j in range(n)]

# Initialize the graph
G = nx.Graph()
G0=nx.Graph()

# Add nodes
G.add_nodes_from(P)

def Adj(i,j):
    if i==0:
        if j==0:
           return [(0,1),(1,0)]
        elif j==n-1:
            return [(0,n-2),(1,n-1)]
        else:
            return [(0,j+1),(0,j-1),(1,j)]
    
    elif i==m-1:
        if j==0:
           return [(m-2,0),(m-1,1)]
        elif j==n-1:
            return [(m-1,n-2),(m-2,n-1)]
        else:
            return [(m-1,j+1),(m-1,j-1),(m-2,j)]
        
    else:
        if j==0:
           return [(i-1,0),(i+1,0),(i,1)]
        elif j==n-1:
            return [(i-1,n-1),(i+1,n-1),(i,n-2)]
        else:
            return [(i,j+1),(i,j-1),(i-1,j),(i+1,j)]





Worm_size=[]

It=100
En=10000

I=list(np.arange(0,It,1))

for ensemble in range(En):

    Dimer=[]

    for i in range(int(m/2)):
       for j in range(n):
           Dimer.append([(2*i,j), (2*i+1,j)])
           Dimer.append([(2*i+1,j),(2*i,j)])

    i0,j0=random.randint(0,m), random.randint(0,n)
    a0,b0=i0,j0

    for l in Adj(i0,j0):
        if [l,(i0,j0)] in Dimer:
            l0=list(l)

    x0,y0=l0[0],l0[1]

    for it in range(It):

        a,b,c,d=i0,j0,x0,y0

        Dimer.remove([(c,d),(a,b)])
        Dimer.remove([(a,b),(c,d)])

        Neig=Adj(c,d)
        Neig.remove((a,b))
        v=len(Neig)


        t=random.randint(0,v)
        s0=Neig[t]

        Dimer.append([(c,d),s0])
        Dimer.append([s0,(c,d)])

        S=list(s0)
        x,y=S[0],S[1]

        Neig2=Adj(x,y)
        
        for s1 in Neig2:
            if [s1,s0] in Dimer: 
              S1=list(s1)
              x1,y1=S1[0],S1[1]

        i0,j0=x,y
        x0,y0=x1,y1

        if (i0,j0)==(a0,b0):
            sz=it+1
            break
        else:
            sz=It
        
    Worm_size.append(sz)

Dist=[]

for i in range(It):
    c=0
    for j in Worm_size:
        if i==j:
            c+=1
    Dist.append(c)


plt.plot(I, Dist, 'bo')
plt.show()

# for l in Dimer:
#     G.add_edge(l[0],l[1])


# for l0 in Dimer0:
#     G0.add_edge(l0[0],l0[1])

# # Draw nodes
# pos = {node: node for node in G.nodes()}
# nx.draw_networkx_nodes(G, pos, node_size=30, node_color='red')

# # Draw edges
# nx.draw_networkx_edges(G, pos, edge_color='black')


# plt.axis('equal')
# plt.show()
