import numpy as np 
import matplotlib.pyplot as plt 
import networkx as nx 

m0,n0=2,3
m,n=2*m0+1,2*n0

A=[]
for y in range(m+1):
    if 0<y<m:  #Inner section
        for x in range(n+1):
            A.append([[x,y,1],[x,y,2]])
            A.append([[x,y,2],[x,y,3]])
            A.append([[x,y,3],[x,y,1]])

            if 1<=x<=n-1:
                if (x+y)%2==0:
                   A.append([[x,y,2],[x,y-1,2]])
                else: 
                   A.append([[x,y,2],[x,y+1,2]])
                

                if x==1:
                    A.append([[x,y,1],[x-1,y,2]])
                    A.append([[x,y,3],[x+1,y,1]])

                elif x==n-1:
                    A.append([[x,y,3],[x+1,y,2]])
                    A.append([[x,y,1],[x-1,y,3]])

                else:
                    A.append([[x,y,1],[x-1,y,3]])
                    A.append([[x,y,3],[x+1,y,1]])
               
            elif x==0:

                A.append([[x,y,1],[x,y+1,3]])
                A.append([[x,y,2],[x+1,y,1]])
                A.append([[x,y,3],[x,y-1,1]])

            else:

                A.append([[x,y,1],[x,y+1,3]])
                A.append([[x,y,2],[x-1,y,3]])
                A.append([[x,y,3],[x,y-1,1]])

        
    elif y==0:
        for x in range(1,n0+1):

           A.append([[x,y,1],[x,y,2]])
           A.append([[x,y,2],[x,y,3]])
           A.append([[x,y,3],[x,y,1]])
           
           if 1<x<n0:
                A.append([[x,y,1],[x-1,y,3]])
                A.append([[x,y,2],[2*x-1,y+1,2]])
                A.append([[x,y,3],[x+1,y,1]])
  
           elif x==1:
               A.append([[x,y,1],[0,y+1,3]])
               A.append([[x,y,2],[2*x-1,y+1,2]])
               A.append([[x,y,3],[x+1,y,1]])

           else:
               A.append([[x,y,1],[x-1,y,3]])
               A.append([[x,y,2],[2*x-1,y+1,2]])
               A.append([[x,y,3],[n,y+1,3]])
        
    else:

        for x in range(1,n0+1):

            A.append([[x,y,1],[x,y,2]])
            A.append([[x,y,2],[x,y,3]])
            A.append([[x,y,3],[x,y,1]])

            if 1<x<n0:
                A.append([[x,y,1],[x-1,y,3]])
                A.append([[x,y,2],[2*x-1,y-1,2]])
                A.append([[x,y,3],[x+1,y,1]])
  
            elif x==1:
               A.append([[x,y,1],[0,1,3]])
               A.append([[x,y,2],[2*x-1,y-1,2]])
               A.append([[x,y,3],[x+1,y,1]])

            else:
               
               A.append([[x,y,1],[x-1,y,3]])
               A.append([[x,y,2],[2*x-1,y-1,2]])
               A.append([[x,y,3],[n,y-1,1]])
           
def Enum(x,y,k):
    if y==0: 
       a=3*(x-1)+k
    elif y==n:
       a=3*n0+(2*m)*(n)*3+(x-1)*3+k
    else: 
       a=3*n0+(y-1)*(n)*3+3*x+k
    
    return a

E=[]

for ed in A:
    a1,a2=ed[0],ed[1]
    c1,c2=Enum(a1[0],a1[1],a1[2]),Enum(a2[0],a2[1],a2[2])
    E.append((c1,c2))

G=nx.Graph(E)
G.add_edges_from(E)
nx.draw(G, with_labels=True)
plt.show()
    