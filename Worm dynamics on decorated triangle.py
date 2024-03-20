import numpy as np 
import matplotlib.pyplot as plt 
import networkx as nx
from matplotlib.widgets import Button
from matplotlib.backend_bases import MouseEvent
import matplotlib.pyplot as plt
import netgraph 
from pyvis.network import Network

m0,n0=2,3
m,n=2*m0+1,2*n0

A=[]
for y in range(m+1):
    if 0<y<m:  #Inner section
        for x in range(n+1):
            A.append([[x,y,1],[x,y,2]])
            A.append([[x,y,2],[x,y,3]])
            A.append([[x,y,3],[x,y,1]])

            A.append([[x,y,2],[x,y,1]])
            A.append([[x,y,3],[x,y,2]])
            A.append([[x,y,1],[x,y,3]])

            if 1<=x<=n-1:
                if (x+y)%2==0:
                   if 1<y<=m-1:
                     A.append([[x,y,2],[x,y-1,2]])
                   else:
                     A.append([[x,y,2],[int((x+1)/2),y-1,2]])

                else: 
                    if 1<=y<m-1:
                      A.append([[x,y,2],[x,y+1,2]])
                    else:
                      A.append([[x,y,2],[int((x+1)/2),y+1,2]])
                

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

                if y==1:
                    A.append([[x,y,1],[x,y+1,3]])
                    A.append([[x,y,2],[x+1,y,1]])
                    A.append([[x,y,3],[1,y-1,1]])

                elif y==m-1:
                    A.append([[x,y,1],[1,y+1,1]])
                    A.append([[x,y,2],[x+1,y,1]])
                    A.append([[x,y,3],[x,y-1,1]])

                else:
                    A.append([[x,y,1],[x,y+1,3]])
                    A.append([[x,y,2],[x+1,y,1]])
                    A.append([[x,y,3],[x,y-1,1]])


            else:

                if y==1:
                    A.append([[x,y,1],[x,y+1,3]])
                    A.append([[x,y,2],[x-1,y,3]])
                    A.append([[x,y,3],[n0,y-1,3]])

                elif y==m-1:
                    A.append([[x,y,1],[n0,y+1,3]])
                    A.append([[x,y,2],[x-1,y,3]])
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


           A.append([[x,y,2],[x,y,1]])
           A.append([[x,y,3],[x,y,2]])
           A.append([[x,y,1],[x,y,3]])
            
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


            A.append([[x,y,2],[x,y,1]])
            A.append([[x,y,3],[x,y,2]])
            A.append([[x,y,1],[x,y,3]])

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
    elif y==m:
       a=3*n0+(2*m0)*(n+1)*3+(x-1)*3+k
    else: 
       a=3*n0+(y-1)*(n+1)*3+3*x+k

    return a


E=[]

for ed in A:
    a1,a2=ed[0],ed[1]
    c1,c2=Enum(a1[0],a1[1],a1[2]),Enum(a2[0],a2[1],a2[2])
    E.append((c1,c2))


E_list=[]
for e in E:
    E_list.append(list(e))

E1=E #Duplicating

def Adj(c):
    Ad=[]
    for e in E_list:
        if c==e[0]: 
           Ad.append(e[1])

    return Ad

Dim0=[]

for y in range(0,m+1):

    if y==0:
      for x in range(1,n0+1):
            c1,c2,c3=Enum(x,y,1),Enum(x,y,2),Enum(x,y,3)
            
            E1.remove((c1,c2))
            E1.remove((c2,c3))
            E1.remove((c3,c1))
            E1.remove((c3,c2))
            E1.remove((c2,c1))
            E1.remove((c1,c3))
    elif y==m:
        for x in range(1,n0+1):
            c1,c2,c3=Enum(x,y,1),Enum(x,y,2),Enum(x,y,3)

            E1.remove((c1,c2))
            E1.remove((c2,c3))
            E1.remove((c3,c1))
            E1.remove((c3,c2))
            E1.remove((c2,c1))
            E1.remove((c1,c3))
                
    else:
        for x in range(0,n+1):
            c1,c2,c3=Enum(x,y,1),Enum(x,y,2),Enum(x,y,3)
            
            E1.remove((c1,c2))
            E1.remove((c2,c3))
            E1.remove((c3,c1))
            E1.remove((c3,c2))
            E1.remove((c2,c1))
            E1.remove((c1,c3))

        
G=nx.Graph(E1)
G.add_edges_from(E1)

nx.draw(G,  with_labels=True)
plt.show() 
   
print(E1)