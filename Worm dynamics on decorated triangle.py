import numpy as np 
import matplotlib.pyplot as plt 
import networkx as nx
from matplotlib.widgets import Button
from matplotlib.backend_bases import MouseEvent
import matplotlib.pyplot as plt
import netgraph 
from pyvis.network import Network
from numpy import random

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
               A.append([[x,y,1],[0,y-1,1]])
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

def Rev_Enum(c):

    def K(c0):
       if c0%3==0: k=3
       else: k=c0%3

       return k

    def Place(c0):
        if c0<=n0:
           return c0,0
        
        elif c0>(m-1)*(n+1)+n0: 
           return c0-n0-(n-1)*(n+1)+n0 , m
        
        else: 
            y=int((c0-n0)/(n+1))+1
            x=c0-n0-(y-1)*(n+1)-1

            return x,y

    return [Place(int((c-1)/3)+1)[0],Place(int((c-1)/3)+1)[1], K(c)]

    
E=[]

for ed in A:
    a1,a2=ed[0],ed[1]
    c1,c2=Enum(a1[0],a1[1],a1[2]),Enum(a2[0],a2[1],a2[2])
    E.append((c1,c2))


E_list=[]
for e in E:
    E_list.append(list(e))

E1=[] #Duplicating

for e in E:
    E1.append(e)

def Adj(c):
    Ad=[]
    for e in E_list:
        if c==e[0]: 
           Ad.append(e[1])

    return Ad


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


#Worm size 
            
#Put weights to the dimers
            
def W(c1,c2,J,b0):
    a,b=Rev_Enum(c1),Rev_Enum(c2)

    if a[0]==b[0] and a[1]==b[1]:
        return np.exp(-b0*J)

    elif c2 in Adj(c1):
        return np.exp(b0*J)
    
    else: 
        return 0

b0,J=0.43,1

It=500
en=1000

Worm_size=[]

for e in range(en):

    Dimer=[]

    for e in E1:
     Dimer.append(list(e))

    x,y,k=random.randint(1,n-1), random.randint(1,m-1), random.randint(1,3)
    c0=2 #Enum(x,y,k)

    for l in Adj(c0):
        if [l,c0] in Dimer:
            l0=l

    c1=l0

    a,b=c0,c1

    for it in range(It):

        Dimer.remove([b,a])
        Dimer.remove([a,b])

        Neig=Adj(b)
        Neig.remove(a)
        v=len(Neig)
        
        if v>0:
            
            Weight,We=[],[]
            for s in Neig:
                Weight.append([W(s,b,b0,J),s])
                We.append(W(s,b,b0,J))
            
            Norm=np.sum(We)

            P=[]
            for w in We:
                P.append(w/Norm)


            P1=[0]
            sum=0
            for vals in P:
                sum+=vals
                P1.append(sum)
            
            ra=random.random()

            for t in range(len(P)):
                if P1[t]<ra<P1[t+1]:
                   s0=Neig[t]
                   break
                else:
                    continue
            
            #t=random.randint(0,v) #Need to give weights here
            s0=Neig[t]

            Dimer.append([b,s0])
            Dimer.append([s0,b])

            Neig2=Adj(s0)
            Neig2.remove(b)
            
            for s1 in Neig2:
                if [s1,s0] in Dimer: 
                    x1=s1

            a=s0
            b=x1

            if a==c0:
                sz=it+1
                break
            else:
                sz=It
        
        else:
            sz=it
            break

    Worm_size.append(sz)

I=list(np.arange(0,It+1,1))

Dist=[]

for i in range(It+1):
    c=0
    for j in Worm_size:
        if i==j:
            c+=1
    Dist.append(c/en)

plt.plot(I, Dist, 'b+')
plt.show()

    
# G=nx.Graph()

# Dim=[]
# for l in Dimer:
#     Dim.append((l[0],l[1]))

# G.add_edges_from(Dim)

# nx.draw(G, with_labels=True, node_color='skyblue', node_size=100, font_size=12, font_weight='bold')

# plt.show()
    
