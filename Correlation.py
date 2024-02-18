import numpy as np 
import math as m 
import networkx as nx 
import matplotlib.pyplot as plt 
from numpy.linalg import eig, det
#from Hyperbolic_Pfaffian import A

gc=1

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
        
def Hyp_corr(n,g,gc,d):

    #Modifications for finding the correlation
    x0,y0=cmax(n,gc-1)+2,cmax(n,gc-1)+2+d+1

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
        
        for v in Ver:
            c0=v[-1]
            if v[-2]!=0:
                if Ver[f(n,c0,j-1)-cmax(n,j-2)-1][-2]!=0:
                    for i in range(0, n-3):
                        c += 1
                        if i==0:  
                           E.append((c0,c))
                           if x0<c0<y0:
                              Eo.append([(c0,c),-1])
                           else:
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
                            if x0<c0<y0:
                               E.append((c0,c))
                               Eo.append([(c0,c),-1])
                            else:
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
# E,Eo=Hyp_corr(6,a)[1],Hyp_corr(6,a)[2]
# G=nx.Graph(E)
# G.add_edges_from(E)
# nx.draw(G, with_labels=True)
# print(Eo)
# #print(Hyp(7,a)[0])
# plt.show()

#Building the Pfaffian matrix and calculate the eigenspectrum (May be Sturm sequnce method)

def ACorr(b, I , n , g , gc, d):

    x0,y0=cmax(n,gc-1)+2,cmax(n,gc-1)+2+d+1

    def K(a):
        return 1/np.tanh(b*I) # 1/np.tanh(b*I*Metric(a))
    
    v=cmax(n,g)
    L=Hyp_corr(n,g ,gc,d)
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

            if i==x0:  #X point
               
                B1,B2,B3,B4,B5,B6=[],[],[],[],[],[]

                for j in range(1, v + 1):
                    if j <= cmax(n, g - 1):  # i,j both not last layer
                        if j==fl(i, gen):
                            w=K(gen)*Eo[Edges.index((j, i))][1]
                            B1 += [0, 0, -w]
                            B2 += [0, 0, 0]
                            B3 += [0, 0, 0]
                            B4 += [0, 0, 0]
                            B5 += [0, 0, 0]
                            B6 += [0, 0, 0]
                        
                        elif j==fr(i, gen):
                            w=K(gen)*Eo[Edges.index((i, j))][1]
                            B1 += [0, 0, 0]
                            B2 += [0, 0, 0]
                            B3 += [w, 0, 0]
                            B4 += [0, 0, 0]
                            B5 += [0, 0, 0]
                            B6 += [0, 0, 0]
                        
                        else:
                            if (i, j) in Edges:
                                w=K(gen)*Eo[Edges.index((i, j))][1]
                                B1 += [0, 0, 0]
                                B2 += [0, w, 0]
                                B3 += [0, 0, 0]
                                B4 += [0, 0, 0]
                                B5 += [0, 0, 0]
                                B6 += [0, 0, 0]

                            elif (j,i) in Edges:
                                w=K(gen)*Eo[Edges.index((j, i))][1]
                                B1 += [0, 0, 0]
                                B2 += [0, -w, 0]
                                B3 += [0, 0, 0]
                                B4 += [0, 0, 0]
                                B5 += [0, 0, 0]
                                B6 += [0, 0, 0]

                            else:
                                if j == i:
                                    B1 += [0, 0, 0, 1, 0, 0]
                                    B2 += [0, 0, 0, 0, -1, 0]
                                    B3 += [0, 0, 0, 0, 0, -1]
                                    B4 += [-1, 0, 0, 0, 1, -1]
                                    B5 += [0, 1, 0, -1, 0, 1]
                                    B6 += [0, 0, 1, 1, -1, 0]
                                    
                                else:
                                    if j==y0:
                                        B1 += [0, 0, 0, 0, 0, 0]
                                        B2 += [0, 0, 0, 0, 0, 0]
                                        B3 += [0, 0, 0, 0, 0, 0]
                                        B4 += [0, 0, 0, 0, 0, 0]
                                        B5 += [0, 0, 0, 0, 0, 0]
                                        B6 += [0, 0, 0, 0, 0, 0]
                                    else:
                                        B1 += [0, 0, 0]
                                        B2 += [0, 0, 0]
                                        B3 += [0, 0, 0]
                                        B4 += [0, 0, 0]
                                        B5 += [0, 0, 0]
                                        B6 += [0, 0, 0]
                                        
                                        
                    else:  # j is in the last layer
                        if Vertices[j - 1][-2]!= 0:
                            B1 += [0, 0]
                            B2 += [0, 0]
                            B3 += [0, 0]
                            B4 += [0, 0]
                            B5 += [0, 0]
                            B6 += [0, 0]

                            # B1 += [0, 0 ,0]
                            # B2 += [0, 0 ,0]
                            # B3 += [0, 0 ,0]
                        else:
                            if (i, j) in Edges:
                                w=K(gen)*Eo[Edges.index((i,j))][1]
                                B1 += [0, 0, 0]
                                B2 += [0, w, 0]
                                B3 += [0, 0, 0]
                                B4 += [0, 0, 0]
                                B5 += [0, 0, 0]
                                B6 += [0, 0, 0]
                            else:
                                B1 += [0, 0, 0]
                                B2 += [0, 0, 0]
                                B3 += [0, 0, 0]
                                B4 += [0, 0, 0]
                                B5 += [0, 0, 0]
                                B6 += [0, 0, 0]
                    
                A += [B1, B2, B3, B4, B5, B6]


            elif i==y0:  #Y point 

                B1,B2,B3,B4,B5,B6=[],[],[],[],[],[]
                for j in range(1, v + 1):
                    if j <= cmax(n, g - 1):  # i,j both not last layer
                        if j==fl(i, gen):
                            w=K(gen)*Eo[Edges.index((j, i))][1]
                            B1 += [0, 0, -w]
                            B2 += [0, 0, 0]
                            B3 += [0, 0, 0]
                            B4 += [0, 0, 0]
                            B5 += [0, 0, 0]
                            B6 += [0, 0, 0]
                        
                        elif j==fr(i, gen):
                            w=K(gen)*Eo[Edges.index((i, j))][1]
                            B1 += [0, 0, 0]
                            B2 += [0, 0, 0]
                            B3 += [w, 0, 0]
                            B4 += [0, 0, 0]
                            B5 += [0, 0, 0]
                            B6 += [0, 0, 0]
                        
                        else:
                            if (i, j) in Edges:
                                w=K(gen)*Eo[Edges.index((i, j))][1]
                                B1 += [0, 0, 0]
                                B2 += [0, w, 0]
                                B3 += [0, 0, 0]
                                B4 += [0, 0, 0]
                                B5 += [0, 0, 0]
                                B6 += [0, 0, 0]

                            elif (j, i) in Edges:
                                w=K(gen)*Eo[Edges.index((j, i))][1]
                                B1 += [0, 0, 0]
                                B2 += [0, -w, 0]
                                B3 += [0, 0, 0]
                                B4 += [0, 0, 0]
                                B5 += [0, 0, 0]
                                B6 += [0, 0, 0]

                            else:
                                if j == i:
                                    B1 += [0, 0, 0, 1, 0, 0]
                                    B2 += [0, 0, 0, 0, 1, 0]
                                    B3 += [0, 0, 0, 0, 0, -1]
                                    B4 += [-1, 0, 0, 0, 1, -1]
                                    B5 += [0, -1, 0, -1, 0, 1]
                                    B6 += [0, 0, 1, 1, -1, 0]
                                    
                                else:
                                    if j==x0:
                                        B1 += [0, 0, 0, 0, 0, 0]
                                        B2 += [0, 0, 0, 0, 0, 0]
                                        B3 += [0, 0, 0, 0, 0, 0]
                                        B4 += [0, 0, 0, 0, 0, 0]
                                        B5 += [0, 0, 0, 0, 0, 0]
                                        B6 += [0, 0, 0, 0, 0, 0]
                                    else:
                                        B1 += [0, 0, 0]
                                        B2 += [0, 0, 0]
                                        B3 += [0, 0, 0]
                                        B4 += [0, 0, 0]
                                        B5 += [0, 0, 0]
                                        B6 += [0, 0, 0]
                                        
                                        
                    else:  # j is in the last layer
                        if Vertices[j - 1][-2]!= 0:
                            B1 += [0, 0]
                            B2 += [0, 0]
                            B3 += [0, 0]
                            B4 += [0, 0]
                            B5 += [0, 0]
                            B6 += [0, 0]
                            # B1 += [0, 0 ,0]
                            # B2 += [0, 0 ,0]
                            # B3 += [0, 0 ,0]
                        else:
                            if (i, j) in Edges:
                                w=K(gen)*Eo[Edges.index((i,j))][1]
                                B1 += [0, 0, 0]
                                B2 += [0, w, 0]
                                B3 += [0, 0, 0]
                                B4 += [0, 0, 0]
                                B5 += [0, 0, 0]
                                B6 += [0, 0, 0]
                            else:
                                B1 += [0, 0, 0]
                                B2 += [0, 0, 0]
                                B3 += [0, 0, 0]
                                B4 += [0, 0, 0]
                                B5 += [0, 0, 0]
                                B6 += [0, 0, 0]
                    
                A += [B1, B2, B3, B4, B5, B6]
        
            #Now the other points 
                
            else:
                B1, B2, B3 = [], [], []
            
                for j in range(1, v + 1):
                    if j <= cmax(n, g - 1):  # i,j both not last layers
                        if j==x0 or j==y0: 
                            if Vertices[i - 1][-2]!= 0:
                                if j==fl(i, gen):
                                    w=K(gen)*Eo[Edges.index((j, i))][1]
                                    B1 += [0, 0, -w, 0, 0, 0]
                                    B2 += [0, 0, 0, 0, 0, 0]
                                    B3 += [0, 0, 0, 0, 0, 0]
                                
                                elif j==fr(i, gen):
                                    w=K(gen)*Eo[Edges.index((i, j))][1]
                                    B1 += [0, 0, 0,0,0,0]
                                    B2 += [0, 0, 0,0,0,0]
                                    B3 += [w, 0, 0,0,0,0]
                                
                                else:
                                    if (i, j) in Edges:
                                        w=K(gen)*Eo[Edges.index((i, j))][1]
                                        B1 += [0, 0, 0,0,0,0]
                                        B2 += [0, w, 0,0,0,0]
                                        B3 += [0, 0, 0,0,0,0]

                                    else:                                           
                                        B1 += [0, 0, 0, 0, 0, 0]
                                        B2 += [0, 0, 0, 0, 0, 0]
                                        B3 += [0, 0, 0, 0, 0, 0]
                                            
                            else:
                                if j == fl(i, gen):
                                    w = K(gen)*Eo[Edges.index((j, i))][1]
                                    B1 += [0, 0, -w,0,0,0]
                                    B2 += [0, 0, 0,0,0,0]
                                    B3 += [0, 0, 0,0,0,0]
                                
                                elif j == fr(i, gen):
                                    w = K(gen)*Eo[Edges.index((i, j))][1]
                                    B1 += [0, 0, 0,0,0,0]
                                    B2 += [0, 0, 0,0,0,0]
                                    B3 += [w, 0, 0,0,0,0]
                                    
                                else:   
                                    if (j, i) in Edges:
                                        w = K(gen)*Eo[Edges.index((j,i))][1]
                                        B1 += [0, 0, 0,0,0,0]
                                        B2 += [0, -w, 0,0,0,0]
                                        B3 += [0, 0, 0,0,0,0]

                                    elif (i,j) in Edges:
                                        w = K(gen)*Eo[Edges.index((j,i))][1]
                                        B1 += [0, 0, 0,0,0,0]
                                        B2 += [0, w, 0,0,0,0]
                                        B3 += [0, 0, 0,0,0,0]
                                    
                                    else:
                                        B1 += [0, 0, 0,0,0,0]
                                        B2 += [0, 0, 0,0,0,0]
                                        B3 += [0, 0, 0,0,0,0]
                               

                        else:
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
                            B1 += [0, 0]
                            B2 += [0, 0]
                            B3 += [0, 0]
                            # B1 += [0, 0 ,0]
                            # B2 += [0, 0 ,0]
                            # B3 += [0, 0 ,0]
                        else:
                            if (i, j) in Edges:
                                w=K(gen)*Eo[Edges.index((i,j))][1]
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
                if j==x0 or j==y0:
                    if (j, i) in Edges:
                        w = K(g)*Eo[Edges.index((j, i))][1]
                        B1 += [0, 0, 0,0,0,0]
                        B2 += [0, -w, 0,0,0,0]
                        B3 += [0, 0, 0,0,0,0]

                    else: 
                        B1 += [0, 0, 0,0,0,0]
                        B2 += [0, 0, 0,0,0,0]
                        B3 += [0, 0, 0,0,0,0]
                   
                else:
                    if j == fl(i, g):
                        w = K(g)*Eo[Edges.index((j, i))][1]
                        B1 += [0, -w]
                        B2 += [0, 0]
                        B3 += [0, 0]
                    elif j == fr(i, g):
                        w = K(g)*Eo[Edges.index((i, j))][1]
                        B1 += [0, 0]
                        B2 += [0, 0]
                        B3 += [w, 0]
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

                if j==x0 or j==y0:
                   B1 +=[0,0,0,0,0,0]
                   B2 +=[0,0,0,0,0,0]

                else:
                    if Vertices[j-1][-2]==0:
                        if j==fl(i, g):
                            w=K(g)*Eo[Edges.index((j, i))][1]
                            B1 += [0, 0, -w]
                            B2 += [0, 0,  0]
                        elif j==fr(i, g):
                            w=K(g)*Eo[Edges.index((i, j))][1]
                            B1 += [0, 0, 0]
                            B2 += [w, 0, 0]
                        else:
                            B1+=[0,0,0]
                            B2+=[0,0,0]
                    else:
                        if j==fl(i, g):
                            w=K(g)*Eo[Edges.index((j, i))][1]
                            B1 += [0, -w]
                            B2 += [0, 0]
                        elif j==fr(i,g):
                            w=K(g)*Eo[Edges.index((i, j))][1]
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


#Partition function part

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

# a=int(input("Enter number of generations"))
# E,Eo=Hyp(6,a)[1],Hyp(6,a)[2]
# G=nx.Graph(E)
# G.add_edges_from(E)
# nx.draw(G, with_labels=True)
# print(Eo)
# #print(Hyp(7,a)[0])
# plt.show()

#Building the Pfaffian matrix and calculate the eigenspectrum (May be Sturm sequnce method)

def A(b, I ,n , g):

    def K(a):
        return 1/np.tanh(b*I) # 1/np.tanh(b*I*Metric(a))
    
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
                        B1 += [0, 0]
                        B2 += [0, 0] #Asymmetric boundary condition
                        B3 += [0, 0]
                        
                        # B1 += [0, 0, 0, 0]
                        # B2 += [0, 0, 0, 0]
                        # B3 += [0, 0, 0, 0]
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
                    B1 += [0, -w]
                    B2 += [0, 0]
                    B3 += [0, 0]

                    # B1 += [0, 0, -w, 0]
                    # B2 += [0, 0, 0, 0]
                    # B3 += [0, 0, 0, 0]
                elif j == fr(i, g):
                    w = K(g)*Eo[Edges.index((i, j))][1]
                    B1 += [0, 0]
                    B2 += [0, 0]
                    B3 += [w, 0]

                    # B1 += [0, 0, 0, 0]
                    # B2 += [0, 0, 0, 0]
                    # B3 += [w, 0, 0, 0]
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
                            if j>cmax(n,g-1):
                                if Vertices[j-1][-2]==0:
                                   B1 += [0, 0, 0]
                                   B2 += [0, 0, 0]
                                   B3 += [0, 0, 0]
                                else:
                                    B1 += [0, 0]
                                    B2 += [0, 0]
                                    B3 += [0, 0]

                                    # B1 += [0, 0, 0, 0]
                                    # B2 += [0, 0, 0, 0]
                                    # B3 += [0, 0, 0, 0]
                            else:
                                B1 += [0, 0, 0]
                                B2 += [0, 0, 0]
                                B3 += [0, 0, 0]

            A += [B1, B2, B3]

        else:
            B1, B2 = [], [] #Is to be changed to triangular modification and put in periodic boundary condtition

            #B1,B2,B3,B4=[],[],[],[]

            for j in range(1, v + 1):
                if Vertices[j-1][-2]==0:
                    if j==fl(i, g):
                        w= K(g)*Eo[Edges.index((j, i))][1]
                        B1 += [0, 0, -w]
                        B2 += [0, 0,  0]

                        # B1 += [0, 0, -w]
                        # B2 += [0, 0,  0]
                        # B3 += [0, 0, 0]
                        # B4 += [0, 0, 0]
                    elif j==fr(i, g):
                        w= K(g)*Eo[Edges.index((i, j))][1]
                        B1 += [0, 0, 0]
                        B2 += [w, 0, 0]

                        # B1 += [0, 0, 0]
                        # B2 += [0, 0,  0]
                        # B3 += [w, 0, 0]
                        # B4 += [0, 0, 0]
                    else:
                        B1+=[0,0,0]
                        B2+=[0,0,0]

                        # B1+=[0,0,0]
                        # B2+=[0,0,0]
                        # B3+=[0,0,0]
                        # B4+=[0,0,0]
                else:
                    if j==fl(i, g):
                        w= K(g)*Eo[Edges.index((j, i))][1]
                        B1 += [0, -w]
                        B2 += [0, 0]

                        # B1+=[0,0,-w,0]
                        # B2+=[0,0,0,0]
                        # B3+=[0,0,0,0]
                        # B4+=[0,0,0,0]
                    elif j==fr(i,g):
                        w= K(g)*Eo[Edges.index((i, j))][1]
                        B1 += [0,  0]
                        B2 += [w, 0]

                        # B1+=[0,0,0,0]
                        # B2+=[0,0,0,0]
                        # B3+=[w,0,0,0]
                        # B4+=[0,0,0,0]
                    else:
                        if j==i:
                            B1+=[0,-1]
                            B2+=[1,0]

                            # B1+=[0,1,-1,0]
                            # B2+=[-1,0,1,1]
                            # B3+=[1,-1,0,0]
                            # B4+=[0,-1,0,0]
                        else:
                            if j>cmax(n,g-1):
                                B1 += [0,0]
                                B2 += [0,0]

                                # B1+=[0,0,0,0]
                                # B2+=[0,0,0,0]
                                # B3+=[0,0,0,0]
                                # B4+=[0,0,0,0]
                            else:
                                B1 += [0,0,0]
                                B2 += [0,0,0]

                                # B1 += [0,0,0]
                                # B2 += [0,0,0]
                                # B3 += [0,0,0]
                                # B4 += [0,0,0]

            A +=[B1, B2]
            #A += [B1, B2, B3, B4]

    return A

# E=ACorr(100,20,7,2,1,3) #Almost correct
# Eig=eig(np.array(E))[0]
# E1=np.imag(Eig)

# p=1
# for e in E1: 
#     if e>=0:
#         p=p*e 

# print(p)

#Hermitization of the matrix

# def Ah(b,I,n,g ,gc,d):
#     Ah=[]
#     L=A(b,I,n,g)
#     for i in range(len(L)):
#         Ah1=[]
#         for j in range(len(L)):
#             Ah1.append(complex(0,1)*L[i][j])
#         Ah.append(Ah1)
#     return Ah

#Plotting the eigenspectrum

n,g,I=7,2,1
n1,g1,I1=7,2,1

N=cmax(n,g)+Layers(n,g)[2]
N1=cmax(n1,g1)+Layers(n1,g1)[2]

B=np.linspace(10**-1,1.5,100)
E,LC,Corr=[],[],[]
E1,LC1,Corr1=[],[],[]
B0=[]

for b in B:
    Pfaff=np.array(ACorr(b,I,n,g,1,3))
    Pfaff1=np.array(ACorr(b,I1,n1,g1,1,5))
    e0=eig(Pfaff)[0]
    e01=eig(Pfaff1)[0]
    e0c=np.imag(e0)
    e0c1=np.imag(e01)   

    Z0_log=N*np.log(np.sinh(b*I))+cmax(n,g)*np.log(2)
    Z1_log=N1*np.log(np.sinh(b*I1))+cmax(n1,g1)*np.log(2)

    s0=0
    for e in e0c:
        if e>=0:
           s0+=np.log(e)
    
    s1=0
    for e1 in e0c1:
        if e1>=0:        
           s1+=np.log(e1)
 

    P=np.array(A(b,I,n,g))
    e=eig(P)[0]
    ec=np.imag(e)

    sum=0
    for h in ec:
        if h>=0:
           sum+=np.log(h)

    # s1=0
    # for e1 in e0c1: 
    #     if e1>=0:
    #        s1 +=np.log(e1)
    # Pf1.append(s1+N1*np.log(np.sinh(b*I1))+cmax(n1,g1)*np.log(2))
    
    E.append(e0c)
    E1.append(e0c1)
    LC.append(s0-sum)
    LC1.append(s1-sum)
    Corr.append(np.exp(s0-sum))
    Corr1.append(np.exp(s1-sum))
    
#plt.plot(B,E,'r+')
#plt.plot(B,E1,'g+')
# plt.plot(B,LC,'r+')
# plt.plot(B,LC1,'g+')
plt.plot(B,Corr,'r+')
plt.plot(B,Corr1,'g+')
plt.plot(B,B*0,'b')
plt.xlabel("Beta--->")
plt.ylabel("Correlation--->")
# plt.ylim([-10,10])
#plt.legend()
plt.show()

