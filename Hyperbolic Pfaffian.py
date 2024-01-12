import math as m 
import numpy as np 
import matplotlib.pyplot as plt 
from hypertiling import HyperbolicTiling, HyperbolicGraph
from hypertiling.graphics.plot import quick_plot
import networkx as nx

#2-point correlation using Dimer model
# def Am(b,J,J1):
#     v,v1=np.tanh(J*b),np.tanh(J1*b)

# n,k,g=7,3,5
# Til=HyperbolicTiling(n,k,g)
# quick_plot(Til)

E=[]
for i in range(1,8):
    E.append((i,i%7+1))
    E.append((i,i+7))
    E.append((i+7,i+14))
    E.append((i+14,i+21))
    E.append((i+21,i+28))
    E.append((i+28,i%7+8))




G=nx.Graph(E)
G.add_edges_from(E)
nx.draw(G, with_labels=True)
plt.show()