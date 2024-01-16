import math as m 
import numpy as np 
import matplotlib.pyplot as plt 
from hypertiling import HyperbolicTiling, HyperbolicGraph
from hypertiling.graphics.plot import quick_plot
import networkx as nx
# import mplcursors

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

for i in range(15, 36):
    if 1<=(i-28)<=7:
        E.append((i,i+21))
        E.append((i+21,i+42))
        E.append((i+42,i+63))
        E.append((i+63,(i-28)%7+36))
    else:
        E.append((i,i+21))
        E.append((i+21,i+42))
        E.append((i+42,i+63))
        E.append((i+63,i+84))
        E.append((i+84,i+28))


G=nx.Graph(E)
G.add_edges_from(E)
nx.draw(G, with_labels=True)
plt.show()

# pos = nx.spring_layout(G, seed=42)
# for node in pos:
#     pos[node][1] *= 2

# # Create a plot
# fig, ax = plt.subplots()

# # Draw the graph
# nx.draw(G, pos=pos, with_labels=True, ax=ax)
# plt.show()
# # Data structure to store node dragging state
# # dragging = None


# def on_press(event):
#     global dragging
#     if event.inaxes is not None:
#         for node, (x, y) in pos.items():
#             if x - 0.1 <= event.xdata <= x + 0.1 and y - 0.1 <= event.ydata <= y + 0.1:
#                 dragging = node
#                 break


# def on_release(event):
#     global dragging
#     dragging = None


# def on_motion(event):
#     global dragging
#     if dragging is not None:
#         pos[dragging] = (event.xdata, event.ydata)
#         ax.clear()
#         nx.draw(G, pos=pos, with_labels=True, ax=ax)
#         fig.canvas.draw()


# fig.canvas.mpl_connect('button_press_event', on_press)
# fig.canvas.mpl_connect('button_release_event', on_release)
# fig.canvas.mpl_connect('motion_notify_event', on_motion)

# plt.show()