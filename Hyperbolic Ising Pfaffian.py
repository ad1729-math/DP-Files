import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from scipy.linalg import fractional_matrix_power

# Create a sample graph (you can replace this with your own graph creation)
# G = nx.Graph()
#G.add_edges_from([(1, 2), (1, 3), (3,2), (3, 4),(4,5),(4,6),(5,6)])
# nx.draw(G, with_labels=True, node_color='skyblue', node_size=500, font_weight='bold')

# # Get the adjacency matrix
# A= nx.to_numpy_array(G)

# print("Adjacency Matrix:")
# print(A)


def pfaffian(matrix):
    n = matrix.shape[0]
    if n % 2 != 0:
        raise ValueError("Pfaffian exists only for even-sized matrices.")
    return np.linalg.det(matrix.reshape(n//2, 2, n//2, 2).transpose(0, 2, 1, 3))

# Create a sample graph (replace this with your graph creation)
G = nx.Graph()
G.add_edges_from([(1, 2), (1, 3), (3,2), (3, 4),(4,5),(4,6),(5,6)])

# Get the adjacency matrix
adj_matrix = nx.to_numpy_array(G)

# Compute the Pfaffian using the adjacency matrix
pfaffian_value = pfaffian(adj_matrix)

print("Pfaffian Matrix:")
print(pfaffian_value)
