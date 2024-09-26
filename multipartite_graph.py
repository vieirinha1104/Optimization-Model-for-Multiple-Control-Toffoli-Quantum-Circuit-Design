import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

# Variables
n = 3  # Number of qubits
d = 2  # Maximum number of gates of the circuit

# Initialize arrays
arr = [0] * n  # Aux array for storing qubit states
brr = ['-'] * n  # Aux array for storring Toffoli gates

# Table for the Quantum Cost of the MCT Gates
costs = [
    [1, 1, 5, 13, 29, 62, 125, 253],
    [1, 1, 5, 13, 29, 52, 80, 253],
    [1, 1, 5, 13, 26, 52, 80, 253],
    [1, 1, 5, 13, 26, 38, 80, 253],
    [1, 1, 5, 13, 26, 38, 50, 253]
]

# Dictonary for indexing the Toffoli gates into their costs
costs_dict = {}
# Example: 0: ['o', 'x', 'x'] -> costs_dict[0] =  5 

# Dictionaries for indexing qubit states and Toffoli gates
qubits_dict = {}  # Dictionary for all possible qubit arrays
# Example: {0: [0, 0, 0], 1: [0, 0, 1], 2: [0, 1, 0], 3: [0, 1, 1], ... }

toffoli_dict = {}  # Dictionary for indexing all possible Toffoli gates
# 'x': control qubit, 'o': target qubit, '-': slack qubit

# Functions

# This function set all the elements in the toffoli_dict, O(n*(3^n)) 
def setToffoliDict(k):
    if(k == n):
        v = [0,0,0] # v[0]: # of target qubits, v[1]: # of control qubits, v[2]: # of slack qubits
        for x in brr:
            if(x == 'o'):
                v[0] += 1
            elif(x == 'x'):
                v[1] += 1
            else:
                v[2] += 1
        if(v[0] != 1):
            return
        m = len(toffoli_dict)
        toffoli_dict[m] = brr.copy()
        costs_dict[m] = costs[v[2]][v[1]]
    else:
        for i in ['o','x','-']:
            aux = brr[k]
            brr[k] = i
            setToffoliDict(k+1)
            brr[k] = aux

# This function set all the elements in the qubits_dict, O(2^n) in all cases
def setQubitsDict(k):
    if(k == n):
        m = len(qubits_dict)
        qubits_dict[m] = arr.copy()
    else:
        for i in range(0,2):
            aux = arr[k]
            arr[k] = i
            setQubitsDict(k+1)
            arr[k] = aux


# This function returns the qubits output of f(qi), where qi: qubits input and f: toffoli gate, O(n)
def transition(index_q, index_f):
    #index_q is the index of q in qubits_dict
    #index_f is the index of f in toffoli_dict
    f = toffoli_dict[index_f]
    q_input = qubits_dict[index_q]
    q_output = q_input.copy()
    aux = 1
    target = -1
    for i in range(0,n):
        if(f[i] == 'o'):
            target = i
        elif(f[i] == 'x'):
            aux *= q_input[i]
    if(aux == 1):
        q_output[target] = 1 - q_input[target]
    # but q_output is the qubits, not the index
    prod = 1
    index = 0
    for i in range(0,n):
        index += prod*q_output[n-i-1]
        prod *= 2
    return index

# Build the Graph, O(d*n*2^n))
def multilayered_graph(n,d):
    G = nx.DiGraph()
    for i in range(0, d+1):
        for j in range(0, 2**n):
            node = (j,i)
            G.add_node(node, layer = i)
    for i in range(0, d):
        nodes_per_layer = [node for node, data in G.nodes(data=True) if data.get('layer') == i]
        for node in nodes_per_layer:
            for gate in toffoli_dict.keys():
                vertex = (transition(node[0], gate),i+1)
                # print(node, vertex, costs_dict[gate])
                G.add_edge(node, vertex, weight = costs_dict[gate], layer = i)
    return G

# Main

# print(toffoli_dict)
# print(costs_dict)
# print(qubits_dict)
setQubitsDict(0)
setToffoliDict(0)
G = multilayered_graph(n,d)
# print("Nodes: ",G.nodes())
# print("\n")
# print("Edges: ", G.edges())

# Plot the Graph
pos = {}
layer_spacing = 2  # Vertical spacing between layers
node_spacing = 1.5  # Horizontal spacing between nodes

# Calculate positions
for layer in range(d+1):
    nodes_in_layer = [node for node in G.nodes() if node[1] == layer]
    nodes_in_layer.sort(key=lambda x: x[0])  # Sort nodes by their index in qubits_dict
    for idx, node in enumerate(nodes_in_layer):
        pos[node] = (idx * node_spacing, -layer * layer_spacing)

# Plot 
plt.figure(figsize=(12, 8))
nx.draw(G, pos, with_labels=True, node_size=700, node_color='lightblue', arrows=True, arrowstyle='-|>', arrowsize=20)

plt.axis("equal")
plt.title('Multipartite Graph')
plt.show()

for i in range(0,2**n):
    path = nx.single_source_shortest_path(G,(i,0))
    target = path.keys()
    print(qubits_dict[i]," : ", end = ' ')
    for e in target:
        if(e[1] == d):
            print(qubits_dict[e[0]], end = ' ')
    print("\n")