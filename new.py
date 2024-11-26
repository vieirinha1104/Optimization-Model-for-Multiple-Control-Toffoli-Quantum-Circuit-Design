import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import itertools
import gurobipy as gp
from gurobipy import Model, GRB

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

# Arrays for the Omega Partition
omega_input = []
omega_output = []

# Dictonary for indexing the Toffoli gates into their costs
costs_dict = {}
# Example: 0: ['o', 'x', 'x'] -> costs_dict[0] =  5 

# Dictionaries for indexing qubit states and Toffoli gates
qubits_dict = {}  # Dictionary for all possible qubit arrays
# Example: {0: [0, 0, 0], 1: [0, 0, 1], 2: [0, 1, 0], 3: [0, 1, 1], ... }

toffoli_dict = {}  # Dictionary for indexing all possible Toffoli gates
# 'x': control qubit, 'o': target qubit, '-': slack qubit

q_dict = {} # q_dict[key] = value, value is the position of the flipped qubit by arc key
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
    flip = -1
    for i in range(0,n):
        if(f[i] == 'o'):
            target = i
        elif(f[i] == 'x'):
            aux *= q_input[i]
    if(aux == 1): # flips
        q_output[target] = 1 - q_input[target]
        flip = target+1
    # but q_output is the qubits, not the index
    prod = 1
    index = 0
    for i in range(0,n):
        index += prod*q_output[n-i-1]
        prod *= 2
    return [index,flip]

# Build the Graph, O(d*n*2^n))
def multilayered_graph(n,d):
    G = nx.DiGraph()
    for i in range(1, d+2):
        for j in range(0, 2**n):
            node = (j,i)
            G.add_node(node, layer = i)
    for i in range(1, d+1):
        nodes_per_layer = [node for node, data in G.nodes(data=True) if data.get('layer') == i]
        for node in nodes_per_layer:
            for gate in toffoli_dict.keys():
                new_state = transition(node[0], gate)
                vertex = (new_state[0],i+1)
                # print(node, vertex, costs_dict[gate])
                if(new_state[1] == -1):
                    G.add_edge(node, vertex, weight = costs_dict[gate], layer = i, flip = False, keep = True)
                else:
                    G.add_edge(node, vertex, weight = costs_dict[gate], layer = i, flip = True, keep = False)
                    q_dict[(node, vertex)] = new_state[1]

    return G

def omegaPartition():
    boolean_function = {}
    file_name = 'test.txt'
    data = []
    with open(file_name, 'r') as file:
        for row in file:
            row = row.strip()
            if row:
                data.append(row)
    index = 0
    for row in data:
        column = row.split(' ')
        boolean_function[index] = column[1]
        index += 1
    values = []
    for i in range(0,2**n):
        values.append(boolean_function[i])
    check = [False]*(2**n)
    for i in range(0, 2**n):
        if(check[i] == True):
            continue
        s1 = values[i]
        aux = [i]
        check[i] = True
        for j in range(i+1,2**n):
            if(j >= 2**n):
                break
            s2 = values[j]
            if(s1 == s2):
                aux.append(j)
                check[j] = True
        omega_input.append((aux, 1)) 

    for x in omega_input:
        output = boolean_function[x[0][0]]
        new_output = []
        for s in output:
            if(s != '-'):
                new_output.append(s)
            else:
                new_output.append(['0','1'])
        aux = itertools.product(*new_output)
        omega_output.append(([int(''.join(map(str, i)), 2) for i in aux], d+1))

# Plot Graph
def plotGraph(G):
    # Layout
    pos = {}
    layer_spacing = 2  # Vertical spacing between layers
    node_spacing = 1.5  # Horizontal spacing between nodes

    # Calculate pos for plot
    for layer in range(0,d+3):
        nodes_in_layer = [node for node in G.nodes() if node[1] == layer]
        nodes_in_layer.sort(key=lambda x: x[0])  # Sort nodes by their index in qubits_dict
        for idx, node in enumerate(nodes_in_layer):
            pos[node] = (idx * node_spacing, -layer * layer_spacing)

    # Plot the Graph
    plt.figure(figsize=(12, 8))
    nx.draw(G, pos, with_labels=True, node_size=700, node_color='lightblue', arrows=True, arrowstyle='-|>', arrowsize=20)

    plt.axis("equal")
    plt.title('Multipartite Graph')
    plt.show()

# Prints:
def printEverything(G):
    print("Nodes: ",G.nodes())
    print("Edges: ", G.edges())
    print(omega_input)
    print(omega_output)
    print(toffoli_dict)
    print(costs_dict)
    print(qubits_dict)

# find_zero_bit_positions(a, n) is Q _ {sigma(a),0}
def find_zero_bit_positions(number, N):
    binary_representation = format(number, '0{}b'.format(N))
    zero_positions = [i for i, bit in enumerate(binary_representation) if bit == '0']
    return zero_positions

# Gurobi Code
def flowModel(G, N, D, f):
    # Model
    m = Model("MCT_Quantum_Circuit_Design")
    # Decision Vars
    y = m.addVars(N+1, D+1, vtype = GRB.BINARY, name ="y")
    t = m.addVars(N+1, D+1, vtype = GRB.BINARY, name ="t")
    w = m.addVars(N+1, D+1, vtype = GRB.BINARY, name ="w")
    # Objective Function
    m.setObjective(gp.quicksum(gp.quicksum(f[N-j][j-1]*y[j, d] for j in range(1, N+1)) for d in range(1, D+1)), GRB.MINIMIZE) # objective function 1a, f[N-j][j-1] := N-j slack qubits and j-1 control qubit
    # Constraints
    m.addConstrs(t[q, d] + w[q, d] <= 1 for q in range(1, N+1) for d in range(1, D+1)) # constraints 1b
    m.addConstrs(gp.quicksum(t[q, d] for q in range(1, N+1)) <= 1 for d in range(1, D+1)) # constraint 1c
    m.addConstrs(w[q, d] <= gp.quicksum(t[r, d] for r in range(1,N+1)) for q in range(1, N+1) for d in range(1, D+1)) # constraints 1d
    m.addConstrs(gp.quicksum(j*y[j, d] for j in range(1, N+1)) == (gp.quicksum(t[q, d] for q in range(1, N+1)) + gp.quicksum(w[q, d] for q in range(1, N+1))) for d in range(1, D+1)) # constraint 1e
    m.addConstrs(gp.quicksum(y[j, d] for j in range(1, N+1)) <= 1 for d in range(1, D+1)) # constraint 1f
    k = len(omega_input)

    # Flow Decision Vars
    x = m.addVars(H.edges, k, vtype = GRB.BINARY, name = "x")
    # Constraints 1g
    for i in range(1, k+1):
            m.addConstrs(gp.quicksum(x[v, u, i] for u in H.successors(v)) - gp.quicksum(x[u, v, i] for u in H.predecessors(v)) == 0 for v in H.nodes() if v not in omega_input[i-1] and v not in omega_output[i-1])
            m.addConstrs(gp.quicksum(x[v, u, i] for v in omega_input[i-1] for u in H.successors(v) if G.nodes[v]['layer'] == 1) == 1)
            m.addConstrs(gp.quicksum(x[u, v, i] for v in omega_output[i-1] for u in H.predecessors(v) if G.nodes[v]['layer'] == d+1) == 1)

    # Flow Constraints
    for i in range(1, k+1):
        H = G.copy()
        H.add_node(('s', 0), layer = 0)
        H.add_node(('t', D+2), layer = D+2)
        # Source Arcs:
        for v in omega_input[i-1]:
            H.add_edge(('s', 0), (v, 1), weight = 0, layer = 0, keep = False, flip = False)
        for u in omega_output[i-1]:
            H.add_edge((u, D+1), ('t', D+2), weight = 0, layer = D+1, keep = False, flip = False)
        # Constraints 1g
        m.addConstrs(gp.quicksum(x[(v, u), i] for u in H.successors) - gp.quicksum(x[(u, v), i] for u in H.predecessors) == 0 for v in H.nodes() if v[0] != 's' and v[0] != 't')

        for v in H.nodes():
            a_in = H.in_edges(v)
            a_out = H.out_edges(v)
            if(v[0] == 's'):
                m.addConstr(gp.quicksum(x[a, i] for a in a_in) - gp.quicksum(x[a, i] for a in a_out) == -len(omega_input[i-1]))
            elif(v[0] == 't'):
                m.addConstr(gp.quicksum(x[a, i] for a in a_in) - gp.quicksum(x[a, i] for a in a_out) == len(omega_input[i-1]))
            else:
                m.addConstr(gp.quicksum(x[a, i] for a in a_in) - gp.quicksum(x[a, i] for a in a_out) == 0)
        # Constraints 2a, 2b, 2c
        keep_arcs = [(u, v) for u, v, d in H.edges(data=True) if d.get('keep') is True]
        flip_arcs =  [(u, v) for u, v, d in H.edges(data=True) if d.get('flip') is True]
        for j in range(0,len(flip_arcs)):
            a = flip_arcs[j]
            b = q_dict[flip_arcs[j]]
            c = flip_arcs[j][0][1]
            m.addConstr(x[flip_arcs[j], i] <= t[q_dict[flip_arcs[j]], flip_arcs[j][0][1]]) # 2a
            m.addConstrs(x[flip_arcs[j], i] <= (1 - w[p, flip_arcs[j][0][1]]) for p in find_zero_bit_positions(flip_arcs[j][1][0], n)) # 2b
        for j in range(0,len(keep_arcs)):
            m.addConstr(x[keep_arcs[j], i] <= (1 - gp.quicksum(t[q, keep_arcs[j][0][1]] for q in range(1,n+1)) + gp.quicksum(w[q, keep_arcs[j][0][1]] for q in find_zero_bit_positions(keep_arcs[j][1][0], n)))) # 2c
    m.optimize()


# Main
setQubitsDict(0)
setToffoliDict(0)
G = multilayered_graph(n,d)
omegaPartition()
# flowModel(G, n, d, costs)
printEverything(G)



