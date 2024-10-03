import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import itertools
import gurobipy as gp
from gurobipy import Model, GRB

# Input Variables
n, d = (3, 2)  # N: Number of qubits, d: Maximum number of gates of the circuit

# Aux vars
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
omega_input, omega_output = ([], [])

# Dictonary for indexing the Toffoli gates into their costs
costs_dict = {} # Example: 0: ['+', 'o', 'o'] -> costs_dict[0] =  5 

# Dictionaries for indexing qubit states 
qubits_dict = {} # Example: {0: [0, 0, 0], 1: [0, 0, 1], 2: [0, 1, 0], 3: [0, 1, 1], ... }

# Dictionary for indexing all possible Toffoli gates
toffoli_dict = {} # 'o': control qubit, '+': target qubit, '-': slack qubit

# Dictionary for finding the position of the flipped qubit by arc key
flipQ_dict = {} 

# Functions
# This function set all the elements in the toffoli_dict, O(n*(3^n)) 
def setToffoliDict(k):
    if(k == n):
        v = [0,0,0] # v[0]: # of target qubits, v[1]: # of control qubits, v[2]: # of slack qubits
        for x in brr:
            if(x == '+'):
                v[0] += 1
            elif(x == 'o'):
                v[1] += 1
            else:
                v[2] += 1
        if(v[0] != 1):
            return
        m = len(toffoli_dict)
        toffoli_dict[m] = brr.copy()
        costs_dict[m] = costs[v[2]][v[1]]
    else:
        for i in ['+','o','-']:
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
        if(f[i] == '+'):
            target = i
        elif(f[i] == 'o'):
            aux *= q_input[i]
    if(aux == 1): # flips
        q_output[target] = 1 - q_input[target]
        flip = target + 1
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
            node = (j, i)
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
                    flipQ_dict[(node, vertex)] = new_state[1]
    return G

def omegaPartition():
    boolean_function = {}
    file_name = 'ex2.txt' # ex1: paper's example 1 instance, regular instance (d=3), ex2: paper's example 2 instance, dont care instance (d=2)
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
        omega_input.append(aux)
    for x in omega_input:
        output = boolean_function[x[0]]
        new_output = []
        for s in output:
            if(s != '-'):
                new_output.append(s)
            else:
                new_output.append(['0','1'])
        aux = itertools.product(*new_output)
        omega_output.append([int(''.join(map(str, i)), 2) for i in aux])

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

# Prints for debug:
def printEverything(G):
    print("Nodes: ",G.nodes())
    print("\n")
    print("Edges: ", G.edges())
    print("\n")
    print("Omega Input: ",omega_input)
    print("\n")
    print("Omega Output: ", omega_output)
    print("\n")
    print("Toffoli Dict: ", toffoli_dict)
    print("\n")
    print("Costs Dict: ", costs_dict)
    print("\n")
    print("Qubits_Dict: ", qubits_dict)
    print("\n") 
    print("Flip Qubits Dict: ", flipQ_dict)

# find_zero_bit_positions(a, n) is Q _ {sigma(a),0}
def find_zero_bit_positions(number, N):
    binary_representation = format(number, '0{}b'.format(N))
    zero_positions = [i+1 for i, bit in enumerate(binary_representation) if bit == '0']
    return zero_positions

# Print the Circuit
def printSolution(m, n, d, t, w):
    if (m.status == GRB.OPTIMAL):
        print("Circuit Cost: ", m.ObjVal)
        print("Circuit Design:\n")
        for j in range (1,d+1):
            print("(", end = '')
            for i in range (1,n+1):
                if(t[i, j].X > 0.5):
                    print("+", end = '')
                elif(w[i, j].X > 0.5):
                    print("o", end = '')
                else:
                    print("-", end = '')
                if(i != n):
                    print(", ", end = '')
            print(")", end = '')
            print("\n")
    else:
        print("Model is infeasible.")
        m.computeIIS()
        print("IIS computed. Constraints that are conflicting:")
        for c in m.getConstrs():
            if c.IISConstr:
                print(c.constrName)

# Build the graph G_k for a given commodity k
def set_Gk(G, D, k):
    H = G.copy()
    H.add_node(('s', 0), layer = 0)
    H.add_node(('t', D+2), layer = D+2)
    # Source Arcs:
    for v in omega_input[k-1]:
        H.add_edge(('s', 0), (v, 1), weight = 0, layer = 0, keep = False, flip = False)
    for u in omega_output[k-1]:
        H.add_edge((u, D+1), ('t', D+2), weight = 0, layer = D+1, keep = False, flip = False)
    return H

# k_arcs is the set of all (e,k), where e is an arc of G_k 
def set_Karcs(G, D, k):
    k_arcs = []
    for i in range(1, k+1):
        H = set_Gk(G, D, i)
        for e in H.edges():
            k_arcs.append((e, i))
        # plotGraph(H)
    return k_arcs

# Return the set of flip arcs and keep arcs for a given G_k
def set_keep_and_flip_arcs(H):
    keep_arcs, flip_arcs = ([],[])
    for e in H.edges():
        u, v = (e[0], e[1])
        if(H[u][v]["keep"] == True):
            keep_arcs.append(e)
        if(H[u][v]["flip"] == True):
            flip_arcs.append(e)
    return flip_arcs, keep_arcs

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
    # Constraints 1b-1f
    m.addConstrs(t[q, d] + w[q, d] <= 1 for q in range(1, N+1) for d in range(1, D+1)) # constraints 1b
    m.addConstrs(gp.quicksum(t[q, d] for q in range(1, N+1)) <= 1 for d in range(1, D+1)) # constraint 1c
    m.addConstrs(w[q, d] <= gp.quicksum(t[r, d] for r in range(1,N+1)) for q in range(1, N+1) for d in range(1, D+1)) # constraints 1d
    m.addConstrs(gp.quicksum(j*y[j, d] for j in range(1, N+1)) == (gp.quicksum(t[q, d] for q in range(1, N+1)) + gp.quicksum(w[q, d] for q in range(1, N+1))) for d in range(1, D+1)) # constraint 1e
    m.addConstrs(gp.quicksum(y[j, d] for j in range(1, N+1)) <= 1 for d in range(1, D+1)) # constraint 1f
    # Commodities
    k = len(omega_input) 
    k_arcs = set_Karcs(G, D, k)
    # Flow Decision Vars
    x = m.addVars(k_arcs, vtype = GRB.BINARY, name = "x")
    # Flow Constraints
    for i in range(1, k+1):
        H = set_Gk(G, D, i)
        # Constraints 1g
        for v in H.nodes():
            a_in = H.in_edges(v) # sigma_{-}(v)
            a_out = H.out_edges(v) # sigma_{+}(v)
            if(v[0] == 's'):
                m.addConstr(gp.quicksum(x[a, i] for a in a_out) - gp.quicksum(x[a, i] for a in a_in) == len(omega_input[i-1]))
            elif(v[0] == 't'):
                m.addConstr(gp.quicksum(x[a, i] for a in a_out) - gp.quicksum(x[a, i] for a in a_in) == -len(omega_input[i-1]))
            else:
                m.addConstr(gp.quicksum(x[a, i] for a in a_out) - gp.quicksum(x[a, i] for a in a_in) == 0)
        # Constraints 2a, 2b, 2c
        flip_arcs, keep_arcs = set_keep_and_flip_arcs(H)
        a, b = (len(flip_arcs), len(keep_arcs))
        for j in range(0, a):
            m.addConstr(x[flip_arcs[j], i] <= t[flipQ_dict[flip_arcs[j]], flip_arcs[j][0][1]]) # 2a
            m.addConstrs(x[flip_arcs[j], i] <= (1 - w[p, flip_arcs[j][0][1]]) for p in find_zero_bit_positions(flip_arcs[j][1][0], n)) # 2b
        for j in range(0, b):
            m.addConstr(x[keep_arcs[j], i] <= (1 - gp.quicksum(t[q, keep_arcs[j][0][1]] for q in range(1,N+1)) + gp.quicksum(w[q, keep_arcs[j][0][1]] for q in find_zero_bit_positions(keep_arcs[j][1][0], N)))) # 2c
    m.optimize()
    printSolution(m, n, d, t, w)

# Main
setQubitsDict(0)
setToffoliDict(0)
G = multilayered_graph(n,d)
omegaPartition()
flowModel(G, n, d, costs)
# printEverything(G)



