import gurobipy as gp
from gurobipy import Model, GRB
# Vars
N = 3
D = 2
f = [
    [1, 1, 5, 13, 29, 62, 125, 253],
    [1, 1, 5, 13, 29, 52, 80, 253],
    [1, 1, 5, 13, 26, 52, 80, 253],
    [1, 1, 5, 13, 26, 38, 80, 253],
    [1, 1, 5, 13, 26, 38, 50, 253]
]
# Model
m = Model("MCT Quantum Circuit Design")
# Decision Vars
y = m.addVars(N+1, D+1, vtype = GRB.BINARY, name ="y")
t = m.addVars(N+1, D+1, vtype = GRB.BINARY, name ="t")
w = m.addVars(N+1, D+1, vtype = GRB.BINARY, name ="w")
# Objective Function
m.setObjective(gp.quicksum(gp.quicksum(f[N-j][j-1]*y[j, d] for j in range(1, N+1)) for d in range(1, D+1)), GRB.MINIMIZE) # objective function 1a, f[N-j][j-1] := N-j slack qubits and j-1 control qubit
# Constraints
for d in range(1, D+1):
    m.addConstrs(t[q, d] + w[q, d] <= 1 for q in range(1, N+1)) # constraints 1b
    m.addConstrs(w[q, d] <= gp.quicksum(t[r, d] for r in range(1,N+1)) for q in range(1, N+1)) # constraints 1d
m.addConstrs(gp.quicksum(t[q, d] for q in range(1, N+1)) <= 1 for d in range(1, D+1)) # constraint 1c
m.addConstrs(gp.quicksum(j*y[j, d] for j in range(1, N+1)) == (gp.quicksum(t[q, d] for q in range(1, N+1)) + gp.quicksum(w[q, d] for q in range(1, N+1))) for d in range(1, D+1)) # constraint 1e
m.addConstrs(gp.quicksum(y[j, d] for j in range(1, N+1)) <= 1 for d in range(1, D+1)) # constraint 1f
# Solve
def printSolution():
    for i in range(1, N+1):
        for j in range (1, D+1):
            print("Y_", i,",", j," : ", y[i, j].X)
            print("T_", i,",", j," : ", t[i, j].X)
            print("W_", i,",", j," : ", w[i, j].X)
m.optimize()
printSolution()
