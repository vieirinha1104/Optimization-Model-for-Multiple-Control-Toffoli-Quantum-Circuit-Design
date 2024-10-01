import itertools 

boolean_function = {}
n = 3
file_name = 'test.txt'
data = []
omega_input = []
omega_output = []
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
    omega_output.append([''.join(i) for i in aux])

print("Omega Input:", omega_input)
print("Omega Output:", omega_output)