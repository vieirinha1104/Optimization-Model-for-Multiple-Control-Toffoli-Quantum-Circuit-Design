boolean_function = {}
n = 3
file_name = 'input.txt'
data = []
with open(file_name, 'r') as file:
    for row in file:
        row = row.strip()
        if row:
            data.append(row)

for row in data:
    column = row.split('\t')
    input = 0
    output = 0
    aux = 1
    for i in range(0,n):
        input += aux*int(column[0][n-1-i])
        output += aux*int(column[1][n-1-i])
        aux *= 2
    boolean_function[input] = output
    
print(boolean_function)