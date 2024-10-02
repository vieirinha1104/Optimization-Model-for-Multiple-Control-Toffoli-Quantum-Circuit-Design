def find_zero_bit_positions(number, N):
    # Converte o número para binário, removendo o prefixo '0b' e preenchendo com zeros à esquerda
    binary_representation = format(number, '0{}b'.format(N))
    
    # Encontra as posições dos bits que são zero
    zero_positions = [i for i, bit in enumerate(binary_representation) if bit == '0']
    
    return zero_positions

# Exemplo de uso
number = 10  # Número em base 10
N = 4       # Número de bits significativos
zero_positions = find_zero_bit_positions(number, N)

print("Posições dos bits que são zero:", zero_positions)