import numpy as np
import math

# Parámetros del sistema de Lorenz
# sigma = 10
# rho = 28
# beta = 8/3

# a,b,c = 0.507309, 0.721377, 2.70869
a, b, c = 3.713830535072366246e+01, 4.626939085801102802e-01, 3.343102504306193623e+01

# Condiciones iniciales
x0, y0, z0 = 4.246, 4.728, 13.470

# Parámetros de la simulación
dt = 0.00075  # Paso de tiempo
num_steps = 10000  # Número de pasos

# Inicializar las listas para almacenar las soluciones
x = np.zeros(num_steps + 1)
y = np.zeros(num_steps + 1)
z = np.zeros(num_steps + 1)

x[0], y[0], z[0] = x0, y0, z0

# Variables para almacenar los máximos y mínimos intermedios
max_value = float('-inf')
min_value = float('inf')

# Simulación del sistema de Lorenz usando el método de Forward Euler
for i in range(num_steps):
    # dx = sigma * (y[i] - x[i])
    # dy = x[i] * (rho - z[i]) - y[i]
    # dz = x[i] * y[i] - beta * z[i]
    
    # dx = -y[i] - z[i]
    # dy = x[i] + a*y[i]
    # dz = b + z[i]*(x[i] - c)

    dx = a*(y[i] - x[i])
    dy = (c - a)*x[i] - x[i]*z[i] + c*y[i]
    dz = x[i]*y[i] - b*z[i]

    # Actualizar máximos y mínimos
    max_value = max(max_value, dx, dy, dz)
    min_value = min(min_value, dx, dy, dz)

    x[i + 1] = x[i] + dx * dt
    y[i + 1] = y[i] + dy * dt
    z[i + 1] = z[i] + dz * dt

# Calcular el número de bits necesarios


def calculate_bits_needed(value):
    if value == 0:
        return 1
    return math.floor(math.log2(abs(value))) + 1


max_abs_value = max(abs(max_value), abs(min_value))
total_bits_needed = calculate_bits_needed(max_abs_value)

# Determinar la parte fraccionaria basada en 32 bits totales
frac_bits = 32 - total_bits_needed

print("Máximo valor absoluto en el sistema:", max_abs_value)
print("Número de bits necesarios en la parte entera:", total_bits_needed)
print("Número de bits necesarios en la parte fraccionaria:", frac_bits)
print(f"Total de bits necesarios: {total_bits_needed}")
