import cmath
import numpy as np

x = np.zeros(2, dtype=complex)
# Coeficientes
a = 1
b = -2
c = 2  # Cambiado para obtener soluciones complejas

# Calcular el discriminante
D = 16

# Calcular las soluciones
x[0] = (-b - cmath.sqrt(D)) / (2*a)
x[1] = (-b + cmath.sqrt(D)) / (2*a)

print("Las soluciones son:", x)
print(x.real)  # Imprime: 3.0
print(x.imag)  # Imprime: 4.0