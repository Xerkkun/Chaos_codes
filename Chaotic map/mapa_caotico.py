import numpy as np
import matplotlib.pyplot as plt

# Parámetros
alpha = 10.75
beta = 2.45
n_iter = 10000  # Número de iteraciones

# Condiciones iniciales
x0 = 0.412
y0 = 0.784

# Inicializar las secuencias
x = np.zeros(n_iter)
y = np.zeros(n_iter)

# Condiciones iniciales
x[0] = x0
y[0] = y0

# Iteraciones
for n in range(1, n_iter):
    x[n] = (alpha * np.sin(2 * np.pi * y[n-1])) / (x[n-1] * np.sin(1 - 2 * np.pi * y[n-1]))
    x[n] = x[n] % 1  # Aplicar la función mod 1
    y[n] = (beta * np.cos(2 * np.pi * x[n-1])) / (y[n-1] * np.sin(1 - 2 * np.pi * x[n-1]))
    y[n] = y[n] % 1  # Aplicar la función mod 1

# Graficar el diagrama de fases
plt.figure(figsize=(8, 8))  # Asegura una figura cuadrada
plt.plot(x, y, 'o', markersize=0.5, color='blue')
plt.title("Diagrama de Fases $x_n$ vs $y_n$")
plt.xlabel("$x_n$")
plt.ylabel("$y_n$")
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.grid(True)
plt.gca().set_aspect('equal', adjustable='box')  # Mantener aspecto cuadrado
plt.show()
