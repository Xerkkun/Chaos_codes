import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

# Definir los puntos de interpolación en el rango [0, pi/2]
x_interp = np.array([0, np.pi/6, np.pi/4, np.pi/3, np.pi/2])
y_interp = np.sin(x_interp)

# Crear el spline cúbico
cs = CubicSpline(x_interp, y_interp)

# Definir un rango de valores para evaluar la función original y la aproximación
x_values = np.linspace(0, np.pi/2, 100)
y_true = np.sin(x_values)
y_spline = cs(x_values)

# Calcular el error entre la función original y la aproximación
error = y_true - y_spline

# Obtener los coeficientes del spline cúbico
a = cs.c[3, :]  # Coeficientes a_i
b = cs.c[2, :]  # Coeficientes b_i
c = cs.c[1, :]  # Coeficientes c_i
d = cs.c[0, :]  # Coeficientes d_i

# Imprimir los coeficientes para cada subintervalo
print("Coeficientes del spline cúbico para cada subintervalo:")
for i in range(len(a)):
    print(f"Subintervalo [{x_interp[i]:.2f}, {x_interp[i+1]:.2f}]:")
    print(f"S(x) = {a[i]} + {b[i]}*(x - {x_interp[i]}) + {c[i]}*(x - {x_interp[i]})^2 + {d[i]}*(x - {x_interp[i]})^3")
    print()

# Graficar la función original, la aproximación spline y el error

# Gráfico de la función original y la aproximación spline
plt.figure(figsize=(8, 6))
plt.plot(x_values, y_true, label='Original function $\sin(x)$', color='blue')
plt.plot(x_values, y_spline, '--', label='Cubic spline approximation', color='red')
plt.scatter(x_interp, y_interp, color='green', label='Interpolation points')
plt.xlabel('x (radians)')
plt.ylabel('$\sin(x)$')
plt.legend()
plt.grid(True)
plt.savefig('spline_approximation.png', format='png')
plt.show()

# Gráfico del error
plt.figure(figsize=(8, 6))
plt.plot(x_values, error, label='Error (Original - Approximation)', color='purple')
plt.axhline(0, color='black', linestyle='--', linewidth=1)
plt.xlabel('x (radians)')
plt.ylabel('Error')
plt.legend()
plt.grid(True)
plt.savefig('spline_error.png', format='png')
plt.show()
