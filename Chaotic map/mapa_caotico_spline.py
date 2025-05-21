import numpy as np
import matplotlib.pyplot as plt

# Parámetros
alpha = 10.75
beta = 2.45
n_iter = 10000  # Número de iteraciones

# Condiciones iniciales
x0 = 0.412
y0 = 0.784

# Definición de los polinomios cúbicos para la aproximación de sin(x)
def spline_sin(x):
    pi_6 = np.pi / 6
    pi_4 = np.pi / 4
    pi_3 = np.pi / 3
    pi_2 = np.pi / 2

    if 0 <= x < pi_6:
        return 0 + 1.0049 * x + -0.0201 * x**2 + -0.1440 * x**3
    elif pi_6 <= x < pi_4:
        x_adj = x - pi_6
        return 0.5 + 0.8654 * x_adj + -0.2463 * x_adj**2 + -0.1440 * x_adj**3
    elif pi_4 <= x < pi_3:
        x_adj = x - pi_4
        return 0.7071 + 0.7069 * x_adj + -0.3594 * x_adj**2 + -0.0837 * x_adj**3
    elif pi_3 <= x <= pi_2:
        x_adj = x - pi_3
        return 0.8660 + 0.5014 * x_adj + -0.4252 * x_adj**2 + -0.0837 * x_adj**3
    else:
        raise ValueError("Input x is out of the expected range for this spline approximation.")

# Extensión de la aproximación para el intervalo [0, 2*pi]
def sin_approx(x):
    x = x % (2 * np.pi)
    if 0 <= x < np.pi / 2:
        return spline_sin(x)
    elif np.pi / 2 <= x < np.pi:
        return spline_sin(np.pi - x)
    elif np.pi <= x < 3 * np.pi / 2:
        return -spline_sin(x - np.pi)
    else:
        return -spline_sin(2 * np.pi - x)

# Aproximación para cos(x) utilizando la relación cos(x) = sin(x + pi/2)
def cos_approx(x):
    return sin_approx(x + np.pi / 2)

# Inicializar las secuencias
x = np.zeros(n_iter)
y = np.zeros(n_iter)

# Condiciones iniciales
x[0] = x0
y[0] = y0

# Iteraciones utilizando las aproximaciones polinomiales
for n in range(1, n_iter):
    sin_val = sin_approx(2 * np.pi * y[n-1])
    sin_1_2pi_val = sin_approx(1 - 2 * np.pi * y[n-1])
    cos_val = cos_approx(2 * np.pi * x[n-1])
    
    x[n] = (alpha * sin_val) / (x[n-1] * sin_1_2pi_val)
    x[n] = x[n] % 1  # Aplicar la función mod 1
    
    y[n] = (beta * cos_val) / (y[n-1] * sin_1_2pi_val)
    y[n] = y[n] % 1  # Aplicar la función mod 1

# Graficar el diagrama de fases
plt.figure(figsize=(8, 8))
plt.plot(x, y, 'o', markersize=0.5, color='blue')
plt.title("Diagrama de Fases $x_n$ vs $y_n$ usando Polinomios Spline")
plt.xlabel("$x_n$")
plt.ylabel("$y_n$")
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.grid(True)
plt.gca().set_aspect('equal', adjustable='box')  # Mantener aspecto cuadrado
plt.show()
