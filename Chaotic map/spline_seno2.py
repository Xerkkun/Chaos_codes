import numpy as np
import matplotlib.pyplot as plt

# Definición de los subintervalos y coeficientes para la interpolación spline cúbica
def spline_approx(x):
    # Coeficientes obtenidos para cada subintervalo
    a0, b0, c0, d0 = 0, 1.0049, -0.0201, -0.1440
    a1, b1, c1, d1 = 0.5, 0.8654, -0.2463, -0.1440
    a2, b2, c2, d2 = 0.7071, 0.7069, -0.3594, -0.0837
    a3, b3, c3, d3 = 0.8660, 0.5014, -0.4252, -0.0837

    # Definición de los puntos de los subintervalos
    pi_6 = np.pi / 6
    pi_4 = np.pi / 4
    pi_3 = np.pi / 3
    pi_2 = np.pi / 2

    # Normalización de x al rango [0, pi/2]
    if 0 <= x < pi_6:
        return a0 + b0 * x + c0 * x**2 + d0 * x**3
    elif pi_6 <= x < pi_4:
        x_adj = x - pi_6
        return a1 + b1 * x_adj + c1 * x_adj**2 + d1 * x_adj**3
    elif pi_4 <= x < pi_3:
        x_adj = x - pi_4
        return a2 + b2 * x_adj + c2 * x_adj**2 + d2 * x_adj**3
    elif pi_3 <= x <= pi_2:
        x_adj = x - pi_3
        return a3 + b3 * x_adj + c3 * x_adj**2 + d3 * x_adj**3

# Extensión de la aproximación al intervalo [0, 2*pi]
def sin_approx(x):
    # Ajuste x al rango [0, 2*pi]
    x = x % (2 * np.pi)

    # Aplicación de la simetría del seno
    if 0 <= x < np.pi / 2:
        return spline_approx(x)
    elif np.pi / 2 <= x < np.pi:
        return spline_approx(np.pi - x)
    elif np.pi <= x < 3 * np.pi / 2:
        return -spline_approx(x - np.pi)
    else:
        return -spline_approx(2 * np.pi - x)

# Comparación de la función seno original con la aproximación
x_values = np.linspace(0, 2 * np.pi, 1000)
sin_values = np.sin(x_values)
approx_values = np.array([sin_approx(x) for x in x_values])

# Gráfica de comparación
plt.figure(figsize=(10, 6))
plt.plot(x_values, sin_values, label='Original function $\sin(x)$', color='blue')
plt.plot(x_values, approx_values, '--', label='Spline approximation', color='red')
plt.xlabel('x (radians)')
plt.ylabel('$\sin(x)$')
plt.legend()
plt.grid(True)
plt.show()
