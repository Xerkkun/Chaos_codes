import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma

# Parámetros del sistema de Lorenz
sigma = 10.0
rho = 28.0
beta = 8.0 / 3.0

# Definición del sistema de Lorenz


def lorenz(t, Y):
    x, y, z = Y
    dxdt = sigma * (y - x)
    dydt = x * (rho - z) - y
    dzdt = x * y - beta * z
    return np.array([dxdt, dydt, dzdt])


# Parámetros del método EFORK
alpha = 0.88  # Orden fraccionario
h = 0.005    # Tamaño del paso
t0 = 0.0    # Tiempo inicial
T = 50.0    # Tiempo final
y0 = np.array([1.0, 1.0, 1.0])  # Condiciones iniciales

# Coeficientes del método 
c2 = 0.5
a21 = c2**alpha / gamma(alpha + 1)
w2 = gamma(alpha + 1) / (c2**alpha * gamma(2 * alpha + 1))
w1 = 1 / gamma(alpha + 1) - w2

# Función para el método EFORK de 2 pasos para sistemas de ecuaciones


def efork_2_steps_system(f, t0, T, y0, h, alpha):
    N = int((T - t0) / h)
    t = np.linspace(t0, T, N+1)
    y = np.zeros((N+1, len(y0)))
    y[0] = y0

    for n in range(N):
        tn = t[n]
        yn = y[n]

        K1 = h**alpha * f(tn, yn)
        K2 = h**alpha * f(tn + c2*h, yn + a21*K1)

        y[n+1] = yn + w1*K1 + w2*K2

    return t, y


# Ejecutar el método EFORK de 2 pasos en el sistema de Lorenz
t, Y = efork_2_steps_system(lorenz, t0, T, y0, h, alpha)

# Graficar el atractor de Lorenz
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')
ax.plot(Y[:, 0], Y[:, 1], Y[:, 2])
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Atractor de Lorenz usando el método EFORK de 2 pasos')
plt.show()
