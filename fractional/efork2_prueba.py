import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma
import mittag_leffler as ml

# Definición del sistema de Lorenz

def function(t, Y, alpha):
    y = Y
    dydt = -y + (t**(4-alpha))/gamma(5-alpha)
    return np.array([dydt])

# Parámetros del método EFORK
alpha = 1/3  # Orden fraccionario
h = 1/320    # Tamaño del paso
t0 = 0.0    # Tiempo inicial
T = 50.0    # Tiempo final
y0 = np.array([0.0])  # Condiciones iniciales

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

        K1 = h**alpha * f(tn, yn, alpha)
        K2 = h**alpha * f(tn + c2*h, yn + a21*K1, alpha)

        y[n+1] = yn + w1*K1 + w2*K2

    return t, y


# Ejecutar el método EFORK de 2 pasos en el sistema de Lorenz
t, Y = efork_2_steps_system(function, t0, T, y0, h, alpha)

time = np.linspace(t0, T, int((T-t0)/h))
f = (time**4) * ml.ml(-time**alpha, alpha, 5)
# Graficar el atractor de Lorenz
plt.plot(Y)
plt.plot(f)
plt.show()
