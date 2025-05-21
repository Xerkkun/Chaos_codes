import numpy as np
from scipy.special import gamma
import matplotlib.pyplot as plt
import mittag_leffler as ml


def function(t, Y, alpha):
    y = Y
    dydt = -y + (t**(4-alpha))/gamma(5-alpha)
    return np.array([dydt])


def lorenz(t, Y):
    x, y, z = Y
    sigma = 10.0
    beta = 8.0/3.0
    rho = 28.0
    dxdt = sigma * (y - x)
    dydt = x * (rho - z) - y
    dzdt = x * y - beta * z
    return np.array([dxdt, dydt, dzdt])


# Parámetros del método EFORK
alpha = 0.99  # Orden fraccionario
h = 0.001   # Tamaño del paso
t0 = 0.0    # Tiempo inicial
T = 50.0    # Tiempo final
#y0 = np.array([0.0])  # Condiciones iniciales
y0 = np.array([1.0, 1.0, 1.0])

# Coeficientes del método (ajustados según las fórmulas)
c2 = 0.5
c3 = 0.75
a21 = c2**alpha / gamma(alpha + 1)
a31 = c3**alpha / gamma(alpha + 1)
a32 = (c3 - c2)**alpha / gamma(alpha + 1)
w2 = gamma(alpha + 1) / (c2**alpha * gamma(2 * alpha + 1))
w3 = gamma(alpha + 1) / (c3**alpha * gamma(3 * alpha + 1))
w1 = 1 / gamma(alpha + 1) - w2 - w3

# Función para el método EFORK de 3 pasos para sistemas de ecuaciones


def efork_3_steps_system(f, t0, T, y0, h, alpha):
    N = int((T - t0) / h)
    t = np.linspace(t0, T, N+1)
    y = np.zeros((N+1, len(y0)))
    y[0] = y0

    for n in range(N):
        tn = t[n]
        yn = y[n]

        K1 = h**alpha * f(tn, yn, alpha)
        K2 = h**alpha * f(tn + c2*h, yn + a21*K1, alpha)
        K3 = h**alpha * f(tn + c3*h, yn + a31*K1 + a32*K2, alpha)

        y[n+1] = yn + w1*K1 + w2*K2 + w3*K3

    return t, y


# Ejecutar el método EFORK de 3 pasos en el sistema de Lorenz
t, Y = efork_3_steps_system(function, t0, T, y0, h, alpha)

time = np.linspace(t0, T, int((T-t0)/h))
f = (time**4) * ml.ml(-time**alpha, alpha, 5)
# Graficar el atractor de Lorenz
# plt.plot(Y)
# plt.plot(f)
# plt.show()
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(Y[:,0], Y[:,1], Y[:,2])
plt.show()
