import numpy as np
import math
import matplotlib.pyplot as plt


def fx(x, a, b, fxTau):
    dx = -a * x + b * (fxTau)
    return dx


n = 500000
mu = 0.96
jj = np.linspace(0, n - 1, n)
gamma_coef = np.zeros(n + 1)
omega_coef = np.zeros(n + 1)

tau = 3.5
a = 1.0
b = 4.5

gamma = 1.0
delta = 0.1

h = 0.001
T = h * n
n_tran = 20000

gamma_coef = np.power(jj + 1, 1 - mu) - np.power(jj, 1 - mu)
omega_coef[0] = gamma_coef[0]

x = np.zeros(n + 1)
xTau = np.zeros(n + 1)
nTau = int(tau * n / T)

# Longitud de memoria
Lm = 1.024

# Número de coeficientes binomiales
m = int(Lm / h)

# Principio de memoria corta
if n < m:
    nu, mm = 1, n
else:
    nu, mm = m, m

aux_x = np.zeros(nu)
aux_x[0] = x[0]

# Condición inicial
x[0] = 0.0
jx = n  # Revisar que pasa cuando jx = i
gamma_mu_h = np.power(h, mu) * math.gamma(2 - mu)
for i in range(1, n + 1):
    if jx <= nTau:
        xTau[jx] = 0
        fxTau = 0.55
    else:
        xTau[jx] = x[jx - nTau]
        fxTau = gamma - delta * np.power(xTau[jx], 2)

    omega_coef[1:i] = gamma_coef[1:i] - gamma_coef[:i - 1]
    omega_coef[i] = gamma_coef[i - 1]

    # Se calculan las sumas de los coeficientes binomiales en cada iteracion
    if nu == 1:
        sum_omega_x = np.dot(omega_coef[1:i], x[i - 1:0:-1])
        x[i] = (omega_coef[i] * x[0] - sum_omega_x + gamma_mu_h *
                fx(x[i - 1], a, b, fxTau)) / omega_coef[0]
    else:
        sum_omega_aux = np.dot(omega_coef[1:nu + 1], aux_x[::-1])
        # Se llenan los valores de x hasta la iteracion m (maximo de memoria)
        if i < mm:
            x[i] = (omega_coef[i] * x[0] - sum_omega_aux +
                    gamma_mu_h * fx(x[i - 1], a, b, fxTau)) / omega_coef[0]
            aux_x[i] = x[i]
        else:
            new_x = (omega_coef[i] * x[0] - sum_omega_aux +
                     gamma_mu_h * fx(x[i - 1], a, b, fxTau)) / omega_coef[0]
            aux_x = np.roll(aux_x, -1)
            aux_x[-1] = new_x
            x[i] = new_x

    jx = i

# Guardar resultados
np.savetxt('resultados_prueba2.txt', np.column_stack((np.arange(n_tran, n + 1),
           x[n_tran:], xTau[n_tran:])), fmt='%d\t%f\t%f', header='j\tx\txTau', comments='')
