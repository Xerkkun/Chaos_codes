import numpy as np
import math
import matplotlib.pyplot as plt


def fx(x, a, b, fxTau):
    dx = -a*x + b*(fxTau)
    return dx


n = int(1E5)   # Número de pasos totales
mu = 0.96   # Orden de integración

tau = 3.5   # Tiempo de retardo

a = 1.0     # Parámetros
b = 4.5
gamma = 1.0
delta = 0.1

h = 0.01   # Ancho de paso
T = h*n     # Tiempo total de simulación
n_tran = 20000   # Transitorio

x = np.zeros((n+1, 1))       # Vector de la variable de estado
xTau = np.zeros((n+1, 1))    # Vector de la variable de estado con el retardo
nTau = int(tau/h)           # El retaro expresado en pasos de simulación

# Longitud de memoria
Lm = 1.28

# Número de coeficientes binomiales
m = int(Lm/h)

# Principio de memoria corta
if n < m:   # Si la cantidad de pasos totales es menor a la memoria disponible
    nu, mm = 1, n
# Si la cantidad de pasos totales es mayor a la memoria disponible (FPGA)
else:
    nu, mm = m, m

jj = np.linspace(0, n-1, n)
gamma_coef = np.zeros((mm+1, 1))
omega_coef = np.zeros((mm+1, 1))

gamma_coef = np.power(jj+1, 1-mu) - np.power(jj, 1-mu)
omega_coef[0] = gamma_coef[0]

new_x = np.zeros(1)
aux_x = np.zeros((nu, 1))
aux_x[0] = x[0]

# Condición inicial
x[0] = 0.0
jx = n  # Revisar que pasa cuando jx = i

for i in range(1, n+1):  # Recorre todas las iteraciones de 1 a n
    if (jx <= nTau):
        xTau[jx] = 0
        fxTau = 0.55
    else:
        xTau[jx] = x[jx - nTau]
        fxTau = gamma - delta*np.power(xTau[jx], 2)
# ===============================================================
    # Se calculan las sumas de los coeficientes binomiales en cada iteracion
    sum = 0
    if nu == 1:
        omega_coef[i] = gamma_coef[i-1]

        for j in range(1, i):
            omega_coef[j] = gamma_coef[j] - gamma_coef[j-1]
            sum += omega_coef[j]*x[i-j]

        x[i] = (omega_coef[i]*x[0] - sum + np.power(h, mu) *
                math.gamma(2-mu)*fx(x[i-1], a, b, fxTau))/omega_coef[0]
    else:
        for j in range(1, nu+1):  # Solo la cantidad de memoria disponible
            omega_coef[j] = gamma_coef[j] - gamma_coef[j-1]
            sum += omega_coef[j]*aux_x[nu-j, :]

        # Se llenan los valores de x hasta la iteracion m (maximo de memoria)
        if i < mm:
            omega_coef[i] = gamma_coef[i-1]
            x[i] = (omega_coef[i]*x[0] - sum + np.power(h, mu) *
                    math.gamma(2-mu)*fx(x[i-1], a, b, fxTau))/omega_coef[0]
            aux_x[i] = x[i]

        # Luego se van desechando los primeros valores de x que se habian calculado con anterioridad para ir guardando los nuevos
        # Si el  numero de iteracion (i) ya excede el valor de la memoria (m)
        # Desplazar todos los valores del vecto x una posicion hacia arriba y dejar vacio el último elemento del vector
        else:
            new_x = (omega_coef[nu-1]*x[0] - sum + np.power(h, mu)
                     * math.gamma(2-mu)*fx(x[i-1], a, b, fxTau))/omega_coef[0]
            aux_x = np.concatenate((aux_x[1:], [new_x]), axis=0)
            x[i] = new_x

# ===============================================================
    jx = i

with open('resultados_prueba3.txt', 'w') as file:
    # Escribir los encabezados de las columnas, separados por tabs
    file.write('j\tx\txTau\n')
    for j in range(n_tran, n+1):
        # Escribir j, x[j], y xTau[j] en el archivo, separados por tabs
        file.write(f'{j}\t{x[j][0]}\t{xTau[j][0]}\n')
