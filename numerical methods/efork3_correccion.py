import numpy as np
from scipy.special import gamma
import matplotlib.pyplot as plt
import time

# Parámetros del sistema de Lorenz
sigma = 10.0
rho = 28.0
beta = 8.0 / 3.0

def lorenz(t, Y):
    x, y, z = Y
    dxdt = sigma * (y - x)
    dydt = x * (rho - z) - y
    dzdt = x * y - beta * z
    return np.array([dxdt, dydt, dzdt])

# Parámetros del método EFORK 
alpha = 0.99  # Orden fraccionario
h = 0.005     # Tamaño del paso
t0 = 0.0      # Tiempo inicial
T = 50.0      # Tiempo final
y0 = np.array([1.0, 1.0, 1.0])  # Condiciones iniciales

# Coeficientes del método (según las fórmulas de efork3)
c2 = 0.5
c3 = 0.75
a21 = c2**alpha / gamma(alpha + 1)
a31 = c3**alpha / gamma(alpha + 1)
a32 = (c3 - c2)**alpha / gamma(alpha + 1)
w2 = gamma(alpha + 1) / (c2**alpha * gamma(2 * alpha + 1))
w3 = gamma(alpha + 1) / (c3**alpha * gamma(3 * alpha + 1))
w1 = 1 / gamma(alpha + 1) - w2 - w3

def fn_corr(k, t_current, Y_current, t_arr, y_arr, alpha, h):
    """
    Calcula la corrección histórica (memoria) para la evaluación de la función en tiempo t_current,
    vectorizando la suma para optimizar la ejecución.
    
    Parámetros:
      k         : número de pasos previos (índice actual)
      t_current : tiempo actual
      Y_current : solución en el tiempo actual (vector)
      t_arr     : vector de tiempos ya calculados
      y_arr     : matriz de soluciones ya calculadas
      alpha     : orden fraccionario
      h         : tamaño del paso
      
    Retorna:
      La evaluación corregida: f_corr = lorenz(t_current, Y_current) - memoria.
    """
    indices = np.arange(k)
    t_i = t_arr[indices]
    t_ip1 = t_arr[indices + 1]
    factor = (t_current - t_i)**(1.0 - alpha) - (t_current - t_ip1)**(1.0 - alpha)
    diff = (y_arr[indices + 1] - y_arr[indices]) / (h * gamma(2 - alpha))
    v = np.sum(diff * factor[:, np.newaxis], axis=0)
    return lorenz(t_current, Y_current) - v

def efork_3_steps_system_with_history(f, t0, T, y0, h, alpha):
    """
    Esquema EFORK de 3 pasos que incorpora corrección histórica para sistemas.
    
    f  : función que evalúa el sistema (en este caso, 'lorenz')
    t0 : tiempo inicial
    T  : tiempo final
    y0 : condiciones iniciales (vector)
    h  : tamaño del paso
    alpha: orden fraccionario
    """
    N = int((T - t0) / h)
    t = np.linspace(t0, T, N + 1)
    y = np.zeros((N + 1, len(y0)))
    y[0] = y0

    # Primer paso sin corrección histórica (no hay historia aún)
    K1 = h**alpha * f(t[0], y[0])
    K2 = h**alpha * f(t[0] + c2 * h, y[0] + a21 * K1)
    K3 = h**alpha * f(t[0] + c3 * h, y[0] + a31 * K1 + a32 * K2)
    y[1] = y[0] + w1 * K1 + w2 * K2 + w3 * K3

    # Variables para el progreso
    progress_step = max(1, N // 10)
    start_time = time.time()

    # Para n >= 1 se incorpora la corrección histórica
    for n in range(1, N):
        tn = t[n]
        yn = y[n]
        f_corr = fn_corr(n, tn, yn, t, y, alpha, h)
        K1 = h**alpha * f_corr
        K2 = h**alpha * fn_corr(n, tn + c2 * h, yn + a21 * K1, t, y, alpha, h)
        K3 = h**alpha * fn_corr(n, tn + c3 * h, yn + a31 * K1 + a32 * K2, t, y, alpha, h)
        y[n+1] = yn + w1 * K1 + w2 * K2 + w3 * K3

        if n % progress_step == 0:
            elapsed = time.time() - start_time
            print(f"Iteración {n}/{N} - Tiempo transcurrido: {elapsed:.2f} s")

    print("Cálculo completado.")
    return t, y

# Ejecutar el método EFORK de 3 pasos con corrección histórica en el sistema de Lorenz
t, Y = efork_3_steps_system_with_history(lorenz, t0, T, y0, h, alpha)

# Guardar la salida en un archivo .dat (una columna para el tiempo y una para cada variable)
data = np.column_stack((t, Y))
np.savetxt("lorenz_output.dat", data, fmt="%.6e", header="t\tx\ty\tz", comments="")

print("Datos guardados en 'lorenz_output.dat'.")

# Graficar el atractor de Lorenz en 3D
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')
ax.plot(Y[:, 0], Y[:, 1], Y[:, 2])
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Atractor de Lorenz con corrección histórica (EFORK 3 pasos)')
plt.show()
