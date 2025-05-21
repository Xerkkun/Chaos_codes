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
n_tran = 70000   # Transitorio

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
x[0] = 1.0
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

with open('resultados_prueba2.txt', 'w') as file:
    # Escribir los encabezados de las columnas, separados por tabs
    file.write('j\tx\txTau\n')
    for j in range(n_tran, n+1):
        # Escribir j, x[j], y xTau[j] en el archivo, separados por tabs
        file.write(f'{j}\t{x[j][0]}\t{xTau[j][0]}\n')


# Asumiendo que el cálculo de tu sistema ya se ha realizado y que 'x' y 'xTau' están llenos de datos.

# # Normalizar el eje que se usará para el mapeo de colores
# # Asumiendo que quieres normalizar sobre el rango de x después del transitorio
# norm = plt.Normalize(min(x[n_tran:]), max(x[n_tran:]))
# colors = plt.cm.hsv(norm(x[n_tran:]))

# # Configuración de la figura
# fig, ax = plt.subplots(figsize=(8, 8))
# ax.set_facecolor('black')

# # Graficar usando scatter
# # Aquí, estamos graficando 'x' versus 'xTau', que son equivalentes a las proyecciones 'y' y 'z' del atractor de Lorenz
# sc = ax.scatter(x[n_tran:], xTau[n_tran:], c=colors, s=0.5, alpha=0.6)

# # Ocultar los ejes
# ax.set_xticks([])
# ax.set_yticks([])
# ax.axis('off')  # Apaga los ejes completamente

# # Guardar la figura con un fondo negro
# plt.tight_layout()
# plt.savefig("dynamic_system_projection_98.pdf", format='pdf',
#             bbox_inches='tight', facecolor='black', edgecolor='none')

# # Cerrar la figura
# plt.close()

plt.figure(figsize=(30, 12))
plt.plot(jj[n_tran-1:]*h, x[n_tran:], "m", lw=1)
plt.xlabel(r'$t$', fontsize=24)
plt.ylabel(r'$x(t)$', fontsize=24)

# Cambiar el tamaño de los números de los ejes
# Cambia '12' por el tamaño deseado para el eje x
plt.tick_params(axis='x', labelsize=24)
# Cambia '12' por el tamaño deseado para el eje y
plt.tick_params(axis='y', labelsize=24)

plt.savefig("tmp_96.pdf", format='pdf', dpi=300, bbox_inches='tight')
plt.show()
plt.clf()
