import numpy as np
import sys
import matplotlib.pyplot as plt
import subprocess
from pathlib import Path
from scipy.fftpack import fft, fftfreq
from scipy.signal import find_peaks
import time
# =============================================================
def rossler(X, a, b, c):
    """Ecuaciones del sistema de Rossler"""
    x, y, z = X
    dx = -y - z
    dy = x + a*y
    dz = b + z*(x - c)
    return [dx, dy, dz]

def chen(X, a, b, c):
    x, y, z = X
    dx = a*(y - x)
    dy = (c - a)*x - x*z + c*y
    dz = x*y - b*z
    return [dx, dy, dz]

# =============================================================
def runge_kutta4(X, h, a, b, c):
    """Método de Runge Kutta de 4to orden"""
    k1 = rossler(X, a, b, c)
    k2 = rossler(X + h*0.5*np.array(k1), a, b, c)
    k3 = rossler(X + h*0.5*np.array(k2), a, b, c)
    k4 = rossler(X + h*np.array(k3), a, b, c)
    X = X + (h/6.)*(np.array(k1) + 2*np.array(k2) + 2*np.array(k3) + np.array(k4))
    return X


# =============================================================
# Rossler
a, b, c = 5.125517129056624821e-01,  5.724468047480300026e-01,  2.265189065493023790e+00
# Ajustar la h a partir de los valores propios
# Rossler
if ((c**2)-(4*a*b) >= 0):
    z_eq = (c - np.sqrt((c**2)-(4*a*b)))/(2*a)
    x_eq = a*z_eq
    y_eq = -z_eq
else:
    x_eq = 0
    y_eq = 0
    z_eq = 0

J = np.array([[0, -1, -1], [1, a, 0], [z_eq, 0, x_eq - c]])

# Chen
# x_eq = np.sqrt(2*b*c-a*b)
# y_eq = x_eq
# z_eq = np.power(x_eq,2)/b

# J = np.array([[-a, a, 0], [c-a-z_eq, c, -x_eq], [y_eq, x_eq, -b]])
eigenvalues, _ = np.linalg.eig(J)

vp = np.array([eigenvalues.real, eigenvalues.imag], dtype=float)
vp_non_zero = vp != 0
vp_inverse = 1/np.abs(vp[vp_non_zero])
vp_min = np.amin(vp_inverse)
vp_max = np.amax(vp_inverse)
t_step = np.round(vp_min/10, 5)

transient = int(vp_max*5/t_step)
steady_state = int(17000) #int(500*vp_max) #int(1E5)
num_steps = transient+steady_state
n = int(1E4)

print('Ancho de paso: ', t_step)
print('Número de pasos: ', num_steps)
print('Transitorio:, ', transient)
print('Estacionario:, ', steady_state)

# Condiciones iniciales 
x0, y0, z0 = 1,1,1

sol = np.zeros((num_steps+1, 3))
sol[0] = [x0, y0, z0]

# Resolver el sistema de Lorenz usando Runge-Kutta 4th
for i in range(num_steps):
    sol[i+1] = runge_kutta4(sol[i], t_step, a, b, c)

# WOLF
#lyap_specPATH = './rossler 1 1 1 ' + str(a) + ' ' + str(b) + ' ' + str(c) + ' 0.001 2 5000 0.05'
lyap_specPATH = './rossler ' + str(x0) + ' ' + str(y0) + ' ' + str(z0) + ' ' + str(a) + ' ' + str(b) + ' ' + str(c) + ' ' + str(t_step) + ' ' + str(transient*t_step) + ' 5000 0.05'

# Ejecuta el comando y captura la salida estándar
result = subprocess.run(lyap_specPATH, shell=True,stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)

# Procesa la salida
output = result.stdout.decode('utf-8').strip()
parts = output.split()

# Inicializa lyapunov_array
lyapunov_array = np.array([0.0, 0.0, 0.0])

# Extrae los exponentes de Lyapunov
lyapunov_exponents = [float(parts[-3]), float(parts[-2]), float(parts[-1])]

print(lyapunov_exponents)
lyapunov_array = np.array(lyapunov_exponents)

# Ordena los exponentes de Lyapunov de mayor a menor
lyapunov_sorted = np.sort(lyapunov_array)[::-1]

j = 0
sum_lyapunov = 0
for exponent in lyapunov_sorted:
    if sum_lyapunov + exponent > 0:
        sum_lyapunov += exponent
        j += 1
    else:
        break

# Calcula la dimensión de Kaplan-Yorke
if j < len(lyapunov_sorted):
    DKY = j + sum_lyapunov / abs(lyapunov_sorted[j])
else:
    DKY = j  # En caso de que todos los exponentes sean positivos
    
if DKY == 3.0:
    DKY = 0.0
# print(lyapunov_array, DKY)

# Transformada de Fourier
# Quitar el transitorio
sol2 = sol[transient:, :]
t2 = np.linspace(t_step*transient, t_step*num_steps, num_steps+1-transient)
dt2 = t2[1] - t2[0]

Y2 = fft(sol2[:, 0]) / (num_steps+1-transient)  # Transformada normalizada
Y3 = fft(sol2[:, 1]) / (num_steps+1-transient)  # Transformada normalizada
Y4 = fft(sol2[:, 2]) / (num_steps+1-transient)  # Transformada normalizada
frq2 = fftfreq(num_steps+1-transient, dt2)
sumax = sum(abs(Y2))
sumay = sum(abs(Y3))
sumaz = sum(abs(Y4))

# Encuentra los picos
peaksx, _ = find_peaks(np.abs(Y2[0:len(sol2[:, 0])//2]))
peaksy, _ = find_peaks(np.abs(Y3[0:len(sol2[:, 1])//2]))
peaksz, _ = find_peaks(np.abs(Y4[0:len(sol2[:, 2])//2]))
num_peaksx = len(peaksx)
num_peaksy = len(peaksy)
num_peaksz = len(peaksz)

# Encuentra picos a partir de un umbral
# Calcular la magnitud de la transformada de Fourier
magnitudex = np.abs(Y2)
magnitudey = np.abs(Y3)
magnitudez = np.abs(Y4)

# Encontrar los índices de los picos que sobrepasan un valor específico
threshold = 6
peaks_thrx, _ = find_peaks(magnitudex, height=threshold)
peaks_thry, _ = find_peaks(magnitudey, height=threshold)
peaks_thrz, _ = find_peaks(magnitudez, height=threshold)

num_peaks_thrx = len(peaks_thrx)
num_peaks_thry = len(peaks_thry)
num_peaks_thrz = len(peaks_thrz)

#Calcular el valor más alto que alcanzan los picos
# max_peak_value_x = np.max(magnitudex[peaksx])
# max_peak_value_y = np.max(magnitudey[peaksy])
# max_peak_value_z = np.max(magnitudey[peaksz])

# print("El valor más alto que alcanzan los picos para x es", max_peak_value_x)
# print("El valor más alto que alcanzan los picos para y es", max_peak_value_y)
# print("El valor más alto que alcanzan los picos para z es", max_peak_value_z)

#Graficar
label_potx = 'pot=' + str(np.round(sumax,2)) + ',peaks=' + str(num_peaksx) + ',' + str(num_peaks_thrx)
label_poty = 'pot=' + str(np.round(sumay,2)) + ',peaks=' + str(num_peaksy) + ',' + str(num_peaks_thry)
label_potz = 'pot=' + str(np.round(sumaz,2)) + ',peaks=' + str(num_peaksz) + ',' + str(num_peaks_thrz)
label_eig = str(np.round(eigenvalues.real, 2) + np.round(eigenvalues.imag, 2) * 1j)
label_dky = str(np.round(DKY,4))
label_h = str(t_step)
title = str(np.round([a,b,c],2))
#file_name = str(j) + ".png"

fig = plt.figure(figsize=(10, 8))

ax1 = fig.add_subplot(321)
ax1.plot(sol[:,0], sol[:,2])
plt.legend([label_dky])
plt.xlabel('Eje x')
plt.ylabel('Eje z')
plt.title(title)
plt.grid(True)

ax2 = fig.add_subplot(322)
ax2.scatter(eigenvalues.real, eigenvalues.imag, marker="o")
plt.legend([label_eig])
plt.xlabel('Re')
plt.ylabel('Im')
plt.axhline(0, color="black")
plt.axvline(0, color="black")

ax1 = fig.add_subplot(323)
ax1.plot(t2, sol2[:,0])
plt.legend([label_h])
ax1.set_xlabel('Tiempo (s)')
ax1.set_ylabel('$x(t)$')

ax2 = fig.add_subplot(324)
ax2.vlines(frq2, 0, np.abs(Y2.imag))
plt.legend([label_potx])
plt.xlim(-2, 2)
plt.xlabel('Frecuencia (Hz)')
plt.ylabel('Im($Y_x$)')

ax2 = fig.add_subplot(325)
ax2.vlines(frq2, 0, np.abs(Y3.imag))
plt.legend([label_poty])
plt.xlim(-2, 2)
plt.xlabel('Frecuencia (Hz)')
plt.ylabel('Im($Y_y$)')

ax2 = fig.add_subplot(326)
ax2.vlines(frq2, 0, np.abs(Y4.imag))
plt.legend([label_potz])
plt.xlim(-2, 2)
plt.xlabel('Frecuencia (Hz)')
plt.ylabel('Im($Y_z$)')

fig.tight_layout()
plt.show()
