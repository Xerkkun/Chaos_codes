import numpy as np
import matplotlib.pyplot as plt
import subprocess
from pathlib import Path
from functions import xu
from numerical_methods_v2 import runge_kutta4
from dinamics import equilibrium_points, eigenvalues, time_step

name = 'chen'
dd = 3
#Xu
# a, b, c, d, e, f, g = 4., 1., 1., 1., 1., 4.0, 1.

# Rossler
# a, b, c = 0.2, 0.2, 5.7
a, b, c = 35, 3, 28

# Parámetros del sistema de Lorenz
# sigma = 10.0
# rho = 28.0
# beta = 8./3.  

# Barati 
# a = 0.6
# b = 0.3
# c= 0.5

# =============================================================
# Ajustar la h a partir de los valores propios
# x_eq, y_eq, z_eq = equilibrium_points(a, b, c) 
# eigenvalues_vector = eigenvalues(x_eq, y_eq, z_eq, a, b, c) 
# t_step, transient, num_steps = time_step(eigenvalues_vector)

# z_eq = (c - np.sqrt((c**2)-(4*a*b)))/(2*a)
# x_eq = a*z_eq
# y_eq = -z_eq

# J = np.array([[0, -1, -1], [1, a, 0], [z_eq, 0, x_eq - c]])

x_eq = np.sqrt(2*b*c-a*b)
y_eq = x_eq
z_eq = np.power(x_eq, 2)/b

J = np.array([[-a, a, 0], [c-a-z_eq, c, -x_eq], [y_eq, x_eq, -b]])

eigenvalues, _ = np.linalg.eig(J)
# Verificar que si se genere caos

vp = np.array( [eigenvalues.real,eigenvalues.imag], dtype=float )
vp_non_zero = vp != 0 
vp_inverse = 1/np.abs(vp[vp_non_zero])
vp_min = np.amin(vp_inverse)
vp_max = np.amax(vp_inverse)
t_step = np.round(vp_min/10, 5)

transient = int(vp_max*10/t_step)
steady_state = int(5000*vp_max) #int(1E5)
num_steps = transient+steady_state
n = int(1E4)

# Condiciones iniciales (cercanas a los puntos de equilibrio)
x0, y0, z0 = 4.246, 4.728, 13.470

# Inicializar vector de soluciones
sol = np.zeros((num_steps+1, 3))
sol[0] = [x0, y0, z0] 

# Tiempo de integración
t_start = 0
t_end = num_steps*t_step
t = np.linspace(t_start, t_end, num_steps+1)
x_t = np.zeros((num_steps+1,1))
x_t[:,0] = t

# Resolver el sistema de Lorenz usando Runge-Kutta 4th
sol, _ = runge_kutta4(sol, x_t, num_steps, t_step, a, b, c, dd, name)

# TISEAN        
# Se guarda la evolucion temporal
output_file = name + "_rk4.rnd"
np.savetxt(output_file, sol, delimiter="  ")

# Rossler, no borrar
k_values = (100, 100)
r_values = (0.01, 0.01)

DKY_sum = 0
for i in range (0,2):
   
    # Se llama a Tisean para calcular el espectro de exponentes de lyapunov y la dimension KY
    lyap_specPATH = 'lyap_spec  ' + name + '_rk4.rnd -x' + \
        str(transient) + ' -c1,2,3 -m3,1 -k' + \
        str(k_values[i]) + ' -o' + name + '.lyaps'
    print(lyap_specPATH)
    #check_output(str(lyap_specPATH), shell=True)

    # Ejecuta el comando de lyap_spec y redirige la salida estándar y la salida de error a /dev/null
    result = subprocess.run(lyap_specPATH, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # Se obtiene la D-KY del archivo de salida de TISEAN
    with open(name + '.lyaps', "r") as file:
        file.seek(0,2)
        file.seek(file.tell()-9)
        DKY = float(file.read())
        file.seek(0,2)
        file.seek(file.tell()-309)
        LE = float(file.read(12))
        
    DKY_sum = DKY_sum + DKY
    
DKY_mean = DKY_sum/len(k_values)

print('Ancho de paso: ', t_step)
print('Número de pasos: ', num_steps)
print('Transitorio:, ', transient)
print('Condiciones iniciales: ', x0, y0, z0)
print('Dimensión de Kaplan-Yorke: ', DKY_mean)