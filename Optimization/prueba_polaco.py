import numpy as np
import sys
import matplotlib.pyplot as plt
import subprocess
from pathlib import Path

#sigma, rho, beta = 10.0, 28.0, 8./3.
a, b, c = 5.125517129056624821e-01,  5.724468047480300026e-01 , 2.265189065493023790e+00

lyap_specPATH = './rossler 1 1 1 ' + \
    str(a) + ' ' + str(b) + ' ' + str(c) + ' 0.001 2 5000 0.05'

# Ejecuta el comando y captura la salida estándar
result = subprocess.run(lyap_specPATH, shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)

# Procesa la salida
output = result.stdout.decode('utf-8').strip()
parts = output.split()

# Inicializa lyapunov_array
lyapunov_array = np.array([0.0, 0.0, 0.0])

# Verifica si parts tiene suficientes elementos
if len(parts) >= 3:
    # Extrae los exponentes de Lyapunov
    lyapunov_exponents = [float(parts[-3]), float(parts[-2]), float(parts[-1])]
    lyapunov_array = np.array(lyapunov_exponents)
else:
    print("Error: Formato de salida inesperado.")

print("Exponentes de Lyapunov:", lyapunov_array)
