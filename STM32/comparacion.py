import os
import numpy as np
import matplotlib.pyplot as plt

# Especifica la ruta de la carpeta donde se encuentran los archivos
carpeta = r'D:\INAOE\Doctorado\STM32'  # Cambia esta ruta según corresponda

# Construir las rutas completas para cada archivo
archivo1 = os.path.join(carpeta, 'Matlab.txt')
archivo2 = os.path.join(carpeta, 'EFORK.txt')
archivo3 = os.path.join(carpeta, 'EFORK_short.txt')
archivo4 = os.path.join(carpeta, 'EFORK_short_optimized.txt')
archivo5 = os.path.join(carpeta, 'GL_short.rnd')

# Cargar los datos de cada archivo
# Se asume que cada archivo tiene 4 columnas: tiempo, x, y, z
datos1 = np.loadtxt(archivo1)
datos2 = np.loadtxt(archivo2)
datos3 = np.loadtxt(archivo3)
datos4 = np.loadtxt(archivo4)
datos5 = np.loadtxt(archivo5)

# Extraer las primeras 100 filas: columna 0 para el tiempo y columna 1 para x
t1, x1 = datos1[:5000, 0], datos1[:5000, 1]
t2, x2 = datos2[:5000, 0], datos2[:5000, 1]
t3, x3 = datos3[:5000, 0], datos3[:5000, 1]
t4, x4 = datos4[:5000, 0], datos4[:5000, 1]
t5, x5 = datos5[:5000, 0], datos5[:5000, 1]

# Graficar las soluciones
plt.figure(figsize=(10, 6))
plt.plot(t1, x1, label='Fde12', linestyle='-')
plt.plot(t2, x2, label='EFORK', linestyle='-')
plt.plot(t3, x3, label='EFORK_short', linestyle='-')
plt.plot(t4, x4, label='EFORK_short_optimized', linestyle='-')
plt.plot(t5, x5, label='GL_short', linestyle='-')
plt.xlabel('Tiempo')
plt.ylabel('x')
plt.title('Comparación de los métodos')
plt.legend()
plt.grid(True)
plt.show()
