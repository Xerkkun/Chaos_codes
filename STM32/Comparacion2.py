import os
import numpy as np
import matplotlib.pyplot as plt

# Especifica la ruta de la carpeta donde se encuentran los archivos
carpeta = r'D:\INAOE\Doctorado\STM32'  # Cambia esta ruta según corresponda

# Construir las rutas completas para cada archivo
archivo1 = os.path.join(carpeta, 'GL_short_C.rnd')
archivo2 = os.path.join(carpeta, 'EFORK_c.txt')

# Cargar los datos de cada archivo
# Se asume que cada archivo tiene 4 columnas: tiempo, x, y, z
datos1 = np.loadtxt(archivo1)
datos2 = np.loadtxt(archivo2)

# Extraer las primeras 100 filas: columna 0 para el tiempo y columna 1 para x
t1, x1 = datos1[:10000, 0], datos1[:10000, 1]
t2, x2 = datos2[:10000, 0], datos2[:10000, 1]

# Graficar las soluciones
plt.figure(figsize=(10, 6))
plt.plot(t1, x1, label='GL', linestyle='-')
plt.plot(t2, x2, label='EFORK', linestyle='-')
plt.xlabel('Tiempo')
plt.ylabel('x')
plt.title('Comparación de los métodos')
plt.legend()
plt.grid(True)
plt.show()
