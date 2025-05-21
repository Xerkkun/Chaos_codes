import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Para proyecciones 3D

# --- Graficar la trayectoria del sistema (archivo de salida) ---

# Especifica la ruta del archivo generado por el programa en C
output_file = "lorenz_variables.rnd"  # Ajusta la ruta si es necesario

# Cargar los datos. Se asume que cada línea tiene: tiempo, x, y, z
data = np.loadtxt(output_file)

# Extraer columnas: tiempo, x, y, z
t = data[:, 0]
x = data[:, 1]
y = data[:, 2]
z = data[:, 3]

# Crear figura 3D
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
ax.plot(x, y, z, lw=0.5)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.set_title("Atractor de Lorenz (Método Grünwald-Letnikov)")
plt.show()
plt.clf()