import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Para proyecciones 3D

# --- Graficar la trayectoria del sistema (archivo de salida) ---

# Especifica la ruta del archivo generado por el programa en C
# ...existing code...
output_file = r"C:\Users\moren\Desktop\Chaos_codes\Chaos_codes\STM32\chen_variables_EFORK_C.rnd"
# ...existing code...
# Cargar los datos. Se asume que cada línea tiene: tiempo, x, y, z
data = np.loadtxt(output_file)

# Extraer columnas: tiempo, x, y, z
t = data[:, 0]
x = data[:, 1]
y = data[:, 2]
z = data[:, 3]

# # Crear figura 3D
# fig = plt.figure(figsize=(10, 8))
# ax = fig.add_subplot(111, projection='3d')
# ax.plot(x, y, z, lw=0.5)
# ax.set_xlabel("x")
# ax.set_ylabel("y")
# ax.set_zlabel("z")
# ax.set_title("Atractor de Lorenz (Método Grünwald-Letnikov)")
# plt.show()
# plt.clf()

plt.figure (figsize=(12, 4))
plt.subplot(1, 3, 1)
plt.plot(x, y, "m", lw=0.3)
plt.xlabel("x")
plt.ylabel("y")

plt.subplot(1, 3, 2)
plt.plot(x, z, "m", lw=0.3)
plt.xlabel("x")
plt.ylabel("z")

plt.subplot(1, 3, 3)
plt.plot(y, z, "m", lw=0.3)
plt.xlabel("y")
plt.ylabel("z")

plt.tight_layout()
plt.savefig("chen_atractores_EFORK_C.pdf", dpi=400, bbox_inches='tight')
plt.show()  # Se muestra la gráfica al final
plt.clf()