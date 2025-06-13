import subprocess
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
import glob

# Rutas a los ejecutables (ajusta si están en otro lugar)
executables = {
    "1": "EFORK_final.exe",
    "2": "GL_final.exe"
}

método_nombre = {
    "1": "EFORK",
    "2": "GL"
}

sistemas = ["lorenz", "rossler", "chen"]  # Para buscar archivos posibles

print("Seleccione el método de simulación:")
print("1) EFORK")
print("2) Grünwald-Letnikov")
opcion = input("Opción (1/2): ").strip()

if opcion not in executables:
    print("Opción inválida.")
    exit(1)

exe_name = executables[opcion]
metodo_str = método_nombre[opcion]

try:
    print(f"Ejecutando: {exe_name}")
    subprocess.run([exe_name], check=True)
except FileNotFoundError:
    print(f"No se encontró el archivo ejecutable: {exe_name}")
    exit(1)
except subprocess.CalledProcessError:
    print("Error al ejecutar el programa.")
    exit(1)

# Buscar el archivo de salida según los nombres correctos
archivos = []
for sistema in sistemas:
    archivos += glob.glob(f"EFORK_{sistema}_c.txt")

if not archivos:
    print("No se encontró ningún archivo .txt generado por EFORK.")
    exit(1)

# Si hay varios, tomar el más reciente
archivo_salida = max(archivos, key=os.path.getctime)
print(f"Cargando archivo: {archivo_salida}")

# Cargar datos
data = np.loadtxt(archivo_salida)
t, x, y, z = data[:, 0], data[:, 1], data[:, 2], data[:, 3]

# Obtener el nombre del sistema para guardar PDF
nombre_sistema = archivo_salida.split('_')[1]  # "lorenz", "rossler" o "chen"

# Crear nombres de los pdfs según el método
pdf2d = f"{nombre_sistema}_atractores_{metodo_str}_c.pdf"
pdf3d = f"{nombre_sistema}_atractores_{metodo_str}_3D_c.pdf"

# Gráficas 2D
plt.figure(figsize=(12, 4))
plt.subplot(1, 3, 1)
plt.plot(x, y, 'm', lw=0.3)
plt.xlabel('x'); plt.ylabel('y')
plt.subplot(1, 3, 2)
plt.plot(x, z, 'm', lw=0.3)
plt.xlabel('x'); plt.ylabel('z')
plt.subplot(1, 3, 3)
plt.plot(y, z, 'm', lw=0.3)
plt.xlabel('y'); plt.ylabel('z')
plt.tight_layout()
plt.savefig(pdf2d, dpi=400, bbox_inches='tight')
plt.show()

# Gráfica 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x, y, z, 'm', lw=0.2)
ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
plt.savefig(pdf3d, dpi=400, bbox_inches='tight')
plt.show()
