import subprocess
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
import glob

# Rutas a los ejecutables (ajusta si están en otro lugar)
ejecutables = {
    "1": "EFORK_final.exe",
    "2": "GL_final.exe"
}

# Nombres de los métodos para archivos y etiquetas
metodo_nombre = {
    "1": "EFORK",
    "2": "GL"
}

# Sistemas disponibles
sistemas = ["lorenz", "rossler", "chen"]

print("Seleccione el metodo de simulacion:")
print("1) EFORK")
print("2) Grünwald-Letnikov")
opcion = input("Opcion (1/2): ").strip()

if opcion not in ejecutables:
    print("Opcion invalida.")
    exit(1)

exe_name = ejecutables[opcion]
metodo_str = metodo_nombre[opcion]

# Ejecutar el método seleccionado
try:
    print(f"Ejecutando: {exe_name}")
    subprocess.run([exe_name], check=True)
except FileNotFoundError:
    print(f"No se encontro el archivo ejecutable: {exe_name}")
    exit(1)
except subprocess.CalledProcessError:
    print("Error al ejecutar el programa.")
    exit(1)

# Buscar el archivo de salida según el método y el sistema
archivos = []
pattern = f"{metodo_str}_{{}}_c.dat"
for sistema in sistemas:
    archivos.extend(glob.glob(pattern.format(sistema)))

if not archivos:
    print(f"No se encontró ningún archivo .dat generado por {metodo_str}.")
    exit(1)

# Tomar el archivo más reciente si hay varios
archivo_salida = max(archivos, key=os.path.getctime)
print(f"Cargando archivo: {archivo_salida}")

# Cargar datos
try:
    data = np.loadtxt(archivo_salida)
except Exception as e:
    print(f"Error al cargar datos: {e}")
    exit(1)

# Separar columnas
t, x, y, z = data[:, 0], data[:, 1], data[:, 2], data[:, 3]

# Obtener nombre del sistema para los archivos PDF
nombre_sistema = os.path.basename(archivo_salida).split('_')[1]

# Generar nombres para los PDFs de salida
pdf2d = f"{nombre_sistema}_{metodo_str}_c.pdf"
pdf3d = f"{nombre_sistema}_{metodo_str}_3D_c.pdf"

# Crear gráficas 2D
plt.figure(figsize=(12, 4))
plt.subplot(1, 3, 1)
plt.plot(x, y, lw=0.3)
plt.xlabel('x')
plt.ylabel('y')
plt.subplot(1, 3, 2)
plt.plot(x, z, lw=0.3)
plt.xlabel('x')
plt.ylabel('z')
plt.subplot(1, 3, 3)
plt.plot(y, z, lw=0.3)
plt.xlabel('y')
plt.ylabel('z')
plt.tight_layout()
plt.savefig(pdf2d, dpi=400, bbox_inches='tight')
plt.show()

# Crear gráfica 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x, y, z, lw=0.2)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.savefig(pdf3d, dpi=400, bbox_inches='tight')
plt.show()
