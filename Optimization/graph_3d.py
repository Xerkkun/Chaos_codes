import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def graph_colors_3d(x, y, z, t, name, method):
    # Normalizar el eje x para usarlo en el mapeo de colores
    norm = plt.Normalize(x.min(), x.max())
    colors = plt.cm.hsv(norm(x))

    # Configurar la figura y el eje con fondo negro
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_facecolor('black')

    # Graficar el atractor con un gradiente de color dependiente de la posición x
    sc = ax.scatter(x, y, z, c=colors, s=0.5, alpha=0.6)

    # Ocultar los ejes
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.axis('off')  # Desactivar completamente el eje

    # Ajustar el layout y guardar en PDF con fondo negro
    plt.tight_layout()
    plt.savefig(name + "_" + method + "_xyz_atractor.pdf", dpi=300,
                format='pdf', bbox_inches='tight', facecolor='black', edgecolor='none')

    # Cerrar la gráfica
    plt.close()


# Ejemplo de uso
x = np.random.rand(100)
y = np.random.rand(100)
z = np.random.rand(100)
# Asumiendo que 't' es un parámetro adicional que podrías querer usar
t = np.random.rand(100)
graph_colors_3d(x, y, z, t, "nombre_del_archivo", "metodo")
