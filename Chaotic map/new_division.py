import numpy as np
import matplotlib.pyplot as plt

# Definir los valores de x e y
x = np.linspace(-1, 1, 500)
y = np.linspace(-1, 1, 500)

# Definir las funciones
f1 = np.sin( y)
f2 = np.cos( x)

# Evitar división por cero en el cociente
with np.errstate(divide='ignore', invalid='ignore'):
    quotient = np.where(f2 != 0, f1 / f2, 0)

# Graficar la primera función
plt.figure(figsize=(15, 5))
plt.subplot(1, 3, 1)
plt.plot(y, f1, label='sin(2*pi*y)')
plt.title('sin(2*pi*y)')
plt.grid(True)
plt.legend()

# Graficar la segunda función
plt.subplot(1, 3, 2)
plt.plot(x, f2, label='x*sin(1-2*pi*x)', color='orange')
plt.title('x*sin(1-2*pi*x)')
plt.grid(True)
plt.legend()

# Graficar el cociente de las dos funciones
plt.subplot(1, 3, 3)
plt.plot(x, quotient, label='sin(2*pi*y) / (x*sin(1-2*pi*x))', color='green')
plt.title('Cociente de f1 y f2')
plt.grid(True)
plt.legend()

# Mostrar las gráficas
plt.tight_layout()
plt.show()
