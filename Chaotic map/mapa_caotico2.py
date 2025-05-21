import numpy as np
import matplotlib.pyplot as plt

# Parámetros y rango de exploración
alpha_values = np.linspace(0.1, 10, 500)
beta_values = np.linspace(0.1, 10, 500)
n_iter = 1000  # Número de iteraciones
last_points = 100  # Puntos finales a graficar

# Condiciones iniciales
x0 = 0.412
y0 = 0.784

# Listas para almacenar los resultados
x_alpha = []
y_alpha = []
x_beta = []
y_beta = []

# Simulación para distintos valores de alpha
for alpha in alpha_values:
    x = x0
    y = y0
    for n in range(n_iter):
        x = (alpha * np.sin(2 * np.pi * y)) / (x * np.sin(1 - 2 * np.pi * y))
        x = x % 1  # Aplicar la función mod 1
        y = (beta_values[-1] * np.cos(2 * np.pi * x)) / (y * np.sin(1 - 2 * np.pi * x))
        y = y % 1  # Aplicar la función mod 1
        if n >= n_iter - last_points:
            x_alpha.append((alpha, x))
            y_alpha.append((alpha, y))

# Simulación para distintos valores de beta
for beta in beta_values:
    x = x0
    y = y0
    for n in range(n_iter):
        x = (alpha_values[-1] * np.sin(2 * np.pi * y)) / (x * np.sin(1 - 2 * np.pi * y))
        x = x % 1  # Aplicar la función mod 1
        y = (beta * np.cos(2 * np.pi * x)) / (y * np.sin(1 - 2 * np.pi * x))
        y = y % 1  # Aplicar la función mod 1
        if n >= n_iter - last_points:
            x_beta.append((beta, x))
            y_beta.append((beta, y))

# Convertir las listas en arrays de numpy para facilitar el plot
x_alpha = np.array(x_alpha)
y_alpha = np.array(y_alpha)
x_beta = np.array(x_beta)
y_beta = np.array(y_beta)

# Graficar los resultados
fig, axs = plt.subplots(2, 2, figsize=(10, 10))

# Gráfico de x_n vs alpha
axs[0, 0].plot(x_alpha[:, 0], x_alpha[:, 1], 'o', markersize=0.5, color='blue')
axs[0, 0].set_title("$x_n$ vs $\\alpha$")
axs[0, 0].set_xlabel("$\\alpha$")
axs[0, 0].set_ylabel("$x_n$")
axs[0, 0].set_xlim([0, 10])
axs[0, 0].set_ylim([0, 1])
axs[0, 0].grid(True)

# Gráfico de y_n vs alpha
axs[0, 1].plot(y_alpha[:, 0], y_alpha[:, 1], 'o', markersize=0.5, color='green')
axs[0, 1].set_title("$y_n$ vs $\\alpha$")
axs[0, 1].set_xlabel("$\\alpha$")
axs[0, 1].set_ylabel("$y_n$")
axs[0, 1].set_xlim([0, 10])
axs[0, 1].set_ylim([0, 1])
axs[0, 1].grid(True)

# Gráfico de x_n vs beta
axs[1, 0].plot(x_beta[:, 0], x_beta[:, 1], 'o', markersize=0.5, color='red')
axs[1, 0].set_title("$x_n$ vs $\\beta$")
axs[1, 0].set_xlabel("$\\beta$")
axs[1, 0].set_ylabel("$x_n$")
axs[1, 0].set_xlim([0, 10])
axs[1, 0].set_ylim([0, 1])
axs[1, 0].grid(True)

# Gráfico de y_n vs beta
axs[1, 1].plot(y_beta[:, 0], y_beta[:, 1], 'o', markersize=0.5, color='purple')
axs[1, 1].set_title("$y_n$ vs $\\beta$")
axs[1, 1].set_xlabel("$\\beta$")
axs[1, 1].set_ylabel("$y_n$")
axs[1, 1].set_xlim([0, 10])
axs[1, 1].set_ylim([0, 1])
axs[1, 1].grid(True)

plt.tight_layout()
plt.show()
