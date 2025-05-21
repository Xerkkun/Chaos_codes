import numpy as np
import matplotlib.pyplot as plt

# Parámetros
alpha = 10.75
beta = 2.45
n_iter = 20000  # Número de iteraciones

# Condiciones iniciales
x0 = 0.412
y0 = 0.784

# Definición de los polinomios cúbicos para la aproximación de sin(x)
def spline_sin(x):
    pi_6 = np.pi / 6
    pi_4 = np.pi / 4
    pi_3 = np.pi / 3
    pi_2 = np.pi / 2

    if 0 <= x < pi_6:
        return 0 + 1.0049 * x + -0.0201 * x**2 + -0.1440 * x**3
    elif pi_6 <= x < pi_4:
        x_adj = x - pi_6
        return 0.5 + 0.8654 * x_adj + -0.2463 * x_adj**2 + -0.1440 * x_adj**3
    elif pi_4 <= x < pi_3:
        x_adj = x - pi_4
        return 0.7071 + 0.7069 * x_adj + -0.3594 * x_adj**2 + -0.0837 * x_adj**3
    elif pi_3 <= x <= pi_2:
        x_adj = x - pi_3
        return 0.8660 + 0.5014 * x_adj + -0.4252 * x_adj**2 + -0.0837 * x_adj**3
    else:
        raise ValueError("Input x is out of the expected range for this spline approximation.")

# Extensión de la aproximación para el intervalo [0, 2*pi]
def sin_approx(x):
    x = x % (2 * np.pi)
    if 0 <= x < np.pi / 2:
        return spline_sin(x)
    elif np.pi / 2 <= x < np.pi:
        return spline_sin(np.pi - x)
    elif np.pi <= x < 3 * np.pi / 2:
        return -spline_sin(x - np.pi)
    else:
        return -spline_sin(2 * np.pi - x)

# Aproximación para cos(x) utilizando la relación cos(x) = sin(x + pi/2)
def cos_approx(x):
    return sin_approx(x + np.pi / 2)

# Inicializar las secuencias
x = np.zeros(n_iter)
y = np.zeros(n_iter)

# Condiciones iniciales
x[0] = x0
y[0] = y0

# Iteraciones utilizando las aproximaciones polinomiales
for n in range(1, n_iter):
    sin_val = sin_approx(2 * np.pi * y[n-1])
    sin_1_2pi_val = sin_approx(1 - 2 * np.pi * y[n-1])
    cos_val = cos_approx(2 * np.pi * x[n-1])

    x[n] = (alpha * sin_val) / (x[n-1] * sin_1_2pi_val)
    x[n] = x[n] % 1  # Aplicar la función mod 1

    y[n] = (beta * cos_val) / (y[n-1] * sin_1_2pi_val)
    y[n] = y[n] % 1  # Aplicar la función mod 1

# Función para leer archivos hexadecimales y convertirlos a valores decimales
def read_hex_file_to_decimal(file_path, scale_factor):
    with open(file_path, 'r') as file:
        hex_values = file.read().splitlines()
    decimal_values = []
    for hex_str in hex_values:
        # Ignorar líneas vacías
        if len(hex_str.strip()) == 0:
            continue
        # Convertir cadena hexadecimal a entero
        val = int(hex_str, 16)
        # Ajuste para números negativos en complemento a dos (32 bits)
        if val >= 2**31:
            val -= 2**32
        # Convertir a decimal con escala
        decimal_val = val / scale_factor
        decimal_values.append(decimal_val)
    return decimal_values

# Rutas de los archivos
x_file_path = 'C:/My_Designs/chaotic_map/chaotic_map/Chaotic_map_p_Xn1.txt'
y_file_path = 'C:/My_Designs/chaotic_map/chaotic_map/Chaotic_map_p_Yn1.txt'

# Leer y convertir datos de los archivos
# Asumiendo que el formato de punto fijo es Q4.27, el factor de escala es 2^27
scale_factor = 2**27

x_file = read_hex_file_to_decimal(x_file_path, scale_factor)
y_file = read_hex_file_to_decimal(y_file_path, scale_factor)

# Convertir listas a arrays de NumPy
x_file = np.array(x_file)
y_file = np.array(y_file)

# Graficar los diagramas de fases lado a lado
plt.figure(figsize=(14, 7))

# Diagrama de fases de Python a la izquierda
plt.subplot(1, 2, 1)
plt.plot(x, y, 'o', markersize=0.5, color='blue')
plt.title("Diagrama de Fases $x_n$ vs $y_n$ (Python)")
plt.xlabel("$x_n$")
plt.ylabel("$y_n$")
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.grid(True)
plt.gca().set_aspect('equal', adjustable='box')  # Mantener aspecto cuadrado

# Diagrama de fases desde el archivo a la derecha
plt.subplot(1, 2, 2)
plt.plot(x_file, y_file, 'o', markersize=0.5, color='green')
plt.title("Diagrama de Fases $x_n$ vs $y_n$ (Archivo)")
plt.xlabel("$x_n$")
plt.ylabel("$y_n$")
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.grid(True)
plt.gca().set_aspect('equal', adjustable='box')  # Mantener aspecto cuadrado

plt.tight_layout()
plt.show()

from sklearn.metrics import mean_squared_error

# Asegurar que las longitudes coincidan
min_length = min(len(x), len(x_file), len(y), len(y_file))
x = x[:min_length]
y = y[:min_length]
x_file = x_file[:min_length]
y_file = y_file[:min_length]

# Calcular el Error Cuadrático Medio
mse_x = mean_squared_error(x, x_file)
mse_y = mean_squared_error(y, y_file)

print(f"Error Cuadrático Medio para x_n: {mse_x}")
print(f"Error Cuadrático Medio para y_n: {mse_y}")

# Calcular diferencias
x_diff = x - x_file
y_diff = y - y_file

# Graficar histogramas de las diferencias
plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
plt.hist(x_diff, bins=50, color='blue', alpha=0.7)
plt.title('Histograma de Diferencias en $x_n$')
plt.xlabel('Diferencia')
plt.ylabel('Frecuencia')

plt.subplot(1, 2, 2)
plt.hist(y_diff, bins=50, color='green', alpha=0.7)
plt.title('Histograma de Diferencias en $y_n$')
plt.xlabel('Diferencia')
plt.ylabel('Frecuencia')

plt.tight_layout()
plt.show()

# Graficar las secuencias x_n y y_n de Python y del archivo en función de las primeras 100 iteraciones
def plot_xn_yn_comparison_100(x_python, y_python, x_file, y_file):
    # Limitar a las primeras 100 iteraciones
    iterations = np.arange(100)  # Primeras 100 iteraciones
    x_python_100 = x_python[:100]
    y_python_100 = y_python[:100]
    x_file_100 = x_file[:100]
    y_file_100 = y_file[:100]

    plt.figure(figsize=(12, 6))

    # Comparación de x_n vs Iteraciones
    plt.subplot(2, 1, 1)
    plt.plot(iterations, x_python_100, color='blue', label='Python $x_n$')
    plt.plot(iterations, x_file_100, color='green', linestyle='--', label='Archivo $x_n$')
    plt.title('Comparación de $x_n$ vs Iteraciones (primeras 100 iteraciones)')
    plt.xlabel('Iteración')
    plt.ylabel('$x_n$')
    plt.legend()
    plt.grid(True)

    # Comparación de y_n vs Iteraciones
    plt.subplot(2, 1, 2)
    plt.plot(iterations, y_python_100, color='red', label='Python $y_n$')
    plt.plot(iterations, y_file_100, color='purple', linestyle='--', label='Archivo $y_n$')
    plt.title('Comparación de $y_n$ vs Iteraciones (primeras 100 iteraciones)')
    plt.xlabel('Iteración')
    plt.ylabel('$y_n$')
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.show()

# Llama a la función para graficar solo las primeras 100 iteraciones comparando Python y Archivo
plot_xn_yn_comparison_100(x, y, x_file, y_file)

# Imprimir los primeros 10 resultados de la simulación y del archivo
def print_first_10_results(x_python, y_python, x_file, y_file):
    print(f"{'Iteración':<10} {'Python x_n':<15} {'Archivo x_n':<15} {'Python y_n':<15} {'Archivo y_n':<15}")
    print("-" * 70)
    for i in range(10):
        print(f"{i:<10} {x_python[i]:<15.10f} {x_file[i]:<15.10f} {y_python[i]:<15.10f} {y_file[i]:<15.10f}")

# Llama a la función para mostrar los primeros 10 resultados
print_first_10_results(x, y, x_file, y_file)
