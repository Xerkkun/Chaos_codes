import matplotlib.pyplot as plt
import numpy as np

def read_hex_file(file_path):
    with open(file_path, 'r') as file:
        # Omitir la primera línea del encabezado
        next(file)
        hex_values = file.read().splitlines()
    return hex_values

def hex_to_fixed_point_32bit(hex_values):
    decimals_x_in = []
    decimals_sin_out = []
    decimals_cos_out = []
    for line in hex_values:
        # Separa los valores de x_in, sin_out y cos_out
        x_in_hex, sin_out_hex, cos_out_hex = line.split()

        # Convertir hexadecimal a entero de 32 bits
        x_in_val = int(x_in_hex, 16)
        sin_out_val = int(sin_out_hex, 16)
        cos_out_val = int(cos_out_hex, 16)

        # Ajuste para números negativos en complemento a dos (32 bits con 1 bit de signo)
        if x_in_val >= 2**31:
            x_in_val -= 2**32
        if sin_out_val >= 2**31:
            sin_out_val -= 2**32
        if cos_out_val >= 2**31:
            cos_out_val -= 2**32

        # Convertir los valores a punto flotante dividiendo por 2^27
        x_in_fixed = x_in_val / (2**27)
        sin_out_fixed = sin_out_val / (2**27)
        cos_out_fixed = cos_out_val / (2**27)

        # Agregar los valores convertidos a las listas
        decimals_x_in.append(x_in_fixed)
        decimals_sin_out.append(sin_out_fixed)
        decimals_cos_out.append(cos_out_fixed)

    return decimals_x_in, decimals_sin_out, decimals_cos_out

def plot_comparison(x_in_values, sin_out_values, cos_out_values):
    # Calcular el seno y coseno originales usando numpy
    sine_original = np.sin(x_in_values)
    cosine_original = np.cos(x_in_values)

    # Filtrar los valores para descartar los que están alejados de la solución original
    tolerance = 0.1  # Definir una tolerancia para el error
    filtered_indices = [i for i in range(len(x_in_values)) if abs(sin_out_values[i] - sine_original[i]) <= tolerance and abs(cos_out_values[i] - cosine_original[i]) <= tolerance]

    x_in_filtered = [x_in_values[i] for i in filtered_indices]
    sin_out_filtered = [sin_out_values[i] for i in filtered_indices]
    cos_out_filtered = [cos_out_values[i] for i in filtered_indices]
    sine_original_filtered = [sine_original[i] for i in filtered_indices]
    cosine_original_filtered = [cosine_original[i] for i in filtered_indices]

    # Crear las gráficas
    plt.figure(figsize=(12, 8))

    # Gráfica del seno
    plt.subplot(2, 1, 1)
    plt.plot(x_in_filtered, sine_original_filtered, 'b-', label='Original Sine')
    plt.plot(x_in_filtered, sin_out_filtered, 'r-', label='Approximated Sine (Q4.27)')
    plt.xlabel('x_in (radians)')
    plt.ylabel('Decimal Value')
    plt.legend()
    plt.grid(True)

    # Gráfica del coseno
    plt.subplot(2, 1, 2)
    plt.plot(x_in_filtered, cosine_original_filtered, 'b-', label='Original Cosine')
    plt.plot(x_in_filtered, cos_out_filtered, 'g-', label='Approximated Cosine (Q4.27)')
    plt.xlabel('x_in (radians)')
    plt.ylabel('Decimal Value')
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.show()

# Ruta del archivo de salida generado por el testbench
output_file = 'C:/My_Designs/chaotic_map/chaotic_map/sine_results.txt'

# Leer los valores hexadecimales del archivo
hex_values = read_hex_file(output_file)

# Convertir los valores hexadecimales a punto fijo (32 bits con 27 bits fraccionarios)
x_in_values, sin_out_values, cos_out_values = hex_to_fixed_point_32bit(hex_values)

# Graficar la señal del seno y coseno aproximado y el original
plot_comparison(x_in_values, sin_out_values, cos_out_values)
