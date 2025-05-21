import matplotlib.pyplot as plt
import numpy as np

def read_hex_file(file_path):
    with open(file_path, 'r') as file:
        hex_values = file.read().splitlines()
    return hex_values

def hex_to_fixed_point_32bit(hex_values):
    decimals_x_in = []
    decimals_sin_out = []
    decimals_cos_out = []
    decimals_frac_div_out = []
    decimals_mod1_out = []  # Añadimos una lista para la quinta columna

    for line in hex_values:
        if line.startswith("x_in") or len(line.strip()) == 0:
            # Ignorar encabezados y líneas vacías
            continue

        # Dividir la línea en valores
        values = line.split()

        # Verificar si hay al menos 5 valores en la línea
        if len(values) < 5:
            print(f"Línea incompleta o inválida: {line}")
            continue  # Ignorar líneas con menos de 5 valores

        # Separa los valores de x_in, sin_out, cos_out, frac_div_out, mod1_out
        x_in_hex, sin_out_hex, cos_out_hex, frac_div_out_hex, mod1_out_hex = values[:5]

        # Convertir hexadecimal a entero de 32 bits
        x_in_val = int(x_in_hex, 16)
        sin_out_val = int(sin_out_hex, 16)
        cos_out_val = int(cos_out_hex, 16)
        frac_div_out_val = int(frac_div_out_hex, 16)
        mod1_out_val = int(mod1_out_hex, 16)

        # Ajuste para números negativos en complemento a dos (32 bits con 1 bit de signo)
        if x_in_val >= 2**31:
            x_in_val -= 2**32
        if sin_out_val >= 2**31:
            sin_out_val -= 2**32
        if cos_out_val >= 2**31:
            cos_out_val -= 2**32
        if frac_div_out_val >= 2**31:
            frac_div_out_val -= 2**32
        if mod1_out_val >= 2**31:
            mod1_out_val -= 2**32

        # Convertir los valores a punto flotante dividiendo por 2^27
        x_in_fixed = x_in_val / (2**27)
        sin_out_fixed = sin_out_val / (2**27)
        cos_out_fixed = cos_out_val / (2**27)
        frac_div_out_fixed = frac_div_out_val / (2**27)
        mod1_out_fixed = mod1_out_val / (2**27)

        # Agregar los valores convertidos a las listas
        decimals_x_in.append(x_in_fixed)
        decimals_sin_out.append(sin_out_fixed)
        decimals_cos_out.append(cos_out_fixed)
        decimals_frac_div_out.append(frac_div_out_fixed)
        decimals_mod1_out.append(mod1_out_fixed)  # Agregar la columna mod1

    return decimals_x_in, decimals_sin_out, decimals_cos_out, decimals_frac_div_out, decimals_mod1_out

def plot_division_comparison(x_in_values, mod1_out_values, python_frac_div_values):
    # Graficar la comparación de la parte fraccionaria de la división calculada por Python y la salida del archivo (mod1)
    plt.figure(figsize=(10, 6))
    plt.plot(x_in_values, python_frac_div_values, 'b-', label='Parte Fraccionaria (Python)')
    plt.plot(x_in_values, mod1_out_values, 'g-', label='Parte Fraccionaria del Archivo (mod1)')
    plt.title('Comparación de la Parte Fraccionaria de la División (Archivo mod1 vs Python)')
    plt.xlabel('x_in (rad)')
    plt.ylabel('Parte Fraccionaria')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

# Ruta del archivo de salida generado por el testbench
output_file = 'C:/My_Designs/chaotic_map/chaotic_map/sine_cosine_div_results_2.txt'

# Leer los valores hexadecimales del archivo
hex_values = read_hex_file(output_file)

# Convertir los valores hexadecimales a punto fijo (32 bits con 27 bits fraccionarios)
x_in_values, sin_out_values, cos_out_values, frac_div_out_values, mod1_out_values = hex_to_fixed_point_32bit(hex_values)

# Filtrar los valores de x_in entre 0 y π
filtered_x_in_values = []
filtered_mod1_out_values = []
filtered_python_frac_div_values = []

# Realizar la división seno/coseno y obtener solo la parte fraccionaria
for x_in, sin_out, cos_out, mod1_out in zip(x_in_values, sin_out_values, cos_out_values, mod1_out_values):
    if 0 <= x_in <= np.pi:  # Filtrar valores de x_in entre 0 y pi
        filtered_x_in_values.append(x_in)
        filtered_mod1_out_values.append(mod1_out)
        if cos_out != 0:
            div_result = sin_out / cos_out
            frac_part = div_result - np.floor(div_result)  # Parte fraccionaria de la división
        else:
            frac_part = 0  # Evitar división por cero
        filtered_python_frac_div_values.append(frac_part)

# Graficar solo la comparación de la parte fraccionaria de la división (usando mod1_out_values) entre 0 y pi
plot_division_comparison(filtered_x_in_values, filtered_mod1_out_values, filtered_python_frac_div_values)
