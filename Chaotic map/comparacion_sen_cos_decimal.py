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

    for line in hex_values:
        if line.startswith("x_in") or len(line.strip()) == 0:
            # Ignorar encabezados y líneas vacías
            continue

        # Dividir la línea en valores
        values = line.split()

        # Verificar si hay al menos 4 valores en la línea
        if len(values) < 4:
            print(f"Línea incompleta o inválida: {line}")
            continue  # Ignorar líneas con menos de 4 valores

        # Separa los valores de x_in, sin_out, cos_out, frac_div_out
        x_in_hex, sin_out_hex, cos_out_hex, frac_div_out_hex = values[:4]

        # Convertir hexadecimal a entero de 32 bits
        x_in_val = int(x_in_hex, 16)
        sin_out_val = int(sin_out_hex, 16)
        cos_out_val = int(cos_out_hex, 16)
        frac_div_out_val = int(frac_div_out_hex, 16)

        # Ajuste para números negativos en complemento a dos (32 bits con 1 bit de signo)
        if x_in_val >= 2**31:
            x_in_val -= 2**32
        if sin_out_val >= 2**31:
            sin_out_val -= 2**32
        if cos_out_val >= 2**31:
            cos_out_val -= 2**32
        if frac_div_out_val >= 2**31:
            frac_div_out_val -= 2**32

        # Convertir los valores a punto flotante dividiendo por 2^27
        x_in_fixed = x_in_val / (2**27)
        sin_out_fixed = sin_out_val / (2**27)
        cos_out_fixed = cos_out_val / (2**27)
        frac_div_out_fixed = frac_div_out_val / (2**27)

        # Agregar los valores convertidos a las listas
        decimals_x_in.append(x_in_fixed)
        decimals_sin_out.append(sin_out_fixed)
        decimals_cos_out.append(cos_out_fixed)
        decimals_frac_div_out.append(frac_div_out_fixed)

    return decimals_x_in, decimals_sin_out, decimals_cos_out, decimals_frac_div_out

def plot_division_comparison(x_in_values, frac_div_out_values, python_div_values):
    # Graficar la comparación de la división realizada por Python y la división obtenida del archivo
    plt.figure(figsize=(10, 6))
    plt.plot(x_in_values, python_div_values, 'b-', label='División Realizada por Python')
    plt.plot(x_in_values, frac_div_out_values, 'g-', label='División del Archivo (Q4.27)')
    plt.title('Comparación entre la división del archivo y la división realizada por Python')
    plt.xlabel('x_in (rad)')
    plt.ylabel('Valor Decimal')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

# Ruta del archivo de salida generado por el testbench
output_file = 'C:/My_Designs/chaotic_map/chaotic_map/sine_cosine_div_results_2.txt'

# Leer los valores hexadecimales del archivo
hex_values = read_hex_file(output_file)

# Convertir los valores hexadecimales a punto fijo (32 bits con 27 bits fraccionarios)
x_in_values, sin_out_values, cos_out_values, frac_div_out_values = hex_to_fixed_point_32bit(hex_values)

# Realizar la división con los valores obtenidos por Python
python_div_values = [sin_out / cos_out if cos_out != 0 else 0 for sin_out, cos_out in zip(sin_out_values, cos_out_values)]

# Graficar solo la comparación de la división entre el archivo y Python
plot_division_comparison(x_in_values, frac_div_out_values, python_div_values)
