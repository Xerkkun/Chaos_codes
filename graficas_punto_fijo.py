import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def read_hex_file(file_path):
    with open(file_path, 'r') as file:
        hex_values = file.read().splitlines()
    return hex_values

def hex_to_decimal(hex_values, total_bits=32):
    decimals = []
    for hex_val in hex_values:
        # Convert hexadecimal to integer
        integer_val = int(hex_val, 16)
        # Check if the number is negative in two's complement
        if integer_val >= 2**(total_bits - 1):
            integer_val -= 2**total_bits
        decimals.append(integer_val)
    return decimals

def save_plot(x, y, xlabel, ylabel, title, filename):
    plt.figure(figsize=(10, 6))
    plt.plot(x, y, label=title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.savefig(filename)
    plt.close()

# Usar raw strings para las rutas de archivo
file_x = r'C:\My_Designs\FER_CAOS\FER\FER_CHAOS\Chen_p_Xn1.txt'
file_y = r'C:\My_Designs\FER_CAOS\FER\FER_CHAOS\Chen_p_Yn1.txt'
file_z = r'C:\My_Designs\FER_CAOS\FER\FER_CHAOS\Chen_p_Zn1.txt'

# Leer archivos
x_hex = read_hex_file(file_x)
y_hex = read_hex_file(file_y)
z_hex = read_hex_file(file_z)

# Convertir a decimal (puedes modificar el total_bits aquí)
total_bits = 32
x_dec = hex_to_decimal(x_hex, total_bits)
y_dec = hex_to_decimal(y_hex, total_bits)
z_dec = hex_to_decimal(z_hex, total_bits)

# Guardar las gráficas por separado
save_plot(x_dec, y_dec, 'X', 'Y', 'Gráfica XY del sistema Chen', 'grafica_xy_Chen.png')
save_plot(x_dec, z_dec, 'X', 'Z', 'Gráfica XZ del sistema Chen', 'grafica_xz_Chen.png')
save_plot(y_dec, z_dec, 'Y', 'Z', 'Gráfica YZ del sistema Chen', 'grafica_yz_Chen.png')

# Crear y guardar las gráficas en un archivo PDF
with PdfPages('Chen_plots.pdf') as pdf:
    # Gráfica XY
    plt.figure(figsize=(10, 6))
    plt.plot(x_dec, y_dec, label='XY Plot')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Gráfica XY del sistema Chen')
    plt.legend()
    plt.grid(True)
    pdf.savefig()
    plt.close()

    # Gráfica XZ
    plt.figure(figsize=(10, 6))
    plt.plot(x_dec, z_dec, label='XZ Plot')
    plt.xlabel('X')
    plt.ylabel('Z')
    plt.title('Gráfica XZ del sistema Chen')
    plt.legend()
    plt.grid(True)
    pdf.savefig()
    plt.close()

    # Gráfica YZ
    plt.figure(figsize=(10, 6))
    plt.plot(y_dec, z_dec, label='YZ Plot')
    plt.xlabel('Y')
    plt.ylabel('Z')
    plt.title('Gráfica YZ del sistema Chen')
    plt.legend()
    plt.grid(True)
    pdf.savefig()
    plt.close()

print("Las gráficas se han guardado en 'Chen_plots.pdf' y por separado en archivos PNG.")
