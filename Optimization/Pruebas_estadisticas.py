
def t_test():

    import numpy as np
    from scipy.stats import ttest_ind

    # Datos de ejemplo
    data1 = np.array([20, 21, 22, 23, 24, 25, 26, 27, 28, 29])
    data2 = np.array([30, 31, 32, 33, 34, 35, 36, 37, 38, 39])

    # Realizar el t-test
    stat, p = ttest_ind(data1, data2)

    print('Estadístico t:', stat)
    print('P-valor:', p)

    # Interpretar el resultado
    alpha = 0.05
    if p > alpha:
        print('No hay evidencia suficiente para rechazar la hipótesis nula (H0), por lo que se asume que las medias son iguales.')
    else:
        print('Se rechaza la hipótesis nula (H0), lo que indica que las medias son diferentes.')
    return stat,p

import numpy as np

# Lista para almacenar los valores de la última columna de cada archivo
all_files_values = []

# Iterar sobre los archivos numerados del 0 al 9
for i in range(10):
    # Nombre del archivo basado en el número de iteración
    filename = f'nombre_del_archivo_{i}.txt'
    
    # Lista para almacenar los valores de la última columna de este archivo
    last_column_values = []
    
    # Abrir el archivo de texto y leerlo línea por línea
    with open(filename, 'r') as file:
        for line in file:
            # Ignorar líneas que no contienen datos
            if line.startswith("=") or "n_gen" in line or "Elapsed time" in line:
                continue
            
            # Dividir la línea en sus componentes individuales
            parts = line.split('|')
            
            # Tomar el último valor, eliminar espacios y convertirlo a un número flotante
            value = float(parts[-1].strip())
            
            # Agregar el valor a la lista
            last_column_values.append(value)
            
            # Si ya hemos almacenado 10 líneas de este archivo, rompemos el bucle
            if len(last_column_values) == 10:
                break
    
    # Agregar los valores de este archivo a la lista principal
    all_files_values.append(last_column_values)

# Convertir la lista de listas a un numpy array 2D
all_files_array = np.array(all_files_values)

print(all_files_array)
