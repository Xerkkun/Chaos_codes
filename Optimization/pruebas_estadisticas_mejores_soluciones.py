# Listas para almacenar los últimos valores de la última columna de todos los archivos
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon, ttest_ind

def main():
    all_files_values_de = []
    all_files_values_pso = []
    all_files_values_gwo = []

    xls = pd.ExcelFile(
        '/home/xerk/Escritorio/Códigos para artículo/Artículo PSO ED GW/Analisis de Convergencia Chen.xlsx', engine='openpyxl')
    nombres_hojas = xls.sheet_names

    for hoja in nombres_hojas:
        # Cargar la hoja actual en un DataFrame
        df = pd.read_excel('/home/xerk/Escritorio/Códigos para artículo/Artículo PSO ED GW/Analisis de Convergencia Chen.xlsx',
                        sheet_name=hoja, engine='openpyxl')

        # Extraer los datos de la columna 'Unnamed: 16'
        datos_columna = df['Unnamed: 16'].tolist()

        # Verificar si la columna tiene al menos un dato para evitar errores
        if datos_columna:
            # Agregar solo el último valor de la columna a arreglo_final
            all_files_values_gwo.append(datos_columna[-1])
        else:
            # Opcional: Agregar un valor por defecto o manejar de otra manera si la columna está vacía
            # O cualquier valor que consideres apropiado
            all_files_values_gwo.append(None)

    # Nota: No es necesario transponer el arreglo en este caso, ya que solo se está agregando un valor por hoja

    # Continúa con el procesamiento de los datos según necesites

    for i in range(10):
        # Nombre del archivo basado en el número de iteración
        filename_de = f'progress_de_chen_polaco{i}.txt'
        filename_pso = f'progress_pso_chen_polaco{i}.txt'

        # Lista para almacenar los valores de la última columna de este archivo
        last_column_values_de = []
        last_column_values_pso = []

        # Abrir el archivo de texto y leerlo línea por línea para DE
        with open(filename_de, 'r') as file:
            for line in file:
                # Ignorar líneas que no contienen datos
                if line.startswith("=") or "n_gen" in line or "Elapsed time" in line:
                    continue

                # Dividir la línea en sus componentes individuales
                parts = line.split('|')

                # Tomar el último valor, eliminar espacios y convertirlo a un número flotante
                value = float(parts[-1].strip())

                # Agregar el valor a la lista
                last_column_values_de.append(value)

        # Abrir el archivo de texto y leerlo línea por línea para PSO
        with open(filename_pso, 'r') as file:
            for line in file:
                # Ignorar líneas que no contienen datos
                if line.startswith("=") or "n_gen" in line or "Elapsed time" in line:
                    continue

                # Dividir la línea en sus componentes individuales
                parts = line.split('|')

                # Tomar el último valor, eliminar espacios y convertirlo a un número flotante
                value = float(parts[-1].strip())

                # Agregar el valor a la lista
                last_column_values_pso.append(value)

        # Agregar solo el último valor de la última columna de este archivo a la lista principal
        if last_column_values_de:
            all_files_values_de.append(last_column_values_de[-1])
        if last_column_values_pso:
            all_files_values_pso.append(last_column_values_pso[-1])

    # Ahora, all_files_values_de y all_files_values_pso contienen los últimos valores de la última columna de cada archivo para DE y PSO, respectivamente.
    # print(all_files_values_de)
    # print(all_files_values_pso)
    # print(all_files_values_gwo)


    all_files_array_de = -np.array(all_files_values_de)
    all_files_array_pso = -np.array(all_files_values_pso)

    stat_t = np.zeros(10)
    p_t = np.zeros(10)
    stat_w = np.zeros(10)
    p_w = np.zeros(10)

    stat_t, p_t = ttest_ind(all_files_array_pso, all_files_values_gwo)
    stat_w, p_w = wilcoxon(all_files_array_pso, all_files_values_gwo)

    print(
        f"Resultados Test t: Estadístico promedio = {stat_t}, P-valor promedio = {p_t}")
    print(
        f"Resultados Test de Wilcoxon: Estadístico promedio = {stat_w}, P-valor promedio = {p_w}")

    # Interpretación del test t
    alpha = 0.05
    if p_t > alpha:
        print('Test t: No hay evidencia suficiente para rechazar la hipótesis nula (H0), por lo que se asume que las medias son iguales.')
    else:
        print('Test t: Se rechaza la hipótesis nula (H0), lo que indica que las medias son diferentes.')

    # Interpretación del test de Wilcoxon
    if p_w > alpha:
        print('Test de Wilcoxon: No hay evidencia suficiente para rechazar la hipótesis nula (H0), por lo que se asume que las medianas son iguales.')
    else:
        print('Test de Wilcoxon: Se rechaza la hipótesis nula (H0), lo que indica que las medianas son diferentes.')

    mean_de = np.mean(all_files_array_de)
    mean_gwo = np.mean(all_files_values_gwo)
    mean_pso = np.mean(all_files_array_pso)
    
    print(mean_de, mean_pso, mean_gwo)

    if mean_pso > mean_gwo:
        print("PSO tiene una media mayor que GWO.")
    elif mean_pso < mean_gwo:
        print("GWO tiene una media mayor que PSO.")
    else:
        print("Ambos grupos tienen la misma media.")

    return stat_t, p_t, stat_w, p_w

if __name__ == '__main__':
    main()
