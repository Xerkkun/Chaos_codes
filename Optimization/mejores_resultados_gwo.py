import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def main():
    # Lista para almacenar los valores de la última columna de cada archivo
    all_files_values_de = []
    all_files_values_pso = []
    all_files_values_gwo = []

    xls = pd.ExcelFile(
        '/home/xerk/Escritorio/Códigos para artículo/Artículo PSO ED GW/Analisis de Convergencia Lorenz - WOLF.xlsx', engine='openpyxl')
    nombres_hojas = xls.sheet_names

    arreglo_final = []

    for hoja in nombres_hojas:
        df = pd.read_excel(
            '/home/xerk/Escritorio/Códigos para artículo/Artículo PSO ED GW/Analisis de Convergencia Lorenz - WOLF.xlsx', sheet_name=hoja, engine='openpyxl')
        datos_columna = df['Unnamed: 16'].tolist()
        arreglo_final.append(datos_columna)

    # Transponer el arreglo para que cada hoja contribuya a una columna diferente
    arreglo_transpuesto = list(map(list, zip(*arreglo_final)))

    lista_datos_de = []
    lista_datos_pso = []

    for i in range(10):
        filename_de = f'progress_de_polaco{i}.txt'
        filename_pso = f'progress_pso_polaco{i}.txt'

        last_column_values_de = []
        last_column_values_pso = []

        with open(filename_de, 'r') as file:
            for line in file:
                if line.startswith("=") or "n_gen" in line or "Elapsed time" in line:
                    continue
                parts = line.split('|')
                value = float(parts[-1].strip())
                last_column_values_de.append(value)

        with open(filename_pso, 'r') as file:
            for line in file:
                if line.startswith("=") or "n_gen" in line or "Elapsed time" in line:
                    continue
                parts = line.split('|')
                value = float(parts[-1].strip())
                last_column_values_pso.append(value)

        all_files_values_de.append(last_column_values_de)
        all_files_values_pso.append(last_column_values_pso)

    all_files_array_de = -np.array(all_files_values_de)
    all_files_array_pso = -np.array(all_files_values_pso)
    all_files_array_gwo_trans = np.array(arreglo_transpuesto)
    all_files_array_gwo = all_files_array_gwo_trans.T

    index_de = np.argmax(all_files_array_de)
    index_pso = np.argmax(all_files_array_pso)
    index_gwo = np.argmax(all_files_array_gwo)

    row_de, col_de = np.unravel_index(index_de, all_files_array_de.shape)
    row_pso, col_pso = np.unravel_index(index_pso, all_files_array_pso.shape)
    row_gwo, col_gwo = np.unravel_index(index_gwo, all_files_array_gwo.shape)

    max_de = all_files_array_de[row_de, col_de]
    max_pso = all_files_array_pso[row_pso, col_pso]
    max_gwo = all_files_array_gwo[row_gwo, col_gwo]

    print(row_de, max_de)
    print(row_pso, max_pso)
    print(row_gwo, max_gwo)

    # Obtener los tres valores más altos en la columna K (columna 10) y sus generaciones en la columna G (columna 6)
    hoja_gwo = nombres_hojas[row_gwo]
    df_gwo = pd.read_excel(
        '/home/xerk/Escritorio/Códigos para artículo/Artículo PSO ED GW/Analisis de Convergencia Lorenz - WOLF.xlsx', sheet_name=hoja_gwo, engine='openpyxl')
    top_3_values_k = df_gwo.nlargest(3, 'Unnamed: 10')['Unnamed: 10'].values
    top_3_generations = df_gwo.loc[df_gwo['Unnamed: 10'].isin(
        top_3_values_k), 'Unnamed: 6'].values

    for i, (value, generation) in enumerate(zip(top_3_values_k, top_3_generations)):
        print(
            f"Top {i+1} value in column K: {value}, corresponding generation in column G: {generation}")


if __name__ == '__main__':
    main()
