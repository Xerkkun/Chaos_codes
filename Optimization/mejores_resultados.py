import numpy as np
import pandas as pd


def get_top_three_indices(arr):
    indices = np.argsort(arr)[-3:][::-1]
    return indices, arr[indices]


def save_top_three_to_file(filename, data, indices):
    with open(filename, 'w') as file:
        for idx in indices:
            file.write(' '.join(map(str, data[idx])) + '\n')


def main():
    # Lista para almacenar los valores de la última columna de cada archivo
    all_files_values_de = []
    all_files_values_pso = []
    all_files_values_gwo = []

    xls = pd.ExcelFile(
        '/home/xerk/Escritorio/Códigos para artículo/Artículo PSO ED GW/Analisis de Convergencia Rossler.xlsx', engine='openpyxl')
    nombres_hojas = xls.sheet_names

    arreglo_final = []

    for hoja in nombres_hojas:
        df = pd.read_excel(
            '/home/xerk/Escritorio/Códigos para artículo/Artículo PSO ED GW/Analisis de Convergencia Rossler.xlsx', sheet_name=hoja, engine='openpyxl')
        datos_columna = df['Unnamed: 16'].tolist()
        arreglo_final.append(datos_columna)

    # Transponer el arreglo para que cada hoja contribuya a una columna diferente
    arreglo_transpuesto = list(map(list, zip(*arreglo_final)))

    for i in range(10):
        filename_de = f'/home/xerk/Escritorio/Códigos para artículo/Artículo PSO ED GW/progress_de_rossler_polaco{i}.txt'
        filename_pso = f'/home/xerk/Escritorio/Códigos para artículo/Artículo PSO ED GW/progress_pso_rossler_polaco{i}.txt'

        last_column_values_de = []
        last_column_values_pso = []

        with open(filename_de, 'r') as file:
            for line in file:
                if line.startswith("=") or "n_gen" in line or "Elapsed time" in line:
                    continue
                parts = line.split('|')
                try:
                    value = abs(float(parts[-1].strip()))
                    last_column_values_de.append(value)
                except ValueError as e:
                    print(
                        f"Error processing line in DE file {i}: {line} - {e}")

        with open(filename_pso, 'r') as file:
            for line in file:
                if line.startswith("=") or "n_gen" in line or "Elapsed time" in line:
                    continue
                parts = line.split('|')
                try:
                    value = abs(float(parts[-1].strip()))
                    last_column_values_pso.append(value)
                except ValueError as e:
                    print(
                        f"Error processing line in PSO file {i}: {line} - {e}")

        all_files_values_de.append(last_column_values_de)
        all_files_values_pso.append(last_column_values_pso)

    all_files_array_de = np.array(all_files_values_de)
    all_files_array_pso = np.array(all_files_values_pso)
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

    print(f"DE - Fila: {row_de}, Valor máximo: {max_de}")
    print(f"PSO - Fila: {row_pso}, Valor máximo: {max_pso}")
    print(f"GWO - Hoja: {nombres_hojas[row_gwo]}, Valor máximo: {max_gwo}")

    # Procesar los archivos last_gen_de y last_gen_pso para obtener los valores más altos
    with open(f'/home/xerk/Escritorio/Códigos para artículo/Artículo PSO ED GW/last_gen_de_rossler_polaco{row_de}.dat', 'r') as file:
        de_data = [line.split() for line in file.readlines()]
        de_values = np.array([abs(float(row[-1])) for row in de_data])
        top_three_de_indices, top_three_de_values = get_top_three_indices(
            de_values)

    with open(f'/home/xerk/Escritorio/Códigos para artículo/Artículo PSO ED GW/last_gen_pso_rossler_polaco{row_pso}.dat', 'r') as file:
        pso_data = [line.split() for line in file.readlines()]
        pso_values = np.array([abs(float(row[-1])) for row in pso_data])
        top_three_pso_indices, top_three_pso_values = get_top_three_indices(
            pso_values)

    save_top_three_to_file('/home/xerk/Escritorio/Códigos para artículo/Artículo PSO ED GW/top_3_de_rossler.txt',
                           de_data, top_three_de_indices)
    save_top_three_to_file('/home/xerk/Escritorio/Códigos para artículo/Artículo PSO ED GW/top_3_pso_rossler.txt',
                           pso_data, top_three_pso_indices)

    # Obtener los tres valores más altos en la columna K (columna 10) y sus generaciones en la columna G (columna 6)
    hoja_gwo = nombres_hojas[row_gwo]
    df_gwo = pd.read_excel(
        '/home/xerk/Escritorio/Códigos para artículo/Artículo PSO ED GW/Analisis de Convergencia Rossler.xlsx', sheet_name=hoja_gwo, engine='openpyxl')
    top_3_values_k = df_gwo.nlargest(3, 'Unnamed: 10')['Unnamed: 10'].values
    top_3_generations = df_gwo.loc[df_gwo['Unnamed: 10'].isin(
        top_3_values_k), 'Unnamed: 6'].values

    # Guardar los valores de las columnas H, I, J para los tres valores más altos de GWO
    with open('/home/xerk/Escritorio/Códigos para artículo/Artículo PSO ED GW/top_3_gwo_rossler.txt', 'w') as file:
        for value, generation in zip(top_3_values_k, top_3_generations):
            row = df_gwo.loc[df_gwo['Unnamed: 10'] == value].index[0]
            file.write(
                f"{df_gwo.loc[row, 'Unnamed: 7']} {df_gwo.loc[row, 'Unnamed: 8']} {df_gwo.loc[row, 'Unnamed: 9']} {df_gwo.loc[row, 'Unnamed: 10']}\n")

    print(f"Top 3 DE values: {top_three_de_values}")
    print(f"Top 3 PSO values: {top_three_pso_values}")
    print(f"Top 3 GWO values: {top_3_values_k}")


if __name__ == '__main__':
    main()
