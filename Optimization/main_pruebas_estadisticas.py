import numpy as np
from scipy.stats import wilcoxon, ttest_ind
import pandas as pd


def main():
    lista_datos_de = []
    lista_datos_pso = []

    for i in range(10):
        nombres_archivos_de = f"last_gen_de_polaco{i}.dat"
        nombres_archivos_pso = f"last_gen_pso_polaco{i}.dat"
        df_de = pd.read_csv(nombres_archivos_de, sep='\s+', header=None)
        df_pso = pd.read_csv(nombres_archivos_pso, sep='\s+', header=None)
        datos_columna_de = df_de.iloc[:, 3].tolist()
        datos_columna_pso = df_pso.iloc[:, 3].tolist()
        lista_datos_de.append(datos_columna_de)
        lista_datos_pso.append(datos_columna_pso)

    all_files_array_de = -np.array(lista_datos_de)
    all_files_array_pso = -np.array(lista_datos_pso)

    stat_t = np.zeros(10)
    p_t = np.zeros(10)
    stat_w = np.zeros(10)
    p_w = np.zeros(10)

    for j in range(9):
        stat_t[j], p_t[j] = ttest_ind(
            all_files_array_pso[j, :], all_files_array_de[j, :])
        stat_w[j], p_w[j] = wilcoxon(
            all_files_array_pso[j, :], all_files_array_de[j, :])

    mean_stat_t = np.mean(stat_t)
    mean_p_t = np.mean(p_t)
    mean_stat_w = np.mean(stat_w)
    mean_p_w = np.mean(p_w)

    print(
        f"Resultados Test t: Estadístico promedio = {mean_stat_t}, P-valor promedio = {mean_p_t}")
    print(
        f"Resultados Test de Wilcoxon: Estadístico promedio = {mean_stat_w}, P-valor promedio = {mean_p_w}")

    # Interpretación del test t
    alpha = 0.05
    if mean_p_t > alpha:
        print('Test t: No hay evidencia suficiente para rechazar la hipótesis nula (H0), por lo que se asume que las medias son iguales.')
    else:
        print('Test t: Se rechaza la hipótesis nula (H0), lo que indica que las medias son diferentes.')

    # Interpretación del test de Wilcoxon
    if mean_p_w > alpha:
        print('Test de Wilcoxon: No hay evidencia suficiente para rechazar la hipótesis nula (H0), por lo que se asume que las medianas son iguales.')
    else:
        print('Test de Wilcoxon: Se rechaza la hipótesis nula (H0), lo que indica que las medianas son diferentes.')

    mean_de = np.mean(all_files_array_de)
    mean_pso = np.mean(all_files_array_pso)
    print(mean_de, mean_pso)

    if mean_de > mean_pso:
        print("DE tiene una media mayor que PSO.")
    elif mean_de < mean_pso:
        print("PSO tiene una media mayor que el DE.")
    else:
        print("Ambos grupos tienen la misma media.")

    return mean_stat_t, mean_p_t, mean_stat_w, mean_p_w


if __name__ == '__main__':
    main()
