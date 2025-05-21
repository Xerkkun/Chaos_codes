import numpy as np
import pandas as pd
from scipy.stats import wilcoxon, levene

# Cargar el archivo ODS
file_path = '/home/xerk/Escritorio/Soluciones.ods'
xls = pd.read_excel(file_path, sheet_name='Hoja1', engine='odf')

# Limpiar los datos y extraer las columnas relevantes
# Renombrar las columnas para facilitar el acceso
xls.columns = ['Ejecución', 'Lorenz_DE', 'Lorenz_PSO', 'Lorenz_GWO', 'NaN1',
               'Rossler_DE', 'Rossler_PSO', 'Rossler_GWO', 'NaN2', 'Chen_DE', 'Chen_PSO', 'Chen_GWO']

# Eliminar las columnas innecesarias
xls = xls[['Lorenz_DE', 'Lorenz_PSO', 'Lorenz_GWO', 'Rossler_DE',
           'Rossler_PSO', 'Rossler_GWO', 'Chen_DE', 'Chen_PSO', 'Chen_GWO']]

# Eliminar filas con valores no numéricos
xls = xls.apply(pd.to_numeric, errors='coerce')
xls = xls.dropna()

# Convertir los datos a float
xls = xls.astype(float)

# Función para realizar la prueba de Wilcoxon


def prueba_wilcoxon(v1, v2):
    s, p = wilcoxon(v1, v2, zero_method='wilcox',
                    correction=True, alternative='two-sided')
    if p >= 0.05:
        resultado = 'Los resultados son iguales'
    else:
        s, p = wilcoxon(v1, v2, zero_method='wilcox',
                        correction=True, alternative='less')
        if p >= 0.05:
            resultado = 'El segundo algoritmo es mejor'
        else:
            s, p = wilcoxon(v1, v2, zero_method='wilcox',
                            correction=True, alternative='greater')
            if p >= 0.05:
                resultado = 'El primer algoritmo es mejor'
            else:
                resultado = 'No se pudo determinar un ganador claro'
    return resultado

# Función para realizar la prueba de Levene


def prueba_levene(v1, v2):
    s, p = levene(v1, v2)
    if p >= 0.05:
        resultado = 'Las varianzas son iguales'
    else:
        resultado = 'Las varianzas son diferentes'
    return resultado


# Aplicar las pruebas para cada sistema caótico
sistemas = ['Lorenz', 'Rossler', 'Chen']
resultados = {}

for sistema in sistemas:
    de = xls[f'{sistema}_DE']
    pso = xls[f'{sistema}_PSO']
    gwo = xls[f'{sistema}_GWO']

    # Imprimir los vectores de datos
    print(f'Datos para el sistema {sistema}:')
    print(f'DE: {de.values}')
    print(f'PSO: {pso.values}')
    print(f'GWO: {gwo.values}')
    print('\n')

    resultado_de_vs_pso = prueba_wilcoxon(de, pso)
    resultado_de_vs_gwo = prueba_wilcoxon(de, gwo)
    resultado_pso_vs_gwo = prueba_wilcoxon(pso, gwo)

    varianza_de_vs_pso = prueba_levene(de, pso)
    varianza_de_vs_gwo = prueba_levene(de, gwo)
    varianza_pso_vs_gwo = prueba_levene(pso, gwo)

    resultados[sistema] = {
        'DE vs PSO': (resultado_de_vs_pso, varianza_de_vs_pso),
        'DE vs GWO': (resultado_de_vs_gwo, varianza_de_vs_gwo),
        'PSO vs GWO': (resultado_pso_vs_gwo, varianza_pso_vs_gwo)
    }

# Mostrar los resultados
for sistema, resultado in resultados.items():
    print(f'Resultados para el sistema {sistema}:')
    for comparacion, res in resultado.items():
        print(f'{comparacion}:')
        print(f'  Wilcoxon: {res[0]}')
        print(f'  Levene: {res[1]}')
    print('\n')
